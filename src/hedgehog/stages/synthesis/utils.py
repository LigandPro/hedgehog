import io
import os
import sys 
import json
import logging
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any

import pandas as pd
import numpy as np
from rdkit import Chem

from hedgehog.configs.logger import logger

_sascorer = None
_syba_model = None
_rascore_available = None

def process_path(folder_to_save, key_word=None):
    """Ensure path ends with '/' and create directory if needed."""
    if not folder_to_save.endswith('/'):
        folder_to_save += '/'
    
    if key_word:
        folder_to_save += f'{key_word}/'
    
    os.makedirs(folder_to_save, exist_ok=True)
    return folder_to_save


def prepare_input_smiles(input_df, output_file):
    """Prepare input SMILES file for aizynthfinder.
    
    Args:
        input_df: DataFrame with 'smiles' column
        output_file: Path to save SMILES file
        
    Returns:
        Number of molecules written
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    smiles_list = input_df['smiles'].dropna().tolist()
    
    with open(output_file, 'w') as f:
        for smi in smiles_list:
            f.write(f"{smi}\n")
    
    return len(smiles_list)


def run_aizynthfinder(input_smiles_file, output_json_file, config_file):
    """Run AiZynthFinder CLI to perform retrosynthesis analysis.
    
    Args:
        input_smiles_file: Path to input SMILES file
        output_json_file: Path to output JSON file
        config_file: Path to AiZynthFinder config file
        
    Returns:
        True if successful, False otherwise
    """
    output_json_file.parent.mkdir(parents=True, exist_ok=True)
    
    aizynthfinder_dir = config_file.parent.parent
    input_abs = input_smiles_file.resolve()
    output_abs = output_json_file.resolve()
    config_abs = config_file.resolve()
    
    cmd = f"cd {aizynthfinder_dir} && uv run aizynthcli --config {config_abs} --smiles {input_abs} --output {output_abs}"
    
    try:
        logger.info("Running retrosynthesis analysis...")
        logger.debug(f"Command: {cmd}")
        
        subprocess.run(cmd,
                       shell=True,
                       capture_output=True,
                       text=True,
                       check=True
                      )
        
        logger.info("Retrosynthesis analysis completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if e.stderr else e.stdout
        logger.error(f"Retrosynthesis analysis failed with exit code {e.returncode}")
        if error_msg:
            logger.error(f"Error output: {error_msg[:1000]}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error running retrosynthesis analysis: {e}")
        return False


def parse_retrosynthesis_results(json_file):
    """Parse retrosynthesis JSON results into a DataFrame.
    
    Args:
        json_file: Path to JSON output from aizynthfinder
        
    Returns:
        DataFrame with columns: index, SMILES, solved, search_time
    """
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        if 'data' not in data:
            logger.warning(f"No 'data' key found in JSON file {json_file}")
            return pd.DataFrame(columns=['index', 'SMILES', 'solved', 'search_time'])
        
        results = []
        for item in data['data']:
            results.append({'index': item.get('index', -1),
                            'SMILES': item.get('target', ''),
                            'solved': 1 if item.get('is_solved', False) else 0,
                            'search_time': item.get('search_time', 0.0)
                          })
        
        return pd.DataFrame(results)
    except Exception as e:
        logger.error(f"Error parsing retrosynthesis results: {e}")
        return pd.DataFrame(columns=['index', 'SMILES', 'solved', 'search_time'])


def merge_retrosynthesis_results(input_df, retrosynth_df):
    """Merge retrosynthesis results with input DataFrame.

    Args:
        input_df: Original input DataFrame with molecules (may have duplicate SMILES)
        retrosynth_df: DataFrame with retrosynthesis results indexed by position
        
    Returns:
        Merged DataFrame with retrosynthesis information, preserving all input rows
    """
    merged = input_df.reset_index(drop=True).copy()
    retrosynth_df_copy = retrosynth_df.reset_index(drop=True)
    
    merged['solved'] = 0
    merged['search_time'] = 0.0
    
    for idx, row in retrosynth_df_copy.iterrows():
        if idx < len(merged):
            merged.loc[idx, 'solved'] = row.get('solved', 0)
            merged.loc[idx, 'search_time'] = row.get('search_time', 0.0)
    
    return merged


def get_input_path(config: Dict[str, Any], folder_to_save: str) -> str:
    """Determine input path for synthesis stage.

    Checks for StructFilters output first, then falls back to other sources.
    Supports both new hierarchical structure and legacy flat structure.
    """
    base_folder = process_path(folder_to_save)

    # New structure paths
    candidates = [
        os.path.join(base_folder, 'stages', '03_structural_filters_post', 'filtered_molecules.csv'),
        os.path.join(base_folder, 'stages', '01_descriptors_initial', 'filtered', 'filtered_molecules.csv'),
        os.path.join(base_folder, 'input', 'sampled_molecules.csv'),
        # Legacy structure paths
        base_folder + 'StructFilters/passStructFiltersSMILES.csv',
        base_folder + 'Descriptors/passDescriptorsSMILES.csv',
        base_folder + 'sampledMols.csv',
    ]

    for candidate in candidates:
        if os.path.exists(candidate):
            logger.debug(f"Using input file: {candidate}")
            return candidate

    logger.warning("No processed data found, using molecules from config")
    return config.get('generated_mols_path', '')


def _load_sascorer():
    global _sascorer
    if _sascorer is None:
        try:
            from rdkit.Chem import RDConfig
            import sys
            sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
            import sascorer
            _sascorer = sascorer
        except Exception as e:
            logger.warning(f"Failed to load SA Score module: {e}")
            _sascorer = False
    return _sascorer


def _load_syba_model():
    """Load SYBA model"""
    global _syba_model
    if _syba_model is None:
        logger.debug("Loading SYBA model...")
        try:
            from syba import syba
            _syba_model = syba.SybaClassifier()
            _syba_model.fitDefaultScore()
            logger.debug("SYBA model loaded successfully")
        except Exception as e:
            _syba_model = False
            logger.error(f"Failed to load SYBA model: {e}")
            logger.warning("SYBA scores will be set to np.nan")
    return _syba_model


def _load_rascore():
    """Check if RAScore model and rascore-env are available."""
    global _rascore_available
    if _rascore_available is None:
        molscore_path = Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
        model_path = molscore_path / "molscore" / "data" / "models" / "RAScore" / "XGB_chembl_ecfp_counts" / "model.pkl"
        
        if model_path.exists():
            _rascore_available = True
            logger.debug("RAScore model available for calculation")
        else:
            _rascore_available = False
            logger.warning(f"RAScore model not found at {model_path}. RA scores will be set to np.nan")
    return _rascore_available


def _calculate_sa_score(smiles):
    """Calculate Synthetic Accessibility score (1=easy to 10=hard).
    
    Uses RDKit's SA Score calculator from contrib.
    Returns np.nan if calculation fails or module not available.
    
    Args:
        smiles: SMILES string
        
    Returns:
        SA score (1-10, lower is better) or np.nan
    """
    sascorer = _load_sascorer()
    if not sascorer:
        return np.nan
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
        return sascorer.calculateScore(mol)
    except Exception as e:
        logger.debug(f"Failed to calculate SA score for {smiles}: {e}")
        return np.nan


def _calculate_syba_score(smiles):
    """Calculate SYBA (SYnthetic Bayesian Accessibility) score.
    
    Uses SYBA model. Higher scores indicate more synthetically accessible.
    Returns np.nan if calculation fails or module not available.
    
    Args:
        smiles: SMILES string
        
    Returns:
        SYBA score or np.nan
    """
    syba_model = _load_syba_model()
    if not syba_model:
        return np.nan
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
        return syba_model.predict(mol=mol)
    except Exception as e:
        logger.debug(f"Failed to calculate SYBA score for {smiles}: {e}")
        return np.nan


def _calculate_ra_scores_batch(smiles_list, config=None):
    """Calculate RA scores for multiple molecules in a single batch.
    Args:
        smiles_list: List of SMILES strings
        config: Configuration dict (can contain 'rascore_conda_prefix' for manual override)
        
    Returns:
        List of RA scores (0-1, higher is better), with np.nan for failed calculations
    """
    if not _load_rascore():
        return [np.nan] * len(smiles_list)
    
    try:
        conda_prefix = None
        if config and isinstance(config, dict):
            conda_prefix = config.get('rascore_conda_prefix')
        
        if not conda_prefix:
            common_locations = [
                Path.home() / 'miniconda3' / 'envs' / 'rascore-env',
                Path.home() / 'anaconda3' / 'envs' / 'rascore-env',
                Path.home() / 'conda' / 'envs' / 'rascore-env',
                Path('/opt/conda/envs/rascore-env'),
            ]
            for location in common_locations:
                if location.exists():
                    conda_prefix = str(location)
                    break
        
        if not conda_prefix or not Path(conda_prefix).exists():
            logger.warning('Could not find rascore-env. RA scores will be set to np.nan')
            return [np.nan] * len(smiles_list)
        
        molscore_path = Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
        model_path = molscore_path / "molscore" / "data" / "models" / "RAScore" / "XGB_chembl_ecfp_counts" / "model.pkl"
        
        if not model_path.exists():
            logger.warning(f"RAScore model not found at {model_path}")
            return [np.nan] * len(smiles_list)
        
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.smi', delete=False) as f:
            smiles_file = f.name
            for smi in smiles_list:
                f.write(f"{smi}\n")
        
        script_content = f"""
import sys
import pickle as pkl
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def get_ecfp6_counts(smiles, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprint(mol, 3, useFeatures=False, useCounts=True)
    arr = np.zeros((nBits,), dtype=np.int32)
    for idx, v in fp.GetNonzeroElements().items():
        arr[idx % nBits] += v
    return arr

with open('{model_path}', 'rb') as f:
    clf = pkl.load(f)

with open('{smiles_file}', 'r') as f:
    smiles_list = [line.strip() for line in f]

for smiles in smiles_list:
    fp = get_ecfp6_counts(smiles)
    if fp is None:
        print('nan')
    else:
        try:
            prob = clf.predict_proba(np.array([fp]))[0, 1]
            print(prob)
        except Exception:
            print('nan')
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            temp_script = f.name
            f.write(script_content)
        
        try:
            python_exec = Path(conda_prefix) / 'bin' / 'python'
            if not python_exec.exists():
                logger.warning(f"Python not found in rascore-env: {python_exec}")
                return [np.nan] * len(smiles_list)
            
            result = subprocess.run(
                [str(python_exec), temp_script],
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if result.returncode == 0:
                scores = []
                for line in result.stdout.strip().split('\n'):
                    try:
                        score = float(line.strip())
                        scores.append(score if not np.isnan(score) else np.nan)
                    except ValueError:
                        scores.append(np.nan)
                
                if len(scores) != len(smiles_list):
                    logger.warning(f"Expected {len(smiles_list)} scores but got {len(scores)}")
                    return [np.nan] * len(smiles_list)
                return scores
            else:
                logger.debug(f"RA score batch calculation failed: {result.stderr}")
                return [np.nan] * len(smiles_list)
        finally:
            try:
                os.unlink(temp_script)
                os.unlink(smiles_file)
            except:
                pass
                
    except Exception as e:
        logger.debug(f"Failed to calculate RA scores in batch: {e}")
        return [np.nan] * len(smiles_list)


def calculate_synthesis_scores(df, folder_to_save=None, config=None):
    """Calculate SA, SYBA, and RA scores for all molecules in DataFrame.
    
    Args:
        df: DataFrame with 'smiles' column
        folder_to_save: Optional folder to save RA score outputs
        config: Optional config dict (can contain 'rascore_conda_prefix' for manual override)
        
    Returns:
        DataFrame with added columns: sa_score, syba_score, ra_score
    """
    logger.info("Calculating synthetic accessibility scores...")
    _load_sascorer()
    _load_syba_model()
    _load_rascore()
    
    result_df = df.copy()
    result_df['sa_score'] = result_df['smiles'].apply(_calculate_sa_score)
    result_df['syba_score'] = result_df['smiles'].apply(_calculate_syba_score)
    
    smiles_list = result_df['smiles'].tolist()
    ra_scores = _calculate_ra_scores_batch(smiles_list, config)
    result_df['ra_score'] = ra_scores
    
    for score_name in ['sa_score', 'syba_score', 'ra_score']:
        valid_scores = result_df[score_name].dropna()
        if len(valid_scores) > 0:
            logger.info(f"  {score_name}: calculated for {len(valid_scores)}/{len(df)} molecules "
                       f"(mean={valid_scores.mean():.2f}, std={valid_scores.std():.2f})")
        else:
            logger.debug(f"  {score_name}: could not be calculated (module not available)")
    
    return result_df


def apply_synthesis_score_filters(df, config):
    """Apply filters based on synthesis score thresholds.
    
    Each molecule is checked independently against each criterion (SA, RA, SYBA).
    A molecule must pass ALL filters for which it has valid scores.
    Molecules with NaN scores for a criterion are not filtered by that criterion.
    Only molecules that pass all applicable score filters are returned.
    
    Args:
        df: DataFrame with synthesis scores
        config: Configuration with threshold settings
        
    Returns:
        Filtered DataFrame with molecules that pass all applicable filters
    """
    filtered_df = df.copy()
    initial_count = len(filtered_df)
    pass_mask = pd.Series([True] * len(filtered_df), index=filtered_df.index)
    
    sa_min = config.get('sa_score_min', 0)
    sa_max = config.get('sa_score_max', 'inf')
    if 'sa_score' in filtered_df.columns:
        valid_scores = filtered_df['sa_score'].dropna()
        if len(valid_scores) > 0:
            if sa_max != 'inf':
                sa_mask = (filtered_df['sa_score'].isna()) | \
                         ((filtered_df['sa_score'] >= sa_min) & (filtered_df['sa_score'] <= sa_max))
            else:
                sa_mask = (filtered_df['sa_score'].isna()) | (filtered_df['sa_score'] >= sa_min)
            pass_mask = pass_mask & sa_mask
        else:
            logger.info(f"SA score filter: skipped (no valid scores calculated)")
    
    ra_min = config.get('ra_score_min', 0)
    ra_max = config.get('ra_score_max', 'inf')
    if 'ra_score' in filtered_df.columns:
        valid_scores = filtered_df['ra_score'].dropna()
        if len(valid_scores) > 0:
            if ra_max != 'inf':
                ra_mask = (filtered_df['ra_score'].isna()) | \
                         ((filtered_df['ra_score'] >= ra_min) & (filtered_df['ra_score'] <= ra_max))
            else:
                ra_mask = (filtered_df['ra_score'].isna()) | (filtered_df['ra_score'] >= ra_min)
            pass_mask = pass_mask & ra_mask
        else:
            logger.info(f"RA score filter: skipped (no valid scores calculated)")
    
    syba_min = config.get('syba_score_min', 0)
    syba_max = config.get('syba_score_max', 'inf')
    if 'syba_score' in filtered_df.columns:
        valid_scores = filtered_df['syba_score'].dropna()
        if len(valid_scores) > 0:
            if syba_max != 'inf':
                syba_mask = (filtered_df['syba_score'].isna()) | \
                           ((filtered_df['syba_score'] >= syba_min) & (filtered_df['syba_score'] <= syba_max))
            else:
                syba_mask = (filtered_df['syba_score'].isna()) | (filtered_df['syba_score'] >= syba_min)
            pass_mask = pass_mask & syba_mask
        else:
            logger.info(f"SYBA score filter: skipped (no valid scores calculated)")
    
    filtered_df = filtered_df[pass_mask]
    
    return filtered_df
