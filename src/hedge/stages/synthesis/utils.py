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

from hedge.configs.logger import logger

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
    """
    base_folder = process_path(folder_to_save)
    
    candidates = [base_folder + 'StructFilters/passStructFiltersSMILES.csv',
                  base_folder + 'Descriptors/passDescriptorsSMILES.csv',
                  base_folder + 'sampledMols.csv',
                 ]
            
    for candidate in candidates:
        if os.path.exists(candidate):
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
    """Load RAscore calculator - uses MolScore implementation."""
    global _rascore_available
    if _rascore_available is None:
        rascore_config = Path(__file__).parent / "RAScore.json"
        molscore_path = Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
        
        if (molscore_path / "molscore").exists() and rascore_config.exists():
            _rascore_available = True
            logger.debug("MolScore implementation available for RA score calculation")
        else:
            _rascore_available = False
            if not rascore_config.exists():
                logger.warning(f"RAScore config not found at {rascore_config}. RA scores will be set to np.nan")
            if not (molscore_path / "molscore").exists():
                logger.warning(f"MolScore module not found at {molscore_path}. RA scores will be set to np.nan")
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


def _calculate_ra_score(smiles, output_dir=None):
    """Calculate Retrosynthetic Accessibility score.
    
    Uses MolScore implementation.
    Returns probability (0-1) that molecule is synthesizable.
    Returns np.nan if calculation fails or module not available.
    
    Args:
        smiles: SMILES string
        
    Returns:
        RA score (0-1, higher is better) or np.nan
    """
    if not _load_rascore():
        return np.nan
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
    
        molscore_path = Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
        if str(molscore_path) not in sys.path:
            sys.path.insert(0, str(molscore_path))
        
        rascore_config_template = Path(__file__).parent / "RAScore.json"
        
        if not rascore_config_template.exists():
            logger.warning(f"RAScore config not found at {rascore_config_template}")
            return np.nan
        
        try:
            import contextlib
            with open(rascore_config_template, 'r') as f:
                rascore_config_data = json.load(f)
            
            if not output_dir:
                logger.warning('No output_dir provided for RA score, using default location')
                default_output = Path(__file__).parent.parent.parent.parent.parent / 'results' / 'MolScore'
                output_dir = str(default_output)
            
            rascore_config_data['output_dir'] = str(output_dir)
            
            output_dir_path = Path(output_dir)
            output_dir_path.mkdir(parents=True, exist_ok=True)
            rascore_config = output_dir_path / 'RAScore.json'
            with open(rascore_config, 'w') as f:
                json.dump(rascore_config_data, f, indent=2)
            
            molscore_logger = logging.getLogger("molscore")
            old_level = molscore_logger.level
            molscore_logger.setLevel(logging.ERROR)
            
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                from molscore.manager import MolScore
                smiles_list = [smiles]
                ms = MolScore(model_name="hedge",
                              task_config=str(rascore_config),
                              budget=len(smiles_list)
                              )
            molscore_logger.setLevel(old_level)
            scores = None
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                while not ms.finished:
                    scores = ms.score(smiles_list)
            if scores is not None:
                if hasattr(scores, '__len__') and len(scores) > 0:
                    score_value = scores[0]
                    if isinstance(score_value, dict):
                        score_value = score_value.get('RAScore_pred_proba') or score_value.get('score') or score_value.get('total_score')
                    try:
                        return float(score_value)
                    except (ValueError, TypeError):
                        logger.warning(f"MolScore returned non-numeric score for {smiles}: {score_value}")
                        return np.nan
                else:
                    logger.warning(f"MolScore returned empty scores for {smiles}")
                    return np.nan
            else:
                logger.warning(f"MolScore returned None for {smiles}")
                return np.nan
        except Exception as molscore_error:
            logger.warning(f"Failed to calculate RA score with MolScore for {smiles}: {molscore_error}")
            import traceback
            logger.debug(traceback.format_exc())
            return np.nan
    except Exception as e:
        logger.warning(f"Failed to calculate RA score for {smiles}: {e}")
        import traceback
        logger.warning(traceback.format_exc())
        return np.nan


def calculate_synthesis_scores(df, folder_to_save=None):
    """Calculate SA, SYBA, and RA scores for all molecules in DataFrame.
    
    Args:
        df: DataFrame with 'smiles' column
        
    Returns:
        DataFrame with added columns: sa_score, syba_score, ra_score
    """
    logger.info("Calculating synthetic accessibility scores...")
    _load_sascorer()
    _load_syba_model()
    _load_rascore()
    
    result_df = df.copy()
    ra_output_dir = None
    if folder_to_save:
        from pathlib import Path
        ra_output_dir = Path(folder_to_save) / 'Synthesis' / 'MolScore'
        ra_output_dir.mkdir(parents=True, exist_ok=True)
    
    result_df['sa_score'] = result_df['smiles'].apply(_calculate_sa_score)
    result_df['syba_score'] = result_df['smiles'].apply(_calculate_syba_score)
    result_df['ra_score'] = result_df['smiles'].apply(lambda smi: _calculate_ra_score(smi, ra_output_dir))
    
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
