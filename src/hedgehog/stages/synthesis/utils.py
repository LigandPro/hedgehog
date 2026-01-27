import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import logger
from hedgehog.stages.structFilters.utils import process_path
from hedgehog.utils.input_paths import get_all_input_candidates

# Lazy-loaded module cache
_lazy_cache: dict[str, Any] = {
    "sascorer": None,
    "syba_model": None,
    "rascore_available": None,
}


def _get_cached(key: str, loader: callable) -> Any:
    """Generic lazy loading with caching.

    Args:
        key: Cache key name
        loader: Function that returns the loaded value or False on failure

    Returns:
        Cached value, or False if loading failed
    """
    if _lazy_cache[key] is None:
        _lazy_cache[key] = loader()
    return _lazy_cache[key]


def prepare_input_smiles(input_df, output_file):
    """Prepare input SMILES file for aizynthfinder.

    Args:
        input_df: DataFrame with 'smiles' column
        output_file: Path to save SMILES file

    Returns:
        Number of molecules written
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    smiles_list = input_df["smiles"].dropna().tolist()

    with open(output_file, "w") as f:
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
        logger.debug("Command: %s", cmd)

        subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

        logger.info("Retrosynthesis analysis completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if e.stderr else e.stdout
        logger.error("Retrosynthesis analysis failed with exit code %d", e.returncode)
        if error_msg:
            logger.error("Error output: %s", error_msg[:1000])
        return False
    except Exception as e:
        logger.error("Unexpected error running retrosynthesis analysis: %s", e)
        return False


def parse_retrosynthesis_results(json_file):
    """Parse retrosynthesis JSON results into a DataFrame.

    Args:
        json_file: Path to JSON output from aizynthfinder

    Returns:
        DataFrame with columns: index, SMILES, solved, search_time
    """
    try:
        with open(json_file) as f:
            data = json.load(f)

        if "data" not in data:
            logger.warning("No 'data' key found in JSON file %s", json_file)
            return pd.DataFrame(columns=["index", "SMILES", "solved", "search_time"])

        results = []
        for item in data["data"]:
            results.append(
                {
                    "index": item.get("index", -1),
                    "SMILES": item.get("target", ""),
                    "solved": 1 if item.get("is_solved", False) else 0,
                    "search_time": item.get("search_time", 0.0),
                }
            )

        return pd.DataFrame(results)
    except Exception as e:
        logger.error("Error parsing retrosynthesis results: %s", e)
        return pd.DataFrame(columns=["index", "SMILES", "solved", "search_time"])


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

    merged["solved"] = 0
    merged["search_time"] = 0.0

    for idx, row in retrosynth_df_copy.iterrows():
        if idx < len(merged):
            merged.loc[idx, "solved"] = row.get("solved", 0)
            merged.loc[idx, "search_time"] = row.get("search_time", 0.0)

    return merged


def _get_input_path_candidates(base_folder: str) -> list:
    """Get ordered list of candidate input paths to check."""
    return [str(p) for p in get_all_input_candidates(Path(base_folder))]


def get_input_path(config: dict[str, Any], folder_to_save: str) -> str:
    """Determine input path for synthesis stage.

    Checks for StructFilters output first, then falls back to other sources.
    Supports both new hierarchical structure and legacy flat structure.
    """
    base_folder = process_path(folder_to_save)

    for candidate in _get_input_path_candidates(base_folder):
        if Path(candidate).exists():
            logger.debug("Using input file: %s", candidate)
            return candidate

    logger.warning("No processed data found, using molecules from config")
    return config.get("generated_mols_path", "")


def _load_sascorer_impl():
    """Load SA Score module implementation."""
    try:
        import sys

        from rdkit.Chem import RDConfig

        sys.path.append(str(Path(RDConfig.RDContribDir) / "SA_Score"))
        import sascorer

        return sascorer
    except Exception as e:
        logger.warning("Failed to load SA Score module: %s", e)
        return False


def _load_syba_model_impl():
    """Load SYBA model implementation."""
    logger.debug("Loading SYBA model...")
    try:
        from syba import syba

        model = syba.SybaClassifier()
        model.fitDefaultScore()
        logger.debug("SYBA model loaded successfully")
        return model
    except Exception as e:
        logger.error("Failed to load SYBA model: %s", e)
        logger.warning("SYBA scores will be set to np.nan")
        return False


def _load_rascore_impl():
    """Check if RAScore model is available."""
    model_path = _get_rascore_model_path()
    if model_path.exists():
        logger.debug("RAScore model available for calculation")
        return True
    logger.warning(
        "RAScore model not found at %s. RA scores will be set to np.nan", model_path
    )
    return False


def _load_sascorer():
    """Load SA Score module with caching."""
    return _get_cached("sascorer", _load_sascorer_impl)


def _load_syba_model():
    """Load SYBA model with caching."""
    return _get_cached("syba_model", _load_syba_model_impl)


def _load_rascore():
    """Check if RAScore is available with caching."""
    return _get_cached("rascore_available", _load_rascore_impl)


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
        logger.debug("Failed to calculate SA score for %s: %s", smiles, e)
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
        logger.debug("Failed to calculate SYBA score for %s: %s", smiles, e)
        return np.nan


def _get_rascore_model_path() -> Path:
    """Get path to RAScore model file."""
    molscore_path = (
        Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
    )
    return (
        molscore_path
        / "molscore"
        / "data"
        / "models"
        / "RAScore"
        / "XGB_chembl_ecfp_counts"
        / "model.pkl"
    )


def _find_rascore_conda_env(config: dict[str, Any] | None) -> str | None:
    """Find rascore conda environment path.

    Args:
        config: Configuration dict that may contain 'rascore_conda_prefix'

    Returns:
        Path to conda environment or None if not found
    """
    if config and isinstance(config, dict):
        prefix = config.get("rascore_conda_prefix")
        if prefix and Path(prefix).exists():
            return prefix

    common_locations = [
        Path.home() / "miniconda3" / "envs" / "rascore-env",
        Path.home() / "anaconda3" / "envs" / "rascore-env",
        Path.home() / "conda" / "envs" / "rascore-env",
        Path("/opt/conda/envs/rascore-env"),
    ]
    for location in common_locations:
        if location.exists():
            return str(location)
    return None


_RASCORE_SCRIPT_TEMPLATE = """import pickle as pkl
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


def _parse_rascore_output(stdout: str, expected_count: int) -> list:
    """Parse RAScore subprocess output into list of scores.

    Args:
        stdout: Standard output from RAScore script
        expected_count: Expected number of scores

    Returns:
        List of scores (floats or np.nan)
    """
    scores = []
    for line in stdout.strip().split("\n"):
        try:
            score = float(line.strip())
            scores.append(score if not np.isnan(score) else np.nan)
        except ValueError:
            scores.append(np.nan)

    if len(scores) != expected_count:
        logger.warning("Expected %d scores but got %d", expected_count, len(scores))
        return [np.nan] * expected_count
    return scores


def _calculate_ra_scores_batch(
    smiles_list: list, config: dict[str, Any] | None = None
) -> list:
    """Calculate RA scores for multiple molecules in a single batch.

    Args:
        smiles_list: List of SMILES strings
        config: Configuration dict (can contain 'rascore_conda_prefix' for manual override)

    Returns:
        List of RA scores (0-1, higher is better), with np.nan for failed calculations
    """
    nan_list = [np.nan] * len(smiles_list)

    if not _load_rascore():
        return nan_list

    conda_prefix = _find_rascore_conda_env(config)
    if not conda_prefix:
        logger.warning("Could not find rascore-env. RA scores will be set to np.nan")
        return nan_list

    model_path = _get_rascore_model_path()
    if not model_path.exists():
        logger.warning("RAScore model not found at %s", model_path)
        return nan_list

    python_exec = Path(conda_prefix) / "bin" / "python"
    if not python_exec.exists():
        logger.warning("Python not found in rascore-env: %s", python_exec)
        return nan_list

    smiles_file = None
    temp_script = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            smiles_file = f.name
            for smi in smiles_list:
                f.write(f"{smi}\n")

        script_content = _RASCORE_SCRIPT_TEMPLATE.format(
            model_path=model_path, smiles_file=smiles_file
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
            temp_script = f.name
            f.write(script_content)

        result = subprocess.run(
            [str(python_exec), temp_script], capture_output=True, text=True, timeout=300
        )

        if result.returncode != 0:
            logger.debug("RA score batch calculation failed: %s", result.stderr)
            return nan_list

        return _parse_rascore_output(result.stdout, len(smiles_list))

    except Exception as e:
        logger.debug("Failed to calculate RA scores in batch: %s", e)
        return nan_list
    finally:
        for path in [temp_script, smiles_file]:
            if path:
                try:
                    Path(path).unlink()
                except OSError:
                    pass


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
    result_df["sa_score"] = result_df["smiles"].apply(_calculate_sa_score)
    result_df["syba_score"] = result_df["smiles"].apply(_calculate_syba_score)

    smiles_list = result_df["smiles"].tolist()
    ra_scores = _calculate_ra_scores_batch(smiles_list, config)
    result_df["ra_score"] = ra_scores

    for score_name in ["sa_score", "syba_score", "ra_score"]:
        valid_scores = result_df[score_name].dropna()
        if len(valid_scores) > 0:
            logger.info(
                "  %s: calculated for %d/%d molecules (mean=%.2f, std=%.2f)",
                score_name,
                len(valid_scores),
                len(df),
                valid_scores.mean(),
                valid_scores.std(),
            )
        else:
            logger.debug(
                "  %s: could not be calculated (module not available)", score_name
            )

    return result_df


def _build_score_filter_mask(
    df: pd.DataFrame, column: str, min_val: float, max_val
) -> pd.Series | None:
    """Build a boolean mask for filtering by score thresholds.

    Molecules with NaN scores pass the filter (not filtered by that criterion).

    Args:
        df: DataFrame containing the score column
        column: Name of the score column
        min_val: Minimum allowed value
        max_val: Maximum allowed value (can be 'inf' for no upper bound)

    Returns:
        Boolean Series mask, or None if column missing or no valid scores
    """
    if column not in df.columns:
        return None

    valid_scores = df[column].dropna()
    if len(valid_scores) == 0:
        logger.info("%s filter: skipped (no valid scores calculated)", column)
        return None

    is_na = df[column].isna()
    above_min = df[column] >= min_val

    if max_val == "inf":
        return is_na | above_min

    below_max = df[column] <= max_val
    return is_na | (above_min & below_max)


def apply_synthesis_score_filters(
    df: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
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
    score_filters = [
        ("sa_score", "sa_score_min", "sa_score_max"),
        ("ra_score", "ra_score_min", "ra_score_max"),
        ("syba_score", "syba_score_min", "syba_score_max"),
    ]

    pass_mask = pd.Series(True, index=df.index)

    for column, min_key, max_key in score_filters:
        min_val = config.get(min_key, 0)
        max_val = config.get(max_key, "inf")
        mask = _build_score_filter_mask(df, column, min_val, max_val)
        if mask is not None:
            pass_mask = pass_mask & mask

    return df[pass_mask].copy()
