import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem

from hedge.configs.logger import logger

_sascorer: Any = None
_syba_model: Any = None
_rascore_available: bool | None = None


def process_path(folder_to_save: str, key_word: str | None = None) -> str:
    """Ensure path ends with '/' and create directory if needed."""
    if not folder_to_save.endswith("/"):
        folder_to_save += "/"

    if key_word:
        folder_to_save += f"{key_word}/"

    Path(folder_to_save).mkdir(parents=True, exist_ok=True)
    return folder_to_save


def prepare_input_smiles(input_df: pd.DataFrame, output_file: Path) -> int:
    """Prepare input SMILES file for aizynthfinder.

    Args:
        input_df: DataFrame with 'smiles' column
        output_file: Path to save SMILES file

    Returns
    -------
        Number of molecules written
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    smiles_list = input_df["smiles"].dropna().tolist()

    with output_file.open("w") as f:
        f.writelines(f"{smi}\n" for smi in smiles_list)

    return len(smiles_list)


def run_aizynthfinder(
    input_smiles_file: Path, output_json_file: Path, config_file: Path
) -> bool:
    """Run AiZynthFinder CLI to perform retrosynthesis analysis.

    Args:
        input_smiles_file: Path to input SMILES file
        output_json_file: Path to output JSON file
        config_file: Path to AiZynthFinder config file

    Returns
    -------
        True if successful, False otherwise
    """
    output_json_file.parent.mkdir(parents=True, exist_ok=True)

    aizynthfinder_dir = config_file.parent.parent
    input_abs = input_smiles_file.resolve()
    output_abs = output_json_file.resolve()
    config_abs = config_file.resolve()

    cmd = (
        f"cd {aizynthfinder_dir} && uv run aizynthcli "
        f"--config {config_abs} --smiles {input_abs} --output {output_abs}"
    )

    logger.info("Running retrosynthesis analysis...")
    logger.debug(f"Command: {cmd}")

    try:
        subprocess.run(  # noqa: S602
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )

        logger.info("Retrosynthesis analysis completed successfully")
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if e.stderr else e.stdout
        logger.error(f"Retrosynthesis analysis failed with exit code {e.returncode}")
        if error_msg:
            logger.error(f"Error output: {error_msg[:1000]}")
        return False
    except (OSError, ValueError) as e:
        logger.error(f"Unexpected error running retrosynthesis analysis: {e}")
        return False
    else:
        return True


def parse_retrosynthesis_results(json_file: Path) -> pd.DataFrame:
    """Parse retrosynthesis JSON results into a DataFrame.

    Args:
        json_file: Path to JSON output from aizynthfinder

    Returns
    -------
        DataFrame with columns: index, SMILES, solved, search_time
    """
    empty_df = pd.DataFrame(columns=["index", "SMILES", "solved", "search_time"])
    try:
        with json_file.open() as f:
            data = json.load(f)

        if "data" not in data:
            logger.warning(f"No 'data' key found in JSON file {json_file}")
            return empty_df

        results = [
            {
                "index": item.get("index", -1),
                "SMILES": item.get("target", ""),
                "solved": 1 if item.get("is_solved", False) else 0,
                "search_time": item.get("search_time", 0.0),
            }
            for item in data["data"]
        ]

        return pd.DataFrame(results)
    except (OSError, ValueError, KeyError) as e:
        logger.error(f"Error parsing retrosynthesis results: {e}")
        return empty_df


def merge_retrosynthesis_results(
    input_df: pd.DataFrame, retrosynth_df: pd.DataFrame
) -> pd.DataFrame:
    """Merge retrosynthesis results with input DataFrame.

    Args:
        input_df: Original input DataFrame with molecules (may have duplicate SMILES)
        retrosynth_df: DataFrame with retrosynthesis results indexed by position

    Returns
    -------
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


def get_input_path(config: dict[str, Any], folder_to_save: str) -> str:
    """Determine input path for synthesis stage.

    Checks for StructFilters output first, then falls back to other sources.
    """
    base_folder = process_path(folder_to_save)

    candidates = [
        base_folder + "StructFilters/passStructFiltersSMILES.csv",
        base_folder + "Descriptors/passDescriptorsSMILES.csv",
        base_folder + "sampledMols.csv",
    ]

    for candidate in candidates:
        if Path(candidate).exists():
            return candidate

    logger.warning("No processed data found, using molecules from config")
    return config.get("generated_mols_path", "")


def _load_sascorer() -> Any:  # noqa: ANN401
    """Load SA Score module from RDKit contrib."""
    global _sascorer  # noqa: PLW0603
    if _sascorer is None:
        try:
            import sys  # noqa: PLC0415

            from rdkit.Chem import RDConfig  # noqa: PLC0415

            sys.path.append(str(Path(RDConfig.RDContribDir) / "SA_Score"))
            import sascorer  # noqa: PLC0415

            _sascorer = sascorer
        except (ImportError, AttributeError) as e:
            logger.warning(f"Failed to load SA Score module: {e}")
            _sascorer = False
    return _sascorer


def _load_syba_model() -> Any:  # noqa: ANN401
    """Load SYBA model."""
    global _syba_model  # noqa: PLW0603
    if _syba_model is None:
        logger.debug("Loading SYBA model...")
        try:
            from syba import syba  # noqa: PLC0415

            _syba_model = syba.SybaClassifier()
            _syba_model.fitDefaultScore()
            logger.debug("SYBA model loaded successfully")
        except (ImportError, AttributeError) as e:
            _syba_model = False
            logger.error(f"Failed to load SYBA model: {e}")
            logger.warning("SYBA scores will be set to np.nan")
    return _syba_model


def _load_rascore() -> bool:
    """Check if RAScore model and rascore-env are available."""
    global _rascore_available  # noqa: PLW0603
    if _rascore_available is None:
        molscore_path = (
            Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
        )
        model_path = (
            molscore_path
            / "molscore"
            / "data"
            / "models"
            / "RAScore"
            / "XGB_chembl_ecfp_counts"
            / "model.pkl"
        )

        if model_path.exists():
            _rascore_available = True
            logger.debug("RAScore model available for calculation")
        else:
            _rascore_available = False
            msg = (
                f"RAScore model not found at {model_path}. "
                "RA scores will be set to np.nan"
            )
            logger.warning(msg)
    return _rascore_available


def _calculate_sa_score(smiles: str) -> float:
    """Calculate Synthetic Accessibility score (1=easy to 10=hard).

    Uses RDKit's SA Score calculator from contrib.
    Returns np.nan if calculation fails or module not available.

    Args:
        smiles: SMILES string

    Returns
    -------
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
    except (ValueError, AttributeError) as e:
        logger.debug(f"Failed to calculate SA score for {smiles}: {e}")
        return np.nan


def _calculate_syba_score(smiles: str) -> float:
    """Calculate SYBA (SYnthetic Bayesian Accessibility) score.

    Uses SYBA model. Higher scores indicate more synthetically accessible.
    Returns np.nan if calculation fails or module not available.

    Args:
        smiles: SMILES string

    Returns
    -------
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
    except (ValueError, AttributeError) as e:
        logger.debug(f"Failed to calculate SYBA score for {smiles}: {e}")
        return np.nan


def _calculate_ra_scores_batch(
    smiles_list: list[str], config: dict[str, Any] | None = None
) -> list[float]:
    """Calculate RA scores for multiple molecules in a single batch.

    Args:
        smiles_list: List of SMILES strings
        config: Configuration dict (can contain 'rascore_conda_prefix' for
            manual override)

    Returns
    -------
        List of RA scores (0-1, higher is better), with np.nan for failed calculations
    """
    if not _load_rascore():
        return [np.nan] * len(smiles_list)

    conda_prefix = _find_conda_prefix(config)
    if not conda_prefix:
        logger.warning("Could not find rascore-env. RA scores will be set to np.nan")
        return [np.nan] * len(smiles_list)

    molscore_path = (
        Path(__file__).parent.parent.parent.parent.parent / "modules" / "MolScore"
    )
    model_path = (
        molscore_path
        / "molscore"
        / "data"
        / "models"
        / "RAScore"
        / "XGB_chembl_ecfp_counts"
        / "model.pkl"
    )

    if not model_path.exists():
        logger.warning(f"RAScore model not found at {model_path}")
        return [np.nan] * len(smiles_list)

    return _run_rascore_batch(smiles_list, conda_prefix, model_path)


def _find_conda_prefix(config: dict[str, Any] | None) -> str | None:
    """Find conda prefix for rascore-env."""
    conda_prefix = None
    if config and isinstance(config, dict):
        conda_prefix = config.get("rascore_conda_prefix")

    if not conda_prefix:
        common_locations = [
            Path.home() / "miniconda3" / "envs" / "rascore-env",
            Path.home() / "anaconda3" / "envs" / "rascore-env",
            Path.home() / "conda" / "envs" / "rascore-env",
            Path("/opt/conda/envs/rascore-env"),
        ]
        for location in common_locations:
            if location.exists():
                conda_prefix = str(location)
                break

    if conda_prefix and Path(conda_prefix).exists():
        return conda_prefix
    return None


def _run_rascore_batch(
    smiles_list: list[str], conda_prefix: str, model_path: Path
) -> list[float]:
    """Run RA score calculation in batch mode."""
    try:
        import tempfile  # noqa: PLC0415

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
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

        with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
            temp_script = f.name
            f.write(script_content)

        try:
            python_exec = Path(conda_prefix) / "bin" / "python"
            if not python_exec.exists():
                logger.warning(f"Python not found in rascore-env: {python_exec}")
                return [np.nan] * len(smiles_list)

            result = subprocess.run(  # noqa: S603
                [str(python_exec), temp_script],
                check=False,
                capture_output=True,
                text=True,
                timeout=300,
            )

            if result.returncode == 0:
                scores = [
                    float(line.strip()) if line.strip() != "nan" else np.nan
                    for line in result.stdout.strip().split("\n")
                ]

                if len(scores) != len(smiles_list):
                    msg = f"Expected {len(smiles_list)} scores but got {len(scores)}"
                    logger.warning(msg)
                    return [np.nan] * len(smiles_list)
                return scores
            logger.debug(f"RA score batch calculation failed: {result.stderr}")
            return [np.nan] * len(smiles_list)
        finally:
            try:
                Path(temp_script).unlink()
                Path(smiles_file).unlink()
            except OSError:
                pass

    except (OSError, ValueError) as e:
        logger.debug(f"Failed to calculate RA scores in batch: {e}")
        return [np.nan] * len(smiles_list)


def calculate_synthesis_scores(
    df: pd.DataFrame,
    folder_to_save: str | None = None,  # noqa: ARG001
    config: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Calculate SA, SYBA, and RA scores for all molecules in DataFrame.

    Args:
        df: DataFrame with 'smiles' column
        folder_to_save: Optional folder to save RA score outputs
        config: Optional config dict (can contain 'rascore_conda_prefix' for
            manual override)

    Returns
    -------
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
                f"  {score_name}: calculated for {len(valid_scores)}/{len(df)} "
                f"molecules (mean={valid_scores.mean():.2f}, "
                f"std={valid_scores.std():.2f})"
            )
        else:
            logger.debug(
                f"  {score_name}: could not be calculated (module not available)"
            )

    return result_df


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

    Returns
    -------
        Filtered DataFrame with molecules that pass all applicable filters
    """
    filtered_df = df.copy()
    pass_mask = pd.Series([True] * len(filtered_df), index=filtered_df.index)

    # Apply SA score filter
    pass_mask = _apply_score_filter(
        filtered_df, pass_mask, "sa_score", config, "SA score filter"
    )

    # Apply RA score filter
    pass_mask = _apply_score_filter(
        filtered_df, pass_mask, "ra_score", config, "RA score filter"
    )

    # Apply SYBA score filter
    pass_mask = _apply_score_filter(
        filtered_df, pass_mask, "syba_score", config, "SYBA score filter"
    )

    return filtered_df[pass_mask]


def _apply_score_filter(
    df: pd.DataFrame,
    pass_mask: pd.Series,
    score_name: str,
    config: dict[str, Any],
    filter_name: str,
) -> pd.Series:
    """Apply a single score filter to the dataframe."""
    if score_name not in df.columns:
        return pass_mask

    valid_scores = df[score_name].dropna()
    if len(valid_scores) == 0:
        logger.info(f"{filter_name}: skipped (no valid scores calculated)")
        return pass_mask

    score_min = config.get(f"{score_name}_min", 0)
    score_max = config.get(f"{score_name}_max", "inf")

    if score_max != "inf":
        score_mask = df[score_name].isna() | (
            (df[score_name] >= score_min) & (df[score_name] <= score_max)
        )
    else:
        score_mask = df[score_name].isna() | (df[score_name] >= score_min)

    return pass_mask & score_mask
