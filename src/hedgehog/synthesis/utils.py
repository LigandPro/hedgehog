import json
import os
import pickle
import queue
import re
import shutil
import subprocess
import tempfile
import threading
import time
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from hedgehog.configs.logger import logger
from hedgehog.setup._download import resolve_uv_binary
from hedgehog.struct_filters.utils import process_path
from hedgehog.utils.input_paths import get_all_input_candidates
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

# Lazy-loaded module cache
_lazy_cache: dict[str, Any] = {
    "sascorer": None,
    "syba_model": None,
    "rascore_booster": None,
}
STRICT_RASCORE_ENV = "HEDGEHOG_STRICT_RASCORE"

_DEFAULT_AIZYNTH_LOGGING_YML = """version: 1
formatters:
  simple:
    format: '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
handlers:
  console:
    class: logging.StreamHandler
    level: INFO
    formatter: simple
    stream: ext://sys.stdout
  file:
    class: logging.FileHandler
    level: DEBUG
    formatter: simple
    filename: aizynthfinder.log
loggers:
  aizynthfinder:
    level: DEBUG
    handlers: [console, file]
    propagate: no
root:
  level: ERROR
  handlers: [console]
"""

_AIZYNTH_PAIR_PROGRESS_RE = re.compile(r"(?<!\d)(\d+)\s*/\s*(\d+)(?!\d)")
_AIZYNTH_PERCENT_PROGRESS_RE = re.compile(r"(\d{1,3}(?:\.\d+)?)\s*%")


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


def _short_exception_message(exc: Exception) -> str:
    """Return a compact one-line message for noisy exceptions."""
    raw = str(exc).strip()
    if not raw:
        return exc.__class__.__name__
    return raw.splitlines()[0]


def _load_xgboost_booster(json_path: Path) -> Any:
    """Load an XGBoost booster while handling known legacy-model warnings."""
    import xgboost as xgb

    booster = xgb.Booster()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UserWarning)
        booster.load_model(str(json_path))

    legacy_warning_fragments = ("found json model saved", "xgboost 1.6")
    for warning in caught:
        message = str(warning.message)
        normalized_message = " ".join(message.lower().split())
        if all(part in normalized_message for part in legacy_warning_fragments):
            logger.debug(
                "Loaded legacy RAScore JSON model format from %s",
                json_path,
            )
            continue
        warnings.warn(warning.message, warning.category, stacklevel=2)

    return booster


def _ensure_aizynth_logging_config(run_dir: Path) -> None:
    """Ensure AiZynthFinder logging.yml exists for CLI startup."""
    data_dir = run_dir / "aizynthfinder" / "data"
    logging_path = data_dir / "logging.yml"
    if logging_path.exists():
        return

    data_dir.mkdir(parents=True, exist_ok=True)
    source_candidates = [
        Path(__file__).resolve().parent / "logging.yml",
        Path(__file__).resolve().parents[1] / "stages" / "synthesis" / "logging.yml",
    ]
    source = next((p for p in source_candidates if p.exists()), None)
    if source is not None:
        shutil.copy2(source, logging_path)
        logger.info("Copied AiZynthFinder logging config to %s", logging_path)
        return

    logging_path.write_text(_DEFAULT_AIZYNTH_LOGGING_YML, encoding="utf-8")
    logger.warning(
        "AiZynthFinder logging.yml was missing; created default config at %s",
        logging_path,
    )


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


def _count_nonempty_lines(path: Path) -> int:
    """Count non-empty lines in a text file (best-effort)."""
    try:
        with path.open(encoding="utf-8", errors="ignore") as handle:
            return sum(1 for line in handle if line.strip())
    except Exception:
        return 0


def _extract_aizynth_progress(
    line: str, fallback_total: int | None = None
) -> tuple[int, int] | None:
    """Extract progress tuple ``(done, total)`` from an AiZynthFinder output line."""
    best: tuple[int, int] | None = None

    for done_raw, total_raw in _AIZYNTH_PAIR_PROGRESS_RE.findall(line):
        done = int(done_raw)
        total = int(total_raw)
        if total <= 0:
            continue
        if done < 0:
            continue
        done = min(done, total)
        if best is None or done > best[0]:
            best = (done, total)

    if best is not None:
        return best

    if fallback_total is None or fallback_total <= 0:
        return None

    percent_matches = _AIZYNTH_PERCENT_PROGRESS_RE.findall(line)
    if not percent_matches:
        return None

    pct = float(percent_matches[-1])
    if pct < 0:
        pct = 0.0
    if pct > 100:
        pct = 100.0
    done = int(round((pct / 100.0) * fallback_total))
    done = min(max(done, 0), fallback_total)
    return done, fallback_total


def run_aizynthfinder(
    input_smiles_file,
    output_json_file,
    config_file,
    aizynthfinder_dir=None,
    synthesis_config: dict[str, Any] | None = None,
    progress_cb: Callable[[int, int, str | None], None] | None = None,
):
    """Run AiZynthFinder CLI to perform retrosynthesis analysis.

    Args:
        input_smiles_file: Path to input SMILES file
        output_json_file: Path to output JSON file
        config_file: Path to AiZynthFinder config file
        aizynthfinder_dir: Optional working directory for AiZynthFinder CLI
        synthesis_config: Optional synthesis stage config (supports ``n_jobs``)

    Returns:
        True if successful, False otherwise
    """
    output_json_file.parent.mkdir(parents=True, exist_ok=True)

    run_dir = (
        Path(aizynthfinder_dir) if aizynthfinder_dir else config_file.parent.parent
    )
    _ensure_aizynth_logging_config(run_dir)
    input_abs = input_smiles_file.resolve()
    output_abs = output_json_file.resolve()
    config_abs = config_file.resolve()

    try:
        uv_binary = resolve_uv_binary()
    except RuntimeError as exc:
        logger.error(
            "Unable to start retrosynthesis: uv executable could not be resolved (%s)",
            exc,
        )
        return False

    cmd = [
        uv_binary,
        "run",
        "aizynthcli",
        "--config",
        str(config_abs),
        "--smiles",
        str(input_abs),
        "--output",
        str(output_abs),
    ]
    nproc = _resolve_aizynth_nproc(
        os.environ.get("AIZYNTH_NPROC"), synthesis_config=synthesis_config
    )
    cmd.extend(["--nproc", str(nproc)])

    try:
        logger.info("Running retrosynthesis analysis...")
        logger.debug("Command: %s", cmd)
        child_env = os.environ.copy()
        child_env.pop("VIRTUAL_ENV", None)

        if progress_cb is None:
            subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=str(run_dir),
                env=child_env,
            )
            logger.info("Retrosynthesis analysis completed successfully")
            return True

        total_targets = _count_nonempty_lines(input_abs)
        if total_targets <= 0:
            total_targets = 1
        progress_cb(0, total_targets, "starting")

        proc = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=str(run_dir),
            env=child_env,
        )
        assert proc.stdout is not None

        line_queue: queue.Queue[str | None] = queue.Queue()

        def _reader() -> None:
            try:
                for raw_line in proc.stdout:
                    line_queue.put(raw_line.rstrip("\r\n"))
            finally:
                line_queue.put(None)

        reader_thread = threading.Thread(
            target=_reader, name="aizynth-output-reader", daemon=True
        )
        reader_thread.start()

        best_done = 0
        best_total = total_targets
        last_emit = time.monotonic()
        start_time = last_emit
        stream_open = True

        while stream_open:
            emitted = False
            while True:
                try:
                    item = line_queue.get_nowait()
                except queue.Empty:
                    break
                if item is None:
                    stream_open = False
                    break
                line = item.strip()
                if line:
                    logger.debug("AiZynthFinder: %s", line[:300])
                parsed = _extract_aizynth_progress(line, fallback_total=best_total)
                if parsed is not None:
                    done, total = parsed
                    if total > 0:
                        if total == best_total:
                            mapped_done = min(done, best_total)
                        else:
                            ratio = max(0.0, min(float(done) / float(total), 1.0))
                            mapped_done = int(round(ratio * best_total))
                        mapped_done = max(0, min(mapped_done, best_total))
                        if mapped_done > best_done:
                            best_done = mapped_done
                            progress_cb(
                                best_done, best_total, line[:80] if line else None
                            )
                            emitted = True

            now = time.monotonic()
            # Heartbeat fallback for runs that do not expose parseable progress.
            if not emitted and now - last_emit >= 1.0:
                elapsed = now - start_time
                if best_total > 1 and best_done < best_total - 1:
                    synthetic = min(int(elapsed // 2), best_total - 1)
                    if synthetic > best_done:
                        best_done = synthetic
                progress_cb(
                    best_done,
                    best_total,
                    f"elapsed {int(elapsed)}s",
                )
                last_emit = now

            if stream_open:
                time.sleep(0.1)

        return_code = proc.wait()
        reader_thread.join(timeout=1.0)
        if return_code != 0:
            logger.error(
                "Retrosynthesis analysis failed with exit code %d", return_code
            )
            return False

        logger.info("Retrosynthesis analysis completed successfully")
        progress_cb(best_total, best_total, "completed")
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


def _resolve_aizynth_nproc(
    raw_value: str | None, synthesis_config: dict[str, Any] | None = None
) -> int:
    """Resolve AiZynthFinder worker count to a valid positive value.

    Semantics:
      - ``AIZYNTH_NPROC`` positive int: use as-is
      - ``AIZYNTH_NPROC`` ``0`` or negative: resolve as auto (all available CPUs)
      - invalid ``AIZYNTH_NPROC``: ignore and fall back
      - fallback source: ``synthesis_config['n_jobs']`` when provided
      - final fallback: auto (all available CPUs)
    """
    config_n_jobs: int = -1
    if synthesis_config is not None and "n_jobs" in synthesis_config:
        try:
            config_n_jobs = int(synthesis_config["n_jobs"])
        except (TypeError, ValueError):
            logger.warning(
                "Invalid synthesis.n_jobs=%r, falling back to auto",
                synthesis_config["n_jobs"],
            )
            config_n_jobs = -1

    if not raw_value:
        resolved = resolve_n_jobs({"n_jobs": config_n_jobs})
        logger.info(
            "AIZYNTH_NPROC is not set; using %d workers from synthesis.n_jobs/default",
            resolved,
        )
        return resolved

    try:
        parsed = int(raw_value)
    except ValueError:
        logger.warning("Invalid AIZYNTH_NPROC=%r, ignoring", raw_value)
        resolved = resolve_n_jobs({"n_jobs": config_n_jobs})
        logger.info(
            "Using %d workers from synthesis.n_jobs/default after invalid AIZYNTH_NPROC",
            resolved,
        )
        return resolved

    if parsed > 0:
        return parsed

    resolved = resolve_n_jobs({"n_jobs": parsed})
    logger.info("Resolved AIZYNTH_NPROC=%s to %d workers", raw_value, resolved)
    return resolved


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

    Uses SMILES-based matching when a SMILES column is available in both
    DataFrames, falling back to positional merge otherwise.

    Args:
        input_df: Original input DataFrame with molecules (may have duplicate SMILES)
        retrosynth_df: DataFrame with retrosynthesis results

    Returns:
        Merged DataFrame with retrosynthesis information, preserving all input rows
    """
    merged = input_df.reset_index(drop=True).copy()
    retrosynth_df_copy = retrosynth_df.reset_index(drop=True)

    merged["solved"] = 0
    merged["search_time"] = 0.0

    # Determine SMILES column names
    input_smi_col = None
    for col in ("smiles", "SMILES"):
        if col in merged.columns:
            input_smi_col = col
            break

    retro_smi_col = None
    for col in ("SMILES", "smiles"):
        if col in retrosynth_df_copy.columns:
            retro_smi_col = col
            break

    if input_smi_col and retro_smi_col:
        # SMILES-based merge: build lookup from retrosynth results
        retro_lookup: dict[str, dict] = {}
        for _, row in retrosynth_df_copy.iterrows():
            smi = str(row[retro_smi_col])
            if smi not in retro_lookup:
                retro_lookup[smi] = {
                    "solved": row.get("solved", 0),
                    "search_time": row.get("search_time", 0.0),
                }

        for idx, row in merged.iterrows():
            smi = str(row[input_smi_col])
            match = retro_lookup.get(smi)
            if match:
                merged.loc[idx, "solved"] = match["solved"]
                merged.loc[idx, "search_time"] = match["search_time"]
    else:
        # Positional fallback when SMILES column is unavailable
        logger.warning("No SMILES column found for merge; using positional matching")
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
    """Load RAScore model.

    Preferred order:
      1) native XGBoost JSON model (if present and valid)
      2) MolScore upstream pickle model
    """
    json_path = _get_rascore_model_path()
    pkl_path = _get_rascore_pickle_model_path()

    if not json_path.exists():
        ensured_json = _ensure_rascore_json_model()
        if ensured_json is not None:
            json_path = ensured_json

    if json_path.exists():
        try:
            booster = _load_xgboost_booster(json_path)
            logger.debug("RAScore xgboost model loaded successfully")
            return {"kind": "xgboost_booster", "model": booster}
        except Exception as e:
            logger.warning(
                "Failed to load RAScore model from %s: %s",
                json_path,
                _short_exception_message(e),
            )

    if not pkl_path.exists():
        try:
            from hedgehog.setup import ensure_rascore_model

            project_root = Path(__file__).resolve().parents[3]
            pkl_path = ensure_rascore_model(project_root)
        except Exception as e:
            logger.warning(
                "RAScore model is unavailable (%s). RA scores will be set to np.nan",
                e,
            )
            return False

    try:
        with open(pkl_path, "rb") as f:
            model = pickle.load(f)
        logger.debug("RAScore pickle model loaded successfully")
        return {"kind": "pickle_classifier", "model": model}
    except Exception as e:
        logger.warning(
            "Failed to load RAScore model from %s: %s",
            pkl_path,
            _short_exception_message(e),
        )
        return False


def _load_sascorer():
    """Load SA Score module with caching."""
    return _get_cached("sascorer", _load_sascorer_impl)


def _load_syba_model():
    """Load SYBA model with caching."""
    return _get_cached("syba_model", _load_syba_model_impl)


def _load_rascore():
    """Load RAScore model with caching."""
    return _get_cached("rascore_booster", _load_rascore_impl)


def _strict_rascore_enabled() -> bool:
    """Return True when strict RAScore mode is enabled via environment."""
    raw = os.environ.get(STRICT_RASCORE_ENV, "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def _calculate_sa_score_single(smiles: str) -> float | None:
    """Calculate SA Score for a single SMILES string.

    Designed as a top-level picklable function for use with parallel_map.
    With fork context, the parent's _lazy_cache is inherited by workers.
    """
    sascorer = _load_sascorer()
    if not sascorer:
        return np.nan
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
        return sascorer.calculateScore(mol)
    except Exception:
        return np.nan


def _calculate_syba_score_single(smiles: str) -> float | None:
    """Calculate SYBA Score for a single SMILES string.

    Designed as a top-level picklable function for use with parallel_map.
    With fork context, the parent's _lazy_cache is inherited by workers.
    """
    syba_model = _load_syba_model()
    if not syba_model:
        return np.nan
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.nan
        return syba_model.predict(mol=mol)
    except Exception:
        return np.nan


def _get_rascore_model_path() -> Path:
    """Get path to RAScore model file (native xgboost JSON format)."""
    molscore_path = Path(__file__).parent.parent.parent.parent / "modules" / "MolScore"
    return (
        molscore_path
        / "molscore"
        / "data"
        / "models"
        / "RAScore"
        / "XGB_chembl_ecfp_counts"
        / "model.json"
    )


def _get_rascore_pickle_model_path() -> Path:
    """Get path to RAScore model file in MolScore upstream pickle format."""
    molscore_path = Path(__file__).parent.parent.parent.parent / "modules" / "MolScore"
    return (
        molscore_path
        / "molscore"
        / "data"
        / "models"
        / "RAScore"
        / "XGB_chembl_ecfp_counts"
        / "model.pkl"
    )


def _get_ecfp6_counts(smiles: str, n_bits: int = 2048) -> np.ndarray | None:
    """Compute ECFP6 count fingerprint for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprint(mol, 3, useFeatures=False, useCounts=True)
    arr = np.zeros((n_bits,), dtype=np.int32)
    for idx, v in fp.GetNonzeroElements().items():
        arr[idx % n_bits] += v
    return arr


def _predict_with_pickle_classifier(model: Any, fps_array: np.ndarray) -> np.ndarray:
    """Predict RA probabilities from a pickled XGBClassifier-like model.

    Primary path uses ``predict_proba``. When the pickled object is from an older
    xgboost+sklearn combination and ``predict_proba`` is incompatible with the
    current sklearn API, fall back to direct booster prediction.
    """
    try:
        return model.predict_proba(fps_array)[:, 1]
    except Exception as e:
        logger.debug("RAScore predict_proba failed (%s). Trying booster fallback.", e)
        import xgboost as xgb

        if not hasattr(model, "get_booster"):
            raise
        booster = model.get_booster()
        return booster.predict(xgb.DMatrix(fps_array))


def _ensure_rascore_pickle_path() -> Path | None:
    """Ensure pickle RAScore model exists and return its path."""
    model_path = _get_rascore_pickle_model_path()
    if model_path.exists():
        return model_path
    try:
        from hedgehog.setup import ensure_rascore_model

        project_root = Path(__file__).resolve().parents[3]
        return ensure_rascore_model(project_root)
    except Exception as e:
        logger.warning(
            "RAScore pickle model is unavailable (%s). RA scores will be set to np.nan",
            e,
        )
        return None


def _convert_rascore_pickle_to_json(pkl_path: Path, json_path: Path) -> bool:
    """Convert legacy RAScore pickle model to portable XGBoost JSON."""
    try:
        uv_binary = resolve_uv_binary()
    except RuntimeError as e:
        logger.warning(
            "Unable to convert RAScore model to JSON: uv executable could not be "
            "resolved (%s)",
            e,
        )
        return False

    project_root = Path(__file__).resolve().parents[3]
    cmd = [
        uv_binary,
        "run",
        "--with",
        "xgboost==1.5.2",
        "--with",
        "setuptools<81",
        "python",
        "-m",
        "hedgehog.workers.rascore_convert_worker",
        "--model-pkl",
        str(pkl_path),
        "--output-json",
        str(json_path),
    ]
    env = os.environ.copy()
    env.setdefault("UV_NO_SYNC", "1")

    try:
        subprocess.run(
            cmd,
            cwd=str(project_root),
            env=env,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
        )
    except Exception as e:
        logger.warning("Failed to convert RAScore pickle model to JSON: %s", e)
        json_path.unlink(missing_ok=True)
        return False

    if not json_path.exists() or json_path.stat().st_size < 1_024:
        logger.warning(
            "Converted RAScore JSON model is missing or invalid: %s",
            json_path,
        )
        json_path.unlink(missing_ok=True)
        return False

    logger.info("RAScore JSON model generated successfully: %s", json_path)
    return True


def _ensure_rascore_json_model() -> Path | None:
    """Ensure JSON RAScore model exists and return its path."""
    json_path = _get_rascore_model_path()
    if json_path.exists():
        return json_path

    pkl_path = _ensure_rascore_pickle_path()
    if pkl_path is None:
        return None

    if _convert_rascore_pickle_to_json(pkl_path, json_path):
        return json_path
    return None


def _calculate_ra_scores_batch_legacy(
    smiles_list: list[str],
    progress_cb: Callable[[int, int], None] | None = None,
) -> list[float]:
    """Calculate RA scores via legacy worker with xgboost 1.5 compatibility."""
    nan_list = [np.nan] * len(smiles_list)
    if not smiles_list:
        return nan_list
    total = len(smiles_list)
    if progress_cb is not None:
        progress_cb(0, total)

    model_path = _ensure_rascore_pickle_path()
    if model_path is None:
        if progress_cb is not None:
            progress_cb(total, total)
        return nan_list

    try:
        uv_binary = resolve_uv_binary()
    except RuntimeError as e:
        logger.warning(
            "Unable to run legacy RAScore worker: uv executable could not be resolved (%s)",
            e,
        )
        if progress_cb is not None:
            progress_cb(total, total)
        return nan_list

    project_root = Path(__file__).resolve().parents[3]
    with tempfile.TemporaryDirectory(prefix="rascore_worker_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        input_json = tmp_path / "rascore_input.json"
        output_json = tmp_path / "rascore_output.json"
        input_json.write_text(json.dumps({"smiles": smiles_list}), encoding="utf-8")

        cmd = [
            uv_binary,
            "run",
            "--with",
            "xgboost==1.5.2",
            "--with",
            "setuptools<81",
            "python",
            "-m",
            "hedgehog.workers.rascore_worker",
            "--model-pkl",
            str(model_path),
            "--input-json",
            str(input_json),
            "--output-json",
            str(output_json),
        ]
        env = os.environ.copy()
        env.setdefault("UV_NO_SYNC", "1")

        try:
            subprocess.run(
                cmd,
                cwd=str(project_root),
                env=env,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
            )
            payload = json.loads(output_json.read_text(encoding="utf-8"))
            raw_scores = payload.get("scores", [])
        except Exception as e:
            logger.warning("Legacy RAScore worker failed: %s", e)
            if progress_cb is not None:
                progress_cb(total, total)
            return nan_list

    if len(raw_scores) != len(smiles_list):
        logger.warning(
            "Legacy RAScore worker returned %d scores for %d molecules",
            len(raw_scores),
            len(smiles_list),
        )
        if progress_cb is not None:
            progress_cb(total, total)
        return nan_list

    scores: list[float] = []
    for value in raw_scores:
        if value is None:
            scores.append(np.nan)
        else:
            try:
                scores.append(float(value))
            except (TypeError, ValueError):
                scores.append(np.nan)
    if progress_cb is not None:
        progress_cb(total, total)
    return scores


def _calculate_ra_scores_batch(
    smiles_list: list,
    config: dict[str, Any] | None = None,
    progress_cb: Callable[[int, int], None] | None = None,
) -> list:
    """Calculate RA scores for multiple molecules in a single batch.

    Args:
        smiles_list: List of SMILES strings
        config: Not used, kept for API compatibility

    Returns:
        List of RA scores (0-1, higher is better), with np.nan for failed calculations
    """
    nan_list = [np.nan] * len(smiles_list)
    total = len(smiles_list)
    if total <= 0:
        return nan_list
    if progress_cb is not None:
        progress_cb(0, total)

    strict_rascore = _strict_rascore_enabled()
    rascore_model = _load_rascore()
    if not rascore_model:
        if strict_rascore:
            raise RuntimeError(
                "RAScore local model load failed in strict mode "
                f"({STRICT_RASCORE_ENV}=1)."
            )
        logger.info(
            "RAScore local model load failed; falling back to legacy worker backend"
        )
        return _calculate_ra_scores_batch_legacy(smiles_list, progress_cb=progress_cb)

    fps = []
    valid_indices = []
    progress_step = max(1, min(1_000, total // 200)) if total > 0 else 1
    next_progress_log = progress_step
    for i, smi in enumerate(smiles_list, start=1):
        fp = _get_ecfp6_counts(smi)
        if fp is not None:
            fps.append(fp)
            valid_indices.append(i - 1)
        if progress_cb is not None and (i == total or i >= next_progress_log):
            progress_cb(i, total)
            while next_progress_log <= i:
                next_progress_log += progress_step

    if not fps:
        if progress_cb is not None:
            progress_cb(total, total)
        return nan_list

    try:
        fps_array = np.array(fps)
        if rascore_model.get("kind") == "pickle_classifier":
            probs = _predict_with_pickle_classifier(rascore_model["model"], fps_array)
        else:
            import xgboost as xgb

            dmatrix = xgb.DMatrix(fps_array)
            probs = rascore_model["model"].predict(dmatrix)
        scores = [np.nan] * len(smiles_list)
        for idx, prob in zip(valid_indices, probs):
            scores[idx] = float(prob)
        if progress_cb is not None:
            progress_cb(total, total)
        return scores
    except Exception as e:
        if strict_rascore:
            raise RuntimeError(
                "Failed to calculate RA scores in strict mode "
                f"({STRICT_RASCORE_ENV}=1)."
            ) from e
        logger.debug("Failed to calculate RA scores in batch: %s", e)
        if progress_cb is not None:
            progress_cb(total, total)
        return nan_list


def calculate_synthesis_scores(df, folder_to_save=None, config=None, progress_cb=None):
    """Calculate SA, SYBA, and RA scores for all molecules in DataFrame.

    Args:
        df: DataFrame with 'smiles' column
        folder_to_save: Optional folder to save outputs
        config: Optional config dict

    Returns:
        DataFrame with added columns: sa_score, syba_score, ra_score
    """
    logger.info("Calculating synthetic accessibility scores...")
    _load_sascorer()
    _load_syba_model()
    _load_rascore()

    n_jobs = resolve_n_jobs(config)
    result_df = df.copy()
    smiles_list = result_df["smiles"].tolist()

    sa_progress = None
    if progress_cb is not None:

        def _sa_progress(done: int, total: int) -> None:
            progress_cb("sa_score", done, total)

        sa_progress = _sa_progress

    result_df["sa_score"] = parallel_map(
        _calculate_sa_score_single, smiles_list, n_jobs, progress=sa_progress
    )

    # SYBA uses GPU internally — always run sequentially to avoid forking issues
    syba_progress = None
    if progress_cb is not None:

        def _syba_progress(done: int, total: int) -> None:
            progress_cb("syba_score", done, total)

        syba_progress = _syba_progress

    result_df["syba_score"] = parallel_map(
        _calculate_syba_score_single, smiles_list, n_jobs=1, progress=syba_progress
    )

    ra_progress = None
    if progress_cb is not None:

        def _ra_progress(done: int, total: int) -> None:
            progress_cb("ra_score", done, total)

        ra_progress = _ra_progress

    ra_scores = _calculate_ra_scores_batch(smiles_list, config, progress_cb=ra_progress)
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
