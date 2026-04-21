import csv
import inspect
import multiprocessing
import shutil
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from hedgehog.configs.logger import load_config, logger
from hedgehog.descriptors.stage import run as descriptors_main
from hedgehog.docking.main import main as docking_main
from hedgehog.docking_filters.main import docking_filters_main
from hedgehog.molprep.main import main as mol_prep_main
from hedgehog.reporting import ReportGenerator
from hedgehog.struct_filters.main import main as structural_filters_main
from hedgehog.synthesis.main import main as synthesis_main
from hedgehog.utils.input_paths import find_latest_input_source as _find_input

# Directory names
DIR_INPUT = "input"
DIR_STAGES = "stages"
DIR_OUTPUT = "output"
DIR_CONFIGS = "configs"

# Stage subdirectories
DIR_MOL_PREP = "stages/00_mol_prep"
DIR_DESCRIPTORS_INITIAL = "stages/01_descriptors_initial"
DIR_STRUCT_FILTERS_POST = "stages/03_structural_filters_post"
DIR_SYNTHESIS = "stages/04_synthesis"
DIR_DOCKING = "stages/05_docking"
DIR_DOCKING_FILTERS = "stages/06_docking_filters"
DIR_DESCRIPTORS_FINAL = "stages/07_descriptors_final"

# Legacy names for backwards compatibility
DIR_DESCRIPTORS = DIR_DESCRIPTORS_INITIAL
DIR_STRUCT_FILTERS = DIR_STRUCT_FILTERS_POST
DIR_FINAL_DESCRIPTORS = DIR_DESCRIPTORS_FINAL
DIR_RUN_CONFIGS = DIR_CONFIGS

# File names
FILE_SAMPLED_MOLECULES = "sampled_molecules.csv"
FILE_FINAL_MOLECULES = "final_molecules.csv"
FILE_FILTERED_MOLECULES = "filtered_molecules.csv"
FILE_PASS_SMILES_TEMPLATE = "filtered_molecules.csv"
FILE_MASTER_CONFIG = "master_config_resolved.yml"
FILE_GNINA_OUTPUT = "gnina_out.sdf"
FILE_SMINA_OUTPUT = "smina_out.sdf"
FILE_MATCHA_OUTPUT = "matcha_out.sdf"
FILE_DOCKING_COMPLETED_EMPTY = "completed_empty.marker"

DOCKING_SCORE_COLUMNS = [
    "gnina_affinity",
    "gnina_cnnscore",
    "gnina_cnnaffinity",
    "gnina_cnn_vs",
    "smina_affinity",
    "matcha_affinity",
]

# Stage names
STAGE_MOL_PREP = "mol_prep"
STAGE_DESCRIPTORS = "descriptors"
STAGE_STRUCT_FILTERS = "struct_filters"
STAGE_SYNTHESIS = "synthesis"
STAGE_DOCKING = "docking"
STAGE_DOCKING_FILTERS = "docking_filters"
STAGE_FINAL_DESCRIPTORS = "final_descriptors"
STAGE_REPORT = "report"

# Config keys
CONFIG_MOL_PREP = "config_mol_prep"
CONFIG_DESCRIPTORS = "config_descriptors"
CONFIG_STRUCT_FILTERS = "config_structFilters"
CONFIG_SYNTHESIS = "config_synthesis"
CONFIG_DOCKING = "config_docking"
CONFIG_DOCKING_FILTERS = "config_docking_filters"
CONFIG_RUN_KEY = "run"
CONFIG_TOOLS = "tools"
CONFIG_FOLDER_TO_SAVE = "folder_to_save"

# Command-line override keys
OVERRIDE_SINGLE_STAGE = "_run_single_stage_override"
OVERRIDE_STAGE_SELECTION = "_run_stage_selection_override"

# Docking tools
DOCKING_TOOL_SMINA = "smina"
DOCKING_TOOL_GNINA = "gnina"
DOCKING_TOOL_MATCHA = "matcha"
DOCKING_TOOL_ALL = "all"
DOCKING_TOOL_BOTH = "both"
DOCKING_RESULTS_DIR_TEMPLATE = {
    DOCKING_TOOL_SMINA: f"{DIR_DOCKING}/smina",
    DOCKING_TOOL_GNINA: f"{DIR_DOCKING}/gnina",
    DOCKING_TOOL_MATCHA: f"{DIR_DOCKING}/matcha",
}


FILE_RUN_INCOMPLETE = ".RUN_INCOMPLETE"


def _plain_output_enabled() -> bool:
    import os

    return os.environ.get("HEDGEHOG_PLAIN_OUTPUT", "").strip() == "1"


class StageProgressReporter:
    """Emit structured progress events for a single pipeline stage."""

    def __init__(
        self,
        emit_event,
        stage: str,
        stage_index: int,
        total_stages: int,
    ) -> None:
        self._emit_event = emit_event
        self.stage = stage
        self.stage_index = stage_index
        self.total_stages = total_stages
        self._last_progress_signature: tuple[int, int, str | None] | None = None

    def start(
        self,
        message: str | None = None,
        molecules_in: int | None = None,
        molecules_out: int | None = None,
    ) -> None:
        self._last_progress_signature = None
        self._emit_event(
            {
                "type": "stage_start",
                "stage": self.stage,
                "stage_index": self.stage_index,
                "total_stages": self.total_stages,
                "message": message,
                "molecules_in": molecules_in,
                "molecules_out": molecules_out,
            }
        )

    def progress(
        self,
        current: int,
        total: int,
        message: str | None = None,
        molecules_in: int | None = None,
        molecules_out: int | None = None,
    ) -> None:
        signature = (int(current), int(total), message)
        if signature == self._last_progress_signature:
            return
        self._last_progress_signature = signature
        self._emit_event(
            {
                "type": "stage_progress",
                "stage": self.stage,
                "stage_index": self.stage_index,
                "total_stages": self.total_stages,
                "current": signature[0],
                "total": signature[1],
                "message": message,
                "molecules_in": molecules_in,
                "molecules_out": molecules_out,
            }
        )

    def complete(
        self,
        ok: bool = True,
        message: str | None = None,
        molecules_in: int | None = None,
        molecules_out: int | None = None,
        elapsed_seconds: float | None = None,
    ) -> None:
        self._last_progress_signature = None
        self._emit_event(
            {
                "type": "stage_complete",
                "stage": self.stage,
                "stage_index": self.stage_index,
                "total_stages": self.total_stages,
                "ok": bool(ok),
                "message": message,
                "molecules_in": molecules_in,
                "molecules_out": molecules_out,
                "elapsed_seconds": elapsed_seconds,
            }
        )


def _log_stage_header(stage_label: str) -> None:
    """Log a formatted stage header."""
    logger.info("")
    if _plain_output_enabled():
        logger.info("=== %s ===", stage_label)
        logger.info("")
        return

    separator = "[dim]" + "\u2500" * 59 + "[/dim]"
    logger.info(separator)
    logger.info("[bold]  %s[/bold]", stage_label)
    logger.info(separator)
    logger.info("")


def _format_molecule_count(value: int | None) -> str:
    """Return a human-readable molecule count string."""
    if value is None:
        return "unknown"
    return f"{value:,}"


def _format_delta(molecules_in: int | None, molecules_out: int | None) -> str:
    """Return the stage delta in signed format."""
    if molecules_in is None or molecules_out is None:
        return "unknown"
    return f"{(molecules_out - molecules_in):+,}"


def _format_retention(molecules_in: int | None, molecules_out: int | None) -> str:
    """Return stage retention percentage with two decimals."""
    if molecules_in is None or molecules_out is None or molecules_in <= 0:
        return "unknown"
    return f"{(100.0 * molecules_out / molecules_in):.2f}%"


def _build_stage_start_message(stage_name: str, molecules_in: int | None) -> str:
    """Build canonical stage-start log message."""
    return (
        f"Stage {stage_name} started: "
        f"{_format_molecule_count(molecules_in)} molecules in."
    )


def _build_stage_complete_message(
    stage_name: str,
    ok: bool,
    molecules_in: int | None,
    molecules_out: int | None,
    elapsed_seconds: float,
    avg_cpu_percent: float | None = None,
    output_estimated: bool = False,
) -> str:
    """Build canonical stage-complete log message."""
    status = "completed" if ok else "failed"
    parts = [
        f"delta {_format_delta(molecules_in, molecules_out)}",
        f"retained {_format_retention(molecules_in, molecules_out)}",
    ]
    if output_estimated and molecules_out is not None:
        parts.insert(0, "estimated out")
    if avg_cpu_percent is not None:
        parts.append(f"avg CPU {avg_cpu_percent:.1f}%")
    parts.append(f"{elapsed_seconds:.1f}s")

    return (
        f"Stage {stage_name} {status}: "
        f"{_format_molecule_count(molecules_in)} in -> "
        f"{_format_molecule_count(molecules_out)} out "
        f"({', '.join(parts)})."
    )


def _accepts_reporter(runner_func) -> bool:
    """Return whether a stage runner accepts the reporter keyword."""
    signature = inspect.signature(runner_func)
    return "reporter" in signature.parameters or any(
        param.kind == inspect.Parameter.VAR_KEYWORD
        for param in signature.parameters.values()
    )


def _file_exists_and_not_empty(file_path: Path) -> bool:
    """Check if a file exists and is not empty."""
    try:
        return file_path.exists() and file_path.stat().st_size > 0
    except OSError:
        return False


def _csv_has_data_rows(file_path: Path) -> bool:
    """Return True if a CSV file has at least one data row (header-only -> False)."""
    try:
        if not file_path.exists():
            return False
        with file_path.open(encoding="utf-8", errors="ignore") as f:
            # First non-empty line is assumed to be a header
            for line in f:
                if line.strip():
                    break
            # Any subsequent non-empty line counts as data
            for line in f:
                if line.strip():
                    return True
        return False
    except OSError:
        return False


def _directory_has_files(dir_path: Path) -> bool:
    """Check if a directory exists and contains files."""
    try:
        return dir_path.exists() and any(p.is_file() for p in dir_path.iterdir())
    except OSError:
        return False


def _extract_float_prop(mol, prop_names: list[str]) -> float | None:
    """Return the first parseable float property from an SDF molecule."""
    for prop_name in prop_names:
        if not mol.HasProp(prop_name):
            continue
        try:
            return float(mol.GetProp(prop_name))
        except (TypeError, ValueError):
            continue
    return None


def _resolve_gnina_output_sdf(base_path: Path, docking_cfg: dict | None = None) -> Path:
    """Resolve GNINA output SDF path from docking config."""
    cfg = docking_cfg or {}
    gnina_cfg = cfg.get("gnina_config", {}) or {}
    cfg_out_dir = gnina_cfg.get("output_dir") or cfg.get("gnina_output_dir")
    if cfg_out_dir:
        out_dir = Path(str(cfg_out_dir))
        if not out_dir.is_absolute():
            out_dir = base_path / out_dir
        return out_dir / FILE_GNINA_OUTPUT
    return (
        base_path / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_GNINA] / FILE_GNINA_OUTPUT
    )


def _find_docking_sdf(
    base_path: Path, tool_name: str, docking_cfg: dict | None = None
) -> Path | None:
    """Find docking output SDF for a docking tool."""
    docking_dir = base_path / DIR_DOCKING
    if tool_name == DOCKING_TOOL_GNINA:
        candidates = [
            _resolve_gnina_output_sdf(base_path, docking_cfg),
            docking_dir / f"{tool_name}_out.sdf",
        ]
    elif tool_name == DOCKING_TOOL_SMINA:
        candidates = [
            docking_dir / tool_name / FILE_SMINA_OUTPUT,
            docking_dir / f"{tool_name}_out.sdf",
        ]
    elif tool_name == DOCKING_TOOL_MATCHA:
        candidates = [
            docking_dir / tool_name / FILE_MATCHA_OUTPUT,
            docking_dir / f"{tool_name}_out.sdf",
        ]
    else:
        candidates = [
            docking_dir / tool_name / f"{tool_name}_out.sdf",
            docking_dir / f"{tool_name}_out.sdf",
        ]

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _extract_best_tool_scores(sdf_path: Path, tool_name: str) -> dict[str, dict]:
    """Extract best (lowest affinity) pose scores per molecule from docking SDF."""
    try:
        from rdkit import Chem
    except ImportError:
        logger.warning(
            "RDKit is not available. Docking score enrichment for %s is skipped.",
            tool_name,
        )
        return {}

    best_by_mol_idx: dict[str, dict] = {}
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path))
    except Exception as e:
        logger.warning("Could not read %s docking SDF %s: %s", tool_name, sdf_path, e)
        return {}

    for mol in supplier:
        if mol is None or not mol.HasProp("_Name"):
            continue

        mol_idx = str(mol.GetProp("_Name")).strip()
        if not mol_idx:
            continue

        affinity = _extract_float_prop(
            mol, ["minimizedAffinity", "affinity", "score", "docking_score"]
        )
        if affinity is None:
            continue

        if tool_name == DOCKING_TOOL_GNINA:
            record = {
                "gnina_affinity": affinity,
                "gnina_cnnscore": _extract_float_prop(mol, ["CNNscore"]),
                "gnina_cnnaffinity": _extract_float_prop(mol, ["CNNaffinity"]),
                "gnina_cnn_vs": _extract_float_prop(mol, ["CNN_VS"]),
            }
            key = "gnina_affinity"
        elif tool_name == DOCKING_TOOL_MATCHA:
            record = {"matcha_affinity": affinity}
            key = "matcha_affinity"
        else:
            record = {"smina_affinity": affinity}
            key = "smina_affinity"

        current = best_by_mol_idx.get(mol_idx)
        if current is None or record[key] < current[key]:
            best_by_mol_idx[mol_idx] = record

    return best_by_mol_idx


def _collect_docking_scores(
    base_path: Path, docking_cfg: dict | None = None
) -> pd.DataFrame:
    """Collect best docking scores per molecule across configured tools."""
    merged: dict[str, dict] = {}

    for tool in (DOCKING_TOOL_GNINA, DOCKING_TOOL_SMINA, DOCKING_TOOL_MATCHA):
        sdf_path = _find_docking_sdf(base_path, tool, docking_cfg=docking_cfg)
        if not sdf_path:
            continue

        tool_scores = _extract_best_tool_scores(sdf_path, tool)
        for mol_idx, score_map in tool_scores.items():
            entry = merged.setdefault(mol_idx, {})
            entry.update(score_map)

    rows = []
    for mol_idx, score_map in merged.items():
        row = {"mol_idx": mol_idx}
        for col in DOCKING_SCORE_COLUMNS:
            row[col] = score_map.get(col)
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=["mol_idx", *DOCKING_SCORE_COLUMNS])
    return pd.DataFrame(rows)


class PipelineStage:
    """Represents a single stage in the molecular analysis pipeline."""

    def __init__(self, name: str, config_key: str, directory: str):
        self.name = name
        self.config_key = config_key
        self.directory = directory
        self.enabled = False
        self.completed = False


class DataChecker:
    """Checks for existence of stage output data files."""

    # Mapping of stage names to their output file paths (relative to base_path)
    _STAGE_OUTPUT_PATHS = {
        DIR_MOL_PREP: Path(DIR_MOL_PREP) / FILE_FILTERED_MOLECULES,
        DIR_DESCRIPTORS_INITIAL: Path(DIR_DESCRIPTORS_INITIAL)
        / "filtered"
        / FILE_FILTERED_MOLECULES,
        DIR_STRUCT_FILTERS_POST: Path(DIR_STRUCT_FILTERS_POST)
        / FILE_FILTERED_MOLECULES,
        DIR_SYNTHESIS: Path(DIR_SYNTHESIS) / FILE_FILTERED_MOLECULES,
        DIR_DOCKING_FILTERS: Path(DIR_DOCKING_FILTERS) / FILE_FILTERED_MOLECULES,
        # Legacy paths
        "Descriptors": Path("Descriptors") / "passDescriptorsSMILES.csv",
        "StructFilters": Path("StructFilters") / "passStructFiltersSMILES.csv",
        "beforeDescriptors": Path("beforeDescriptors_StructFilters")
        / "passStructFiltersSMILES.csv",
        "Synthesis": Path("Synthesis") / "passSynthesisSMILES.csv",
    }

    def __init__(self, config: dict, progress_callback=None):
        self.config = config
        self.progress_callback = progress_callback
        self.base_path = Path(config[CONFIG_FOLDER_TO_SAVE])

    def check_stage_data(self, stage_name: str) -> bool:
        """Check if data exists for a given stage."""
        path = self._get_stage_output_path(stage_name.strip())
        return _file_exists_and_not_empty(path) if path else False

    def stage_has_molecules(self, stage_name: str) -> bool:
        """Check if a stage output contains at least one molecule row."""
        path = self._get_stage_output_path(stage_name.strip())
        if path is None:
            return False
        if path.suffix.lower() == ".csv":
            return _csv_has_data_rows(path)
        return _file_exists_and_not_empty(path)

    def _get_stage_output_path(self, stage_name: str) -> Path | None:
        """Get the expected output file path for a stage."""
        relative_path = self._STAGE_OUTPUT_PATHS.get(stage_name)
        return self.base_path / relative_path if relative_path else None


class PipelineStageRunner:
    """Executes individual pipeline stages and manages stage data flow."""

    # Local priority list for stage-based data checking (uses directory names)
    DATA_SOURCE_PRIORITY = [
        DIR_DOCKING_FILTERS,
        DIR_SYNTHESIS,
        DIR_STRUCT_FILTERS_POST,
        DIR_DESCRIPTORS_INITIAL,
        DIR_MOL_PREP,
    ]

    def __init__(self, config: dict, data_checker: DataChecker, progress_callback=None):
        self.config = config
        self.progress_callback = progress_callback
        self.data_checker = data_checker

    def find_latest_data_source(self) -> str | None:
        """Find the most recent stage with available output data.

        Uses centralized input path discovery from utils.input_paths module.
        """
        # Use centralized function for file discovery
        base_path = self.data_checker.base_path
        input_path = _find_input(base_path)
        if input_path:
            # Convert path back to stage directory name
            rel_path = input_path.relative_to(base_path)
            for source in self.DATA_SOURCE_PRIORITY:
                if source in str(rel_path):
                    logger.debug("Found data from stage: %s", source)
                    return source
            # Fallback: return None if path doesn't match known stages
            logger.debug("Found data at %s but doesn't match known stages", input_path)
            return None

        # Fallback to original logic if centralized function returns None
        for source in self.DATA_SOURCE_PRIORITY:
            if self.data_checker.check_stage_data(source):
                logger.debug("Found data from stage: %s", source)
                return source
        logger.debug("No processed data found from any stage")
        return None

    def run_mol_prep(
        self,
        data,
        subfolder: str | None = None,
        reporter: StageProgressReporter | None = None,
    ) -> bool:
        """Run Datamol-based Mol Prep stage."""
        try:
            config_path = self.config.get(CONFIG_MOL_PREP)
            if not config_path:
                logger.info("Mol Prep config not specified, skipping")
                return False

            cfg = load_config(config_path)
            if not cfg.get(CONFIG_RUN_KEY, False):
                logger.info("Mol Prep disabled in config")
                return False

            mol_prep_main(data, self.config, subfolder=subfolder, reporter=reporter)
            return True
        except InterruptedError:
            raise
        except Exception as exc:
            logger.error("Error running Mol Prep: %s", exc)
            return False

    def run_descriptors(
        self,
        data,
        subfolder: str | None = None,
        reporter: StageProgressReporter | None = None,
    ) -> bool:
        """Run molecular descriptors calculation."""
        try:
            config_descriptors = load_config(self.config[CONFIG_DESCRIPTORS])
            if not config_descriptors.get(CONFIG_RUN_KEY, False):
                logger.info("Descriptors calculation disabled in config")
                return False
            descriptors_main(data, self.config, subfolder=subfolder, reporter=reporter)
            return True
        except InterruptedError:
            raise
        except Exception as e:
            logger.error("Error running descriptors: %s", e)
            return False

    def run_structural_filters(
        self, stage_dir: str, reporter: StageProgressReporter | None = None
    ) -> bool:
        """Run structural filters on molecules."""
        try:
            config_struct_filters = load_config(self.config[CONFIG_STRUCT_FILTERS])
            if not config_struct_filters.get(CONFIG_RUN_KEY, False):
                logger.info("Structural filters disabled in config")
                return False

            structural_filters_main(self.config, stage_dir, reporter=reporter)
            return True
        except InterruptedError:
            raise
        except Exception as e:
            logger.error("Error running structural filters: %s", e)
            return False

    def run_synthesis(self, reporter: StageProgressReporter | None = None) -> bool:
        """Run synthesis analysis."""
        try:
            if not self._validate_synthesis_input():
                return False

            config_synthesis = load_config(self.config[CONFIG_SYNTHESIS])
            if not config_synthesis.get(CONFIG_RUN_KEY, False):
                logger.info("Synthesis disabled in config")
                return False

            synthesis_main(self.config, reporter=reporter)

            output_path = (
                self.data_checker.base_path / DIR_SYNTHESIS / FILE_FILTERED_MOLECULES
            )
            if not output_path.exists():
                logger.error("Synthesis finished but no output file detected")
                return False
            return True
        except InterruptedError:
            raise
        except Exception as e:
            logger.error("Error running synthesis: %s", e)
            return False

    def _validate_synthesis_input(self) -> bool:
        """Validate that input data exists for synthesis stage."""
        if self.find_latest_data_source():
            return True

        if self.config.get(OVERRIDE_SINGLE_STAGE) != STAGE_SYNTHESIS:
            logger.warning("No data available for synthesis, check provided path")
            return False

        base = self.data_checker.base_path
        sampled_path = base / "input" / FILE_SAMPLED_MOLECULES
        if not sampled_path.exists():
            sampled_path = base / FILE_SAMPLED_MOLECULES

        if sampled_path.exists():
            logger.info("Using %s", FILE_SAMPLED_MOLECULES)
            return True

        logger.warning("No data available for synthesis, check `config.yml` file")
        return False

    def _parse_docking_tools(self, tools_cfg) -> list[str]:
        """Parse the tools configuration to get a list of docking tools."""
        if isinstance(tools_cfg, str):
            tools_list = [t.strip().lower() for t in tools_cfg.split(",")]
        elif isinstance(tools_cfg, (list, tuple)):
            tools_list = [str(t).strip().lower() for t in tools_cfg]
        else:
            tools_list = [DOCKING_TOOL_BOTH]

        if not tools_list:
            return [DOCKING_TOOL_SMINA, DOCKING_TOOL_GNINA]

        selected_tools = []

        def _append(tool_name: str) -> None:
            if tool_name not in selected_tools:
                selected_tools.append(tool_name)

        for tool_name in tools_list:
            if tool_name == DOCKING_TOOL_ALL:
                _append(DOCKING_TOOL_SMINA)
                _append(DOCKING_TOOL_GNINA)
                _append(DOCKING_TOOL_MATCHA)
                continue
            if tool_name == DOCKING_TOOL_BOTH:
                _append(DOCKING_TOOL_SMINA)
                _append(DOCKING_TOOL_GNINA)
                continue
            if tool_name in {
                DOCKING_TOOL_SMINA,
                DOCKING_TOOL_GNINA,
                DOCKING_TOOL_MATCHA,
            }:
                _append(tool_name)

        return selected_tools or [DOCKING_TOOL_SMINA, DOCKING_TOOL_GNINA]

    def docking_results_present(self) -> bool:
        """Check if docking results exist for any configured tools."""
        try:
            cfg = load_config(self.config[CONFIG_DOCKING])
        except Exception:
            return False

        base_folder = self.data_checker.base_path.resolve()
        tools_list = self._parse_docking_tools(cfg.get(CONFIG_TOOLS, DOCKING_TOOL_BOTH))

        for tool in tools_list:
            if tool == DOCKING_TOOL_GNINA:
                gnina_sdf = (
                    self._get_gnina_output_dir(cfg, base_folder) / FILE_GNINA_OUTPUT
                )
                if _file_exists_and_not_empty(gnina_sdf):
                    return True
                logger.debug("GNINA results not found at %s", gnina_sdf)
            elif tool == DOCKING_TOOL_SMINA:
                smina_sdf = (
                    base_folder
                    / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_SMINA]
                    / FILE_SMINA_OUTPUT
                )
                if _file_exists_and_not_empty(smina_sdf):
                    return True
                logger.debug("SMINA results not found at %s", smina_sdf)
            elif tool == DOCKING_TOOL_MATCHA:
                matcha_sdf = (
                    base_folder
                    / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_MATCHA]
                    / FILE_MATCHA_OUTPUT
                )
                if _file_exists_and_not_empty(matcha_sdf):
                    return True
                logger.debug("MATCHA results not found at %s", matcha_sdf)
        return False

    def _get_gnina_output_dir(self, cfg: dict, base_folder: Path) -> Path:
        """Get the GNINA output directory from config."""
        gnina_config = cfg.get("gnina_config", {})
        cfg_out_dir = gnina_config.get("output_dir") or cfg.get("gnina_output_dir")
        if cfg_out_dir:
            out_dir = Path(cfg_out_dir)
            return out_dir if out_dir.is_absolute() else base_folder / out_dir
        return base_folder / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_GNINA]

    def run_docking(self, reporter: StageProgressReporter | None = None) -> bool:
        """Run molecular docking."""
        try:
            config_docking = load_config(self.config[CONFIG_DOCKING])
            if not config_docking.get(CONFIG_RUN_KEY, False):
                logger.info("Docking disabled in config")
                return False

            if not docking_main(self.config, reporter=reporter):
                return False

            if not self.docking_results_present():
                completed_empty = (
                    self.data_checker.base_path
                    / DIR_DOCKING
                    / FILE_DOCKING_COMPLETED_EMPTY
                )
                if _file_exists_and_not_empty(completed_empty):
                    logger.info(
                        "Docking completed with no valid ligands (%s)",
                        completed_empty,
                    )
                    return True
                logger.error(
                    "Docking finished but no results detected in output directories"
                )
                return False
            return True
        except InterruptedError:
            raise
        except Exception as e:
            logger.error("Error running docking: %s", e)
            return False

    def run_docking_filters(
        self, reporter: StageProgressReporter | None = None
    ) -> bool:
        """Run docking filters stage."""
        try:
            # Check if config exists
            config_path = self.config.get(CONFIG_DOCKING_FILTERS)
            if config_path is None:
                logger.info("Docking filters config not specified, skipping")
                return False

            config_filters = load_config(config_path)
            if not config_filters.get(CONFIG_RUN_KEY, False):
                logger.info("Docking filters disabled in config")
                return False

            input_sdf_cfg = config_filters.get("input_sdf")
            if input_sdf_cfg:
                input_sdf_path = Path(str(input_sdf_cfg))
                if not input_sdf_path.is_absolute():
                    input_sdf_path = (
                        self.data_checker.base_path.resolve() / input_sdf_path
                    )
                if _file_exists_and_not_empty(input_sdf_path):
                    logger.info(
                        "Using docking_filters.input_sdf: %s (skipping docking results presence guard)",
                        input_sdf_path,
                    )
                    result = docking_filters_main(self.config, reporter=reporter)
                    return result is not None and len(result) > 0

            # Check if docking results exist before running filters
            if not self.docking_results_present():
                logger.warning(
                    "No docking results found. Skipping docking filters. "
                    "Ensure docking stage completed successfully before running filters."
                )
                return False

            # Run docking filters
            result = docking_filters_main(self.config, reporter=reporter)
            return result is not None and len(result) > 0
        except InterruptedError:
            raise
        except Exception as e:
            logger.error("Error running docking filters: %s", e)
            return False


# ---------------------------------------------------------------------------
# MoleculeCounter -- extracted from MolecularAnalysisPipeline
# ---------------------------------------------------------------------------


class MoleculeCounter:
    """Count molecules in stage input/output CSV files."""

    # Stage name -> output CSV path (relative to base_path)
    _OUTPUT_PATHS: dict[str, tuple[str, ...]] = {
        STAGE_MOL_PREP: (DIR_MOL_PREP, FILE_FILTERED_MOLECULES),
        STAGE_DESCRIPTORS: (
            DIR_DESCRIPTORS_INITIAL,
            "filtered",
            FILE_FILTERED_MOLECULES,
        ),
        STAGE_STRUCT_FILTERS: (DIR_STRUCT_FILTERS_POST, FILE_FILTERED_MOLECULES),
        STAGE_SYNTHESIS: (DIR_SYNTHESIS, FILE_FILTERED_MOLECULES),
        STAGE_DOCKING_FILTERS: (DIR_DOCKING_FILTERS, FILE_FILTERED_MOLECULES),
        STAGE_FINAL_DESCRIPTORS: (
            DIR_DESCRIPTORS_FINAL,
            "filtered",
            FILE_FILTERED_MOLECULES,
        ),
    }

    def __init__(self, base_path: Path) -> None:
        self.base_path = base_path

    def count_csv_rows(self, path: Path) -> int | None:
        """Count data rows in a CSV file (excluding header)."""
        if not path.exists():
            return None
        try:
            with path.open(encoding="utf-8", errors="ignore", newline="") as handle:
                reader = csv.reader(handle)
                try:
                    next(reader)
                except StopIteration:
                    return 0
                return sum(1 for row in reader if any(cell.strip() for cell in row))
        except OSError:
            return None

    def resolve_input_count(
        self,
        stage_name: str,
        args: tuple,
        current_data,
        get_latest_data_fn=None,
    ) -> int | None:
        """Resolve molecule count entering a stage.

        Parameters
        ----------
        stage_name:
            Name of the pipeline stage.
        args:
            Positional arguments passed to the stage runner.
        current_data:
            The current pipeline DataFrame (may be None).
        get_latest_data_fn:
            Optional callable returning the latest data for final_descriptors.
        """
        for arg in args:
            if isinstance(arg, pd.DataFrame):
                return len(arg)

        if stage_name == STAGE_FINAL_DESCRIPTORS and get_latest_data_fn is not None:
            final_data = get_latest_data_fn(
                skip_descriptors=True, fallback_on_empty=False
            )
            if final_data is not None:
                return len(final_data)

        source = _find_input(self.base_path)
        if source is not None and source.suffix.lower() == ".csv":
            counted = self.count_csv_rows(source)
            if counted is not None:
                return counted

        if current_data is not None:
            return len(current_data)
        return None

    def resolve_output_count(
        self,
        stage_name: str,
        completed: bool,
        input_count: int | None,
    ) -> int | None:
        """Resolve molecule count after a stage completes."""
        path_parts = self._OUTPUT_PATHS.get(stage_name)
        if path_parts is not None:
            output_path = self.base_path.joinpath(*path_parts)
            counted = self.count_csv_rows(output_path)
            if counted is not None:
                return counted

        if stage_name == STAGE_DOCKING:
            return input_count

        if not completed:
            return None
        return None


# ---------------------------------------------------------------------------
# PipelineReporter -- extracted from MolecularAnalysisPipeline
# ---------------------------------------------------------------------------

# Skip-condition map: stage -> required upstream directory
_SKIP_CONDITIONS: dict[str, str] = {
    STAGE_STRUCT_FILTERS: DIR_DESCRIPTORS,
    STAGE_SYNTHESIS: DIR_STRUCT_FILTERS,
}


class PipelineReporter:
    """Log summaries, classify outcomes, save outputs, generate reports."""

    def __init__(
        self,
        base_path: Path,
        stages: list[PipelineStage],
        config: dict,
        counter: MoleculeCounter,
        data_checker: DataChecker,
        stage_timings: dict[str, float],
        stage_cpu_avg: dict[str, float],
        find_data_source_fn,
        build_data_path_fn,
    ) -> None:
        self.base_path = base_path
        self.stages = stages
        self.config = config
        self.counter = counter
        self.data_checker = data_checker
        self.stage_timings = stage_timings
        self.stage_cpu_avg = stage_cpu_avg
        self._find_data_source = find_data_source_fn
        self._build_data_path = build_data_path_fn

    # -- outcome classification -------------------------------------------

    def classify_outcome(self, stage: PipelineStage) -> str:
        """Classify a stage outcome as COMPLETED / FAILED / DISABLED / SKIPPED.

        Returns a plain label (no Rich markup).
        """
        if not stage.enabled:
            return "DISABLED"
        if stage.completed:
            return "COMPLETED"

        if stage.name in (STAGE_DOCKING, STAGE_DOCKING_FILTERS):
            source = _find_input(self.base_path)
            if source is None:
                return "SKIPPED"
            if source.suffix.lower() == ".csv" and not _csv_has_data_rows(source):
                return "SKIPPED"

        if stage.name == STAGE_FINAL_DESCRIPTORS:
            sources = [
                DIR_DOCKING_FILTERS,
                DIR_SYNTHESIS,
                DIR_STRUCT_FILTERS_POST,
                DIR_MOL_PREP,
            ]
            latest = self._find_data_source(sources)
            if not latest:
                return "SKIPPED"
            try:
                latest_df = pd.read_csv(self._build_data_path(latest))
            except Exception:
                return "FAILED"
            if len(latest_df) == 0:
                return "SKIPPED"

        required_data = _SKIP_CONDITIONS.get(stage.name)
        if required_data and not self.data_checker.stage_has_molecules(required_data):
            return "SKIPPED"

        return "FAILED"

    def stage_is_failed(self, stage: PipelineStage) -> bool:
        """Return True if an enabled stage is considered a failure."""
        outcome = self.classify_outcome(stage)
        return outcome == "FAILED"

    # -- summary logging --------------------------------------------------

    def log_summary(self) -> None:
        """Log a summary of pipeline execution status and timings."""
        _log_stage_header("Pipeline Execution Summary")

        for stage in self.stages:
            status = self._format_stage_status(stage)
            timing = self.stage_timings.get(stage.name)
            avg_cpu = self.stage_cpu_avg.get(stage.name)
            if timing is not None and avg_cpu is not None:
                logger.info(
                    "%s: %s (%.1fs, avg CPU %.1f%%)",
                    stage.name,
                    status,
                    timing,
                    avg_cpu,
                )
            elif timing is not None:
                logger.info("%s: %s (%.1fs)", stage.name, status, timing)
            else:
                logger.info("%s: %s", stage.name, status)

        if self.stage_timings:
            total_time = sum(self.stage_timings.values())
            logger.info("Total pipeline time: %.1f seconds", total_time)

    def log_molecule_summary(self, initial_count: int, final_count: int) -> None:
        """Log molecule count summary."""
        logger.info("---------> Molecule Count Summary")
        logger.info("Molecules before filtration: %d", initial_count)
        logger.info("Molecules remaining after all stages: %d", final_count)
        if initial_count > 0:
            retention = 100 * final_count / initial_count
            logger.info("Retention rate: %.2f%%", retention)

    # -- output saving ----------------------------------------------------

    def save_final_output(self, final_data, final_count: int) -> None:
        """Save final molecules to output directory."""
        final_output_path = self.base_path / DIR_OUTPUT / FILE_FINAL_MOLECULES
        final_output_path.parent.mkdir(parents=True, exist_ok=True)

        id_cols = ["smiles", "model_name", "mol_idx"]
        stable_empty_cols = [*id_cols, *DOCKING_SCORE_COLUMNS]

        if final_data is None or len(final_data) == 0:
            pd.DataFrame(columns=stable_empty_cols).to_csv(
                final_output_path, index=False
            )
            logger.info("Saved 0 final molecules to %s", final_output_path)
            return

        output_df = final_data.copy()
        if "mol_idx" in output_df.columns:
            docking_cfg = None
            docking_cfg_path = self.config.get(CONFIG_DOCKING)
            if docking_cfg_path:
                try:
                    docking_cfg = load_config(docking_cfg_path)
                except Exception as exc:
                    logger.warning(
                        "Could not load docking config for score collection (%s): %s",
                        docking_cfg_path,
                        exc,
                    )
            docking_scores = _collect_docking_scores(
                self.base_path, docking_cfg=docking_cfg
            )
            if not docking_scores.empty:
                output_df["_join_mol_idx"] = output_df["mol_idx"].astype(str)
                docking_scores = docking_scores.rename(
                    columns={"mol_idx": "_join_mol_idx"}
                )
                output_df = output_df.merge(
                    docking_scores, on="_join_mol_idx", how="left"
                )
                output_df = output_df.drop(columns=["_join_mol_idx"])
        else:
            logger.warning(
                "Final output is missing 'mol_idx'; docking score enrichment was skipped."
            )

        for col in DOCKING_SCORE_COLUMNS:
            if col not in output_df.columns:
                output_df[col] = pd.NA

        output_df.to_csv(final_output_path, index=False)
        logger.info("Saved %d final molecules to %s", final_count, final_output_path)

    # -- HTML report ------------------------------------------------------

    def generate_html_report(
        self,
        initial_count: int,
        final_count: int,
        progress_reporter: StageProgressReporter | None = None,
    ) -> tuple[bool, float]:
        """Generate HTML report for the pipeline run.

        Returns:
            Tuple ``(ok, elapsed_seconds)``.
        """
        started_at = time.perf_counter()
        try:
            logger.info("Generating HTML report (this may take a while)...")
            report_generator = ReportGenerator(
                base_path=self.base_path,
                stages=self.stages,
                config=self.config,
                initial_count=initial_count,
                final_count=final_count,
                progress_callback=(
                    progress_reporter.progress
                    if progress_reporter is not None
                    else None
                ),
            )
            report_path = report_generator.generate()
            logger.info("HTML report generated: %s", report_path)
            return True, time.perf_counter() - started_at
        except Exception as e:
            logger.warning("Failed to generate HTML report: %s", e)
            return False, time.perf_counter() - started_at

    # -- internal helpers -------------------------------------------------

    def _format_stage_status(self, stage: PipelineStage) -> str:
        """Get display status string for a stage (with optional Rich markup)."""
        outcome = self.classify_outcome(stage)

        if _plain_output_enabled():
            if outcome == "SKIPPED":
                return "SKIPPED (no molecules)"
            return outcome

        _RICH_MAP = {
            "DISABLED": "DISABLED",
            "COMPLETED": "[bold]\u2713 COMPLETED[/bold]",
            "SKIPPED": "[dim]\u27c2 SKIPPED (no molecules)[/dim]",
            "FAILED": "[bold]\u2717 FAILED[/bold]",
        }
        return _RICH_MAP.get(outcome, outcome)


# ---------------------------------------------------------------------------
# MolecularAnalysisPipeline -- now focused on stage orchestration
# ---------------------------------------------------------------------------


class MolecularAnalysisPipeline:
    """Orchestrates the execution of the molecular analysis pipeline."""

    # Stage definitions: (name, config_key, directory)
    _STAGE_DEFINITIONS = [
        (STAGE_MOL_PREP, CONFIG_MOL_PREP, DIR_MOL_PREP),
        (STAGE_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_DESCRIPTORS),
        (STAGE_STRUCT_FILTERS, CONFIG_STRUCT_FILTERS, DIR_STRUCT_FILTERS),
        (STAGE_SYNTHESIS, CONFIG_SYNTHESIS, DIR_SYNTHESIS),
        (STAGE_DOCKING, CONFIG_DOCKING, DIR_DOCKING),
        (STAGE_DOCKING_FILTERS, CONFIG_DOCKING_FILTERS, DIR_DOCKING_FILTERS),
        (STAGE_FINAL_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_FINAL_DESCRIPTORS),
    ]

    # Stage display labels for logging
    _STAGE_LABELS = {
        STAGE_MOL_PREP: "Stage 1: Mol Prep (Datamol)",
        STAGE_DESCRIPTORS: "Stage 2: Molecular Descriptors",
        STAGE_STRUCT_FILTERS: "Stage 3: Structural Filters",
        STAGE_SYNTHESIS: "Stage 4: Synthesis Analysis",
        STAGE_DOCKING: "Stage 5: Molecular Docking",
        STAGE_DOCKING_FILTERS: "Stage 6: Docking Filters",
        STAGE_FINAL_DESCRIPTORS: "Stage 7': Final Descriptors Calculation",
    }

    def __init__(self, config: dict, progress_callback=None):
        self.config = config
        self.progress_callback = progress_callback
        self.data_checker = DataChecker(config)
        self.stage_runner = PipelineStageRunner(config, self.data_checker)
        self.current_data = None
        self.stages = [PipelineStage(*defn) for defn in self._STAGE_DEFINITIONS]
        self._stage_by_name = {stage.name: stage for stage in self.stages}
        self.stage_timings: dict[str, float] = {}  # Stage name -> elapsed seconds
        self.stage_cpu_avg: dict[str, float] = {}  # Stage name -> avg CPU percentage

        self.counter = MoleculeCounter(self.data_checker.base_path)
        self.reporter = PipelineReporter(
            base_path=self.data_checker.base_path,
            stages=self.stages,
            config=config,
            counter=self.counter,
            data_checker=self.data_checker,
            stage_timings=self.stage_timings,
            stage_cpu_avg=self.stage_cpu_avg,
            find_data_source_fn=self._find_data_source,
            build_data_path_fn=self._build_data_path,
        )

        self._initialize_stages()

    def _initialize_stages(self) -> None:
        """Initialize stages by loading their configs and determining enabled status."""
        single_stage_override = self.config.get(OVERRIDE_SINGLE_STAGE)
        selected_stage_names: set[str] | None = None
        selection_override = self.config.get(OVERRIDE_STAGE_SELECTION)

        if isinstance(selection_override, str):
            selected_stage_names = {
                part.strip() for part in selection_override.split(",") if part.strip()
            }
        elif isinstance(selection_override, (list, tuple, set)):
            selected_stage_names = {str(part).strip() for part in selection_override}
            selected_stage_names = {part for part in selected_stage_names if part}

        if selected_stage_names:
            if len(selected_stage_names) == 1:
                single_stage_override = next(iter(selected_stage_names))
                logger.info(
                    "Single stage mode: enabling only %s", single_stage_override
                )
            else:
                ordered = [
                    stage_name
                    for stage_name, _, _ in self._STAGE_DEFINITIONS
                    if stage_name in selected_stage_names
                ]
                logger.info(
                    "Stage selection mode: enabling only %s", ", ".join(ordered)
                )

        for stage in self.stages:
            try:
                stage_config = load_config(self.config[stage.config_key])
                stage.enabled = stage_config.get(CONFIG_RUN_KEY, False)

                if selected_stage_names is not None:
                    stage.enabled = stage.name in selected_stage_names
                elif single_stage_override:
                    stage.enabled = stage.name == single_stage_override
                    if stage.enabled:
                        logger.info("Single stage mode: enabling only %s", stage.name)

                logger.debug(
                    "Stage %s: %s",
                    stage.name,
                    "Enabled" if stage.enabled else "Disabled",
                )
            except Exception as e:
                logger.warning("Could not load config for %s: %s", stage.name, e)
                stage.enabled = False

        # In stage-selection mode, always run Mol Prep first (if enabled in config),
        # so downstream stages operate on standardized molecules.
        selection_requires_mol_prep = False
        if selected_stage_names:
            selection_requires_mol_prep = any(
                stage_name != STAGE_MOL_PREP for stage_name in selected_stage_names
            )
        elif single_stage_override:
            selection_requires_mol_prep = single_stage_override != STAGE_MOL_PREP

        if selection_requires_mol_prep:
            mol_prep = self._stage_by_name.get(STAGE_MOL_PREP)
            if mol_prep is not None:
                try:
                    cfg_path = self.config.get(CONFIG_MOL_PREP)
                    if cfg_path:
                        cfg = load_config(cfg_path)
                        if cfg.get(CONFIG_RUN_KEY, False):
                            mol_prep.enabled = True
                            logger.info(
                                "Single stage mode: also enabling %s", STAGE_MOL_PREP
                            )
                except Exception:
                    pass

    def _cancel_requested(self) -> bool:
        """Return True when the progress callback exposes a cancellation signal."""
        if not self.progress_callback:
            return False

        cancel_probe = getattr(self.progress_callback, "is_cancelled", None)
        if cancel_probe is None:
            return False

        if callable(cancel_probe):
            try:
                return bool(cancel_probe())
            except Exception:
                return False

        return bool(cancel_probe)

    def _raise_if_cancel_requested(self, boundary: str) -> None:
        """Stop pipeline execution when cancellation is requested."""
        if self._cancel_requested():
            logger.info("Pipeline cancellation requested at %s", boundary)
            raise InterruptedError("Pipeline cancelled")

    def _emit_progress_event(self, event: dict) -> None:
        """Emit a structured progress event to the configured callback.

        Supports both:
          - New signature: ``callback(event: dict)``
          - Legacy signature: ``callback(stage_name: str, current: int, total: int)``
        """
        if not self.progress_callback:
            return

        try:
            self.progress_callback(event)
            return
        except TypeError:
            # Legacy callback: (stage_name, current, total)
            try:
                event_type = event.get("type")
                stage = event.get("stage")
                if not stage:
                    return

                if event_type == "stage_progress":
                    self.progress_callback(
                        stage, int(event.get("current", 0)), int(event.get("total", 0))
                    )
                elif event_type == "stage_start":
                    self.progress_callback(
                        stage,
                        int(event.get("stage_index", 0)),
                        int(event.get("total_stages", 0)),
                    )
                elif event_type == "stage_complete":
                    total = int(event.get("total_stages", 0))
                    self.progress_callback(stage, total, total)
            except Exception:
                return

    # -- delegated counting (backwards compat) ----------------------------

    def _count_csv_rows(self, path: Path) -> int | None:
        """Count data rows in a CSV file (excluding header)."""
        return self.counter.count_csv_rows(path)

    def _resolve_stage_input_count(self, stage_name: str, args: tuple) -> int | None:
        """Resolve molecule count entering a stage."""
        return self.counter.resolve_input_count(
            stage_name, args, self.current_data, self.get_latest_data
        )

    def _resolve_stage_output_count(
        self,
        stage_name: str,
        completed: bool,
        input_count: int | None,
    ) -> int | None:
        """Resolve molecule count after a stage completes."""
        return self.counter.resolve_output_count(stage_name, completed, input_count)

    # -- data loading -----------------------------------------------------

    def get_latest_data(
        self, skip_descriptors: bool = False, fallback_on_empty: bool = True
    ):
        """Load the most recent stage output data."""
        priority = self.stage_runner.DATA_SOURCE_PRIORITY.copy()
        if skip_descriptors:
            priority = [p for p in priority if p != DIR_DESCRIPTORS]

        latest_source = self._find_data_source(priority)
        if not latest_source:
            return self.current_data

        try:
            path = self._build_data_path(latest_source)
            data = pd.read_csv(path)

            if len(data) == 0 and fallback_on_empty:
                return self._try_fallback_sources(priority, latest_source, data)
            return data
        except Exception as e:
            logger.warning("Could not load data from %s: %s", latest_source, e)
            return self.current_data

    def _find_data_source(self, priority: list[str]) -> str | None:
        """Find first available data source from priority list."""
        for source in priority:
            if self.data_checker.check_stage_data(source):
                logger.debug("Found data from stage: %s", source)
                return source
        return None

    def _build_data_path(self, source: str) -> Path:
        """Build the path to stage output data."""
        base = self.data_checker.base_path
        if source.startswith("stages/"):
            if "descriptors" in source:
                return base / source / "filtered" / FILE_FILTERED_MOLECULES
            return base / source / FILE_FILTERED_MOLECULES
        return base / source / FILE_FILTERED_MOLECULES

    def _try_fallback_sources(
        self, priority: list[str], current_source: str, empty_data
    ):
        """Try fallback sources when current source is empty."""
        logger.warning("No molecules in %s, trying previous step...", current_source)
        current_idx = priority.index(current_source)

        for next_source in priority[current_idx + 1 :]:
            if self.data_checker.check_stage_data(next_source):
                path = self._build_data_path(next_source)
                try:
                    data = pd.read_csv(path)
                    if len(data) > 0:
                        logger.info(
                            "Loaded data from %s: %d molecules (previous step)",
                            next_source,
                            len(data),
                        )
                        return data
                except Exception:
                    continue

        logger.warning("All checked data sources are empty")
        return empty_data

    # -- pipeline execution -----------------------------------------------

    def run_pipeline(self, data) -> bool:
        """Execute the full molecular analysis pipeline."""
        self.current_data = data
        success_count = 0
        total_enabled = sum(1 for s in self.stages if s.enabled)
        last_completed_stage: str | None = None

        steps = [
            (STAGE_MOL_PREP, lambda: self._run_mol_prep(data)),
            (STAGE_DESCRIPTORS, lambda: self._run_descriptors(data)),
            (STAGE_STRUCT_FILTERS, self._run_post_descriptors_filters),
            (STAGE_SYNTHESIS, self._run_synthesis),
            (STAGE_DOCKING, self._run_docking),
            (STAGE_DOCKING_FILTERS, self._run_docking_filters),
            (STAGE_FINAL_DESCRIPTORS, self._run_final_descriptors),
        ]

        for name, run_step in steps:
            self._raise_if_cancel_requested(f"stage boundary before '{name}'")
            completed, early_exit = run_step()
            self._raise_if_cancel_requested(f"stage boundary after '{name}'")
            if completed:
                self._stage_by_name[name].completed = True
                success_count += 1
                last_completed_stage = name
            if early_exit:
                return self._finalize_pipeline(
                    data,
                    success_count,
                    total_enabled,
                    final_stage_name=last_completed_stage,
                )

        logger.info(
            "Pipeline completed: %d/%d stages successful", success_count, total_enabled
        )
        self._raise_if_cancel_requested("pipeline finalization")
        return self._finalize_pipeline(data, success_count, total_enabled)

    def _run_stage(
        self, stage_name: str, runner_func, *args, on_failure=None
    ) -> tuple[bool, bool]:
        """Run a pipeline stage with standard logging. Returns (completed, early_exit).

        Args:
            stage_name: Name of the stage to run.
            runner_func: Callable that performs the stage work. Must return bool.
            *args: Positional arguments passed to runner_func.
            on_failure: Optional callback invoked when runner_func returns False.
                Must return (completed, early_exit) tuple.
        """
        stage = self._stage_by_name[stage_name]
        if not stage.enabled:
            return False, False

        enabled_stages = [s.name for s in self.stages if s.enabled]
        current_idx = (
            enabled_stages.index(stage.name) if stage.name in enabled_stages else 0
        )
        reporter = StageProgressReporter(
            emit_event=self._emit_progress_event,
            stage=stage.name,
            stage_index=current_idx + 1,
            total_stages=len(enabled_stages),
        )
        input_count = self._resolve_stage_input_count(stage.name, args)
        start_message = _build_stage_start_message(stage.name, input_count)
        reporter.start(message=start_message, molecules_in=input_count)

        _log_stage_header(self._STAGE_LABELS[stage.name])
        logger.info(start_message)
        start_time = time.perf_counter()
        cpu_start = time.process_time()
        if _accepts_reporter(runner_func):
            completed = runner_func(*args, reporter=reporter)
        else:
            completed = runner_func(*args)
        elapsed = time.perf_counter() - start_time
        cpu_elapsed = max(0.0, time.process_time() - cpu_start)
        avg_cpu_percent = (100.0 * cpu_elapsed / elapsed) if elapsed > 0 else None
        output_count = self._resolve_stage_output_count(
            stage.name, bool(completed), input_count
        )
        self.stage_timings[stage.name] = elapsed
        if avg_cpu_percent is not None:
            self.stage_cpu_avg[stage.name] = avg_cpu_percent
        output_estimated = stage.name == STAGE_DOCKING and output_count is not None
        complete_message = _build_stage_complete_message(
            stage_name=stage.name,
            ok=bool(completed),
            molecules_in=input_count,
            molecules_out=output_count,
            elapsed_seconds=elapsed,
            avg_cpu_percent=avg_cpu_percent,
            output_estimated=output_estimated,
        )
        logger.info(complete_message)
        reporter.complete(
            ok=bool(completed),
            message=complete_message,
            molecules_in=input_count,
            molecules_out=output_count,
            elapsed_seconds=elapsed,
        )

        if completed:
            return True, False
        if on_failure:
            return on_failure()
        return False, False

    def _run_mol_prep(self, data) -> tuple[bool, bool]:
        """Run MolPrep stage (Datamol standardization + strict filtering)."""

        def _on_failure():
            logger.info("MolPrep failed; ending pipeline early.")
            return False, True

        completed, early_exit = self._run_stage(
            STAGE_MOL_PREP, self.stage_runner.run_mol_prep, data, on_failure=_on_failure
        )
        if not completed or early_exit:
            return completed, early_exit

        out_path = self.data_checker.base_path / DIR_MOL_PREP / FILE_FILTERED_MOLECULES
        if out_path.exists():
            try:
                df_out = pd.read_csv(out_path)
                self.current_data = df_out
                if len(df_out) == 0:
                    logger.info(
                        "No molecules left after MolPrep; ending pipeline early."
                    )
                    return True, True
            except Exception as e:
                logger.warning("Could not load MolPrep output (%s): %s", out_path, e)

        return True, False

    def _run_descriptors(self, data) -> tuple[bool, bool]:
        """Run descriptors calculation stage."""
        descriptors_input = data

        # If MolPrep is enabled, use its output as the input to descriptors.
        mol_prep_stage = self._stage_by_name.get(STAGE_MOL_PREP)
        if mol_prep_stage is not None and mol_prep_stage.enabled:
            prep_path = (
                self.data_checker.base_path / DIR_MOL_PREP / FILE_FILTERED_MOLECULES
            )
            if prep_path.exists():
                try:
                    prep_df = pd.read_csv(prep_path)
                    if len(prep_df) > 0:
                        descriptors_input = prep_df
                        self.current_data = prep_df
                except Exception as e:
                    logger.warning(
                        "Could not load MolPrep output (%s): %s",
                        prep_path,
                        e,
                    )

        return self._run_stage(
            STAGE_DESCRIPTORS, self.stage_runner.run_descriptors, descriptors_input
        )

    def _run_post_descriptors_filters(self) -> tuple[bool, bool]:
        """Run post-descriptors structural filters stage."""

        def _on_failure():
            logger.info(
                "Structural filters did not complete successfully; ending pipeline early."
            )
            return False, True

        return self._run_stage(
            STAGE_STRUCT_FILTERS,
            self.stage_runner.run_structural_filters,
            DIR_STRUCT_FILTERS_POST,
            on_failure=_on_failure,
        )

    def _run_synthesis(self) -> tuple[bool, bool]:
        """Run synthesis analysis stage."""
        return self._run_stage(
            STAGE_SYNTHESIS,
            self.stage_runner.run_synthesis,
            on_failure=self._handle_synthesis_failure,
        )

    def _handle_synthesis_failure(self) -> tuple[bool, bool]:
        """Handle synthesis stage failure. Returns (completed, early_exit)."""
        output_path = (
            self.data_checker.base_path / DIR_SYNTHESIS / FILE_FILTERED_MOLECULES
        )
        if not output_path.exists():
            logger.error(
                "Synthesis stage failed (no output file created). Check logs above for error details."
            )
            logger.info("Continuing pipeline without synthesis results...")
            return False, False

        try:
            df_check = pd.read_csv(output_path)
            if len(df_check) == 0:
                logger.info("No molecules left after synthesis; ending pipeline early.")
                return False, True
            logger.warning(
                "Synthesis stage failed but output file exists. Continuing with available molecules."
            )
        except Exception:
            logger.warning("Synthesis stage failed. Check logs for details.")

        return False, False

    def _run_docking(self) -> tuple[bool, bool]:
        """Run molecular docking stage."""
        source = _find_input(self.data_checker.base_path)
        if source is None:
            logger.info("No docking input found; skipping docking")
            return False, False
        if source.suffix.lower() == ".csv" and not _csv_has_data_rows(source):
            logger.info("No molecules available for docking; skipping docking")
            return False, False
        return self._run_stage(STAGE_DOCKING, self.stage_runner.run_docking)

    def _run_docking_filters(self) -> tuple[bool, bool]:
        """Run docking filters stage."""
        source = _find_input(self.data_checker.base_path)
        if source is None:
            logger.info("No docking input found; skipping docking filters")
            return False, False
        if source.suffix.lower() == ".csv" and not _csv_has_data_rows(source):
            logger.info("No molecules available for docking; skipping docking filters")
            return False, False
        return self._run_stage(
            STAGE_DOCKING_FILTERS, self.stage_runner.run_docking_filters
        )

    def _run_final_descriptors(self) -> tuple[bool, bool]:
        """Run final descriptors calculation stage."""

        def _run_final_desc():
            final_data = self.get_latest_data(
                skip_descriptors=True, fallback_on_empty=False
            )
            if final_data is None:
                logger.info("No data source found for final descriptors (skipping)")
                return False
            if len(final_data) == 0:
                logger.info(
                    "No molecules from previous steps; skipping final descriptors"
                )
                return False
            return self.stage_runner.run_descriptors(
                final_data, subfolder=DIR_FINAL_DESCRIPTORS
            )

        return self._run_stage(STAGE_FINAL_DESCRIPTORS, _run_final_desc)

    # -- finalization (delegates to reporter) ------------------------------

    def _finalize_pipeline(
        self,
        data,
        success_count: int,
        total_enabled: int,
        final_stage_name: str | None = None,
    ) -> bool:
        """Finalize pipeline execution with summary and output."""
        self.reporter.log_summary()

        initial_count = len(data)
        final_data = self._load_stage_output(final_stage_name)
        if final_data is None:
            final_data = self.get_latest_data(
                skip_descriptors=True, fallback_on_empty=False
            )
        final_count = len(final_data) if final_data is not None else 0

        self.reporter.log_molecule_summary(initial_count, final_count)
        self.reporter.save_final_output(final_data, final_count)
        _generate_structure_readme(
            self.data_checker.base_path,
            self.stages,
            initial_count,
            final_count,
            stage_timings=self.stage_timings,
            config=self.config,
        )

        report_progress: StageProgressReporter | None = None
        if self.progress_callback is not None:
            enabled_stage_count = sum(1 for s in self.stages if s.enabled)
            report_progress = StageProgressReporter(
                emit_event=self._emit_progress_event,
                stage=STAGE_REPORT,
                stage_index=enabled_stage_count + 1,
                total_stages=enabled_stage_count + 1,
            )
            report_progress.start(message="Generating HTML report")

        report_ok, report_elapsed = self.reporter.generate_html_report(
            initial_count,
            final_count,
            progress_reporter=report_progress,
        )
        if report_progress is not None:
            report_progress.complete(
                ok=report_ok,
                message=f"{report_elapsed:.1f}s",
                elapsed_seconds=report_elapsed,
            )

        return not any(self.reporter.stage_is_failed(s) for s in self.stages)

    def _load_stage_output(self, stage_name: str | None):
        """Load the explicit output for a completed stage, if one exists."""
        if not stage_name:
            return None

        stage = self._stage_by_name.get(stage_name)
        if stage is None:
            return None

        path = self.data_checker._get_stage_output_path(stage.directory)
        if path is None or not path.exists():
            return None

        try:
            return pd.read_csv(path)
        except Exception as e:
            logger.warning("Could not load %s output (%s): %s", stage.name, path, e)
            return None

    # -- backwards compatibility delegates --------------------------------
    # Tests call these methods directly on the pipeline instance.

    def _save_final_output(self, final_data, final_count: int) -> None:
        """Save final molecules to output directory (delegates to reporter)."""
        self.reporter.save_final_output(final_data, final_count)

    def _log_pipeline_summary(self) -> None:
        """Log pipeline summary (delegates to reporter)."""
        self.reporter.log_summary()

    def _get_stage_status(self, stage: PipelineStage) -> str:
        """Get display status for a stage (delegates to reporter)."""
        return self.reporter._format_stage_status(stage)

    def _stage_is_failed(self, stage: PipelineStage) -> bool:
        """Check if a stage is considered failed (delegates to reporter)."""
        return self.reporter.stage_is_failed(stage)

    def _generate_html_report(self, initial_count: int, final_count: int) -> None:
        """Generate HTML report (delegates to reporter)."""
        self.reporter.generate_html_report(initial_count, final_count)

    def _log_molecule_summary(self, initial_count: int, final_count: int) -> None:
        """Log molecule count summary (delegates to reporter)."""
        self.reporter.log_molecule_summary(initial_count, final_count)


# ---------------------------------------------------------------------------
# Config snapshot
# ---------------------------------------------------------------------------


def _save_config_snapshot(config: dict) -> None:
    """Save a snapshot of configuration files for provenance."""
    try:
        base_path = Path(config[CONFIG_FOLDER_TO_SAVE])
        dest_dir = base_path / DIR_CONFIGS
        dest_dir.mkdir(parents=True, exist_ok=True)

        master_config_path = dest_dir / FILE_MASTER_CONFIG
        with open(master_config_path, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)

        for key in config:
            if not key.startswith("config_"):
                continue
            path_str = config.get(key)
            if not path_str:
                continue
            try:
                src_path = Path(path_str)
                if src_path.exists():
                    shutil.copyfile(src_path, dest_dir / src_path.name)
            except OSError as copy_err:
                logger.warning("Could not copy config file for %s: %s", key, copy_err)

        logger.info("Saved run config snapshot to: %s", dest_dir)
    except Exception as snapshot_err:
        logger.warning("Config snapshot failed: %s", snapshot_err)


# ---------------------------------------------------------------------------
# RUN_INFO.md generation (Task 4: template + data-prep separated)
# ---------------------------------------------------------------------------

# Stage descriptions for README generation
_STAGE_DESCRIPTIONS = {
    STAGE_MOL_PREP: "Datamol molecule preparation (standardization + strict filtering)",
    STAGE_DESCRIPTORS: "Physicochemical descriptors",
    STAGE_STRUCT_FILTERS: "Structural filters",
    STAGE_SYNTHESIS: "Retrosynthesis analysis",
    STAGE_DOCKING: "Molecular docking",
    STAGE_DOCKING_FILTERS: "Post-docking pose & interaction filters",
    STAGE_FINAL_DESCRIPTORS: "Final descriptor calculation",
}

# Directory-level tree templates per stage (static parts)
_STAGE_TREE_TEMPLATES: dict[str, list[str]] = {
    STAGE_MOL_PREP: [
        "|   +-- filtered_molecules.csv     Standardized molecules (passed)",
        "|   +-- failed_molecules.csv       Rejected molecules with reasons",
        "|   +-- metrics.csv                Summary counts and failure breakdown",
        "|   +-- duplicates_removed.csv     Removed duplicates (optional)",
    ],
    STAGE_DESCRIPTORS: [
        "|   +-- metrics/",
        "|   |   +-- descriptors_all.csv    All computed descriptors",
        "|   |   +-- skipped_molecules.csv  Failed to parse",
        "|   +-- filtered/",
        "|   |   +-- filtered_molecules.csv Passed descriptor thresholds",
        "|   |   +-- failed_molecules.csv   Failed descriptor thresholds",
        "|   |   +-- descriptors_passed.csv Detailed metrics (passed)",
        "|   |   +-- descriptors_failed.csv Detailed metrics (failed)",
        "|   |   +-- pass_flags.csv         Pass/fail flags per descriptor",
        "|   +-- plots/",
        "|       +-- descriptors_distribution.png",
    ],
    STAGE_STRUCT_FILTERS: [
        "|   +-- {filter_name}/             Per-filter results",
        "|   |   +-- metrics.csv",
        "|   |   +-- extended.csv",
        "|   |   +-- filtered_molecules.csv",
        "|   +-- filtered_molecules.csv     Combined passed molecules",
        "|   +-- failed_molecules.csv       Combined failed molecules",
        "|   +-- plots/",
        "|       +-- molecule_counts_comparison.png",
        "|       +-- restriction_ratios_comparison.png",
    ],
    STAGE_SYNTHESIS: [
        "|   +-- synthesis_scores.csv       SA Score, RA Score, SYBA",
        "|   +-- synthesis_extended.csv     With retrosynthesis results",
        "|   +-- filtered_molecules.csv     Synthesizable molecules",
        "|   +-- input_smiles.smi           Input for AiZynthFinder",
        "|   +-- retrosynthesis_results.json",
    ],
    STAGE_DOCKING_FILTERS: [
        "|   +-- metrics.csv                Per-molecule filter metrics and pass flags",
        "|   +-- filtered_molecules.csv     Passing molecules (single best pose)",
        "|   +-- filtered_poses.csv         Passing poses with full metrics (1 per molecule)",
        "|   +-- filtered_poses.sdf         Filtered poses (SDF, 1 per molecule)",
        "|   +-- interaction_events.csv     Per-interaction event records",
        "|   +-- interaction_residue_summary.csv  Residue-level interaction counts",
        "|   +-- interaction_type_summary.csv     Interaction-type distribution",
        "|   +-- interaction_matrix.csv     Residue x interaction-type matrix",
        "|   +-- interaction_report_meta.json     Interaction reporting metadata",
    ],
    STAGE_FINAL_DESCRIPTORS: [
        "|   +-- metrics/",
        "|   |   +-- descriptors_all.csv    All computed descriptors",
        "|   |   +-- skipped_molecules.csv  Failed to parse",
        "|   +-- filtered/",
        "|   |   +-- filtered_molecules.csv Passed descriptor thresholds",
        "|   |   +-- failed_molecules.csv   Failed descriptor thresholds",
        "|   |   +-- descriptors_passed.csv Detailed metrics (passed)",
        "|   |   +-- descriptors_failed.csv Detailed metrics (failed)",
        "|   |   +-- pass_flags.csv         Pass/fail flags per descriptor",
        "|   +-- plots/",
        "|       +-- descriptors_distribution.png",
    ],
}


def _build_docking_tree(base_path: Path, config: dict | None) -> list[str]:
    """Build docking stage directory tree based on actual tools used."""
    docking_dir = base_path / DIR_DOCKING
    has_smina = (docking_dir / "smina").is_dir()
    has_gnina = (docking_dir / "gnina").is_dir()
    has_matcha = (docking_dir / "matcha").is_dir()

    # Fall back to config if directories don't exist yet
    if not has_smina and not has_gnina and not has_matcha and config:
        docking_cfg_path = config.get(CONFIG_DOCKING)
        if docking_cfg_path:
            try:
                docking_cfg = load_config(docking_cfg_path)
                raw_tools = docking_cfg.get(CONFIG_TOOLS, "")
                if isinstance(raw_tools, str):
                    tools = {
                        tool.strip().lower()
                        for tool in raw_tools.split(",")
                        if tool.strip()
                    }
                elif isinstance(raw_tools, (list, tuple)):
                    tools = {
                        str(tool).strip().lower()
                        for tool in raw_tools
                        if str(tool).strip()
                    }
                else:
                    tools = set()
                has_smina = (
                    DOCKING_TOOL_SMINA in tools
                    or DOCKING_TOOL_BOTH in tools
                    or DOCKING_TOOL_ALL in tools
                )
                has_gnina = (
                    DOCKING_TOOL_GNINA in tools
                    or DOCKING_TOOL_BOTH in tools
                    or DOCKING_TOOL_ALL in tools
                )
                has_matcha = DOCKING_TOOL_MATCHA in tools or DOCKING_TOOL_ALL in tools
            except Exception:
                has_smina = has_gnina = has_matcha = True

    lines = [
        "|   +-- ligands.csv                Prepared ligands",
        "|   +-- job_meta.json              Job metadata",
    ]
    if has_smina:
        lines.append("|   +-- smina/")
        lines.append(
            "|   |   +-- smina_out.sdf          Aggregated SMINA results (1 pose/molecule)"
        )
    if has_gnina:
        lines.append("|   +-- gnina/")
        lines.append(
            "|   |   +-- gnina_out.sdf          Aggregated GNINA results (1 pose/molecule)"
        )
    if has_matcha:
        lines.append("|   +-- matcha/")
        lines.append(
            "|   |   +-- matcha_out.sdf         Aggregated Matcha results (1 pose/molecule)"
        )

    # _workdir subtree
    lines.append("|   +-- _workdir/                  Intermediate files")
    lines.append("|       +-- molecules/             Per-molecule SDF files")
    lines.append("|       +-- configs/               Per-molecule docking configs")
    if has_smina:
        lines.append(
            "|       +-- smina/                 SMINA per-molecule results & logs"
        )
    if has_gnina:
        lines.append(
            "|       +-- gnina/                 GNINA per-molecule results & logs"
        )
    if has_matcha:
        lines.append("|       +-- run_matcha.sh          Matcha launch script")

    # Batch scripts and config files present in _workdir
    workdir = docking_dir / "_workdir"
    if workdir.is_dir():
        for name in sorted(workdir.iterdir()):
            if name.is_file() and name.suffix in (".sh", ".ini"):
                lines.append(f"|       +-- {name.name}")
    else:
        # Infer from tools
        if has_smina:
            lines.append("|       +-- run_smina.sh           SMINA run script")
        if has_gnina:
            lines.append("|       +-- gnina_config.ini        GNINA configuration")
            lines.append("|       +-- run_gnina.sh           GNINA run script")

    return lines


def _count_stage_molecules(base_path: Path, stage_dir: str) -> int | None:
    """Count molecules in a stage's filtered_molecules.csv.

    Returns the row count (excluding header), or None if the file is missing.
    """
    csv_path = base_path / stage_dir / "filtered" / FILE_FILTERED_MOLECULES
    if not csv_path.exists():
        csv_path = base_path / stage_dir / FILE_FILTERED_MOLECULES
    if not csv_path.exists():
        return None
    try:
        # Count lines minus header for speed (avoid pandas overhead)
        with open(csv_path) as f:
            return max(sum(1 for _ in f) - 1, 0)
    except Exception:
        return None


def _format_duration(seconds: float) -> str:
    """Format seconds into a human-readable duration string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m {int(secs)}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h {int(minutes)}m"


def _build_stage_tree_section(
    stage_name: str, stage_dir: str, base_path: Path, config: dict | None
) -> str:
    """Build a single stage's tree block for the README."""
    desc = _STAGE_DESCRIPTIONS.get(stage_name, "")
    dir_basename = stage_dir.split("/", 1)[1] if "/" in stage_dir else stage_dir
    header = f"|   +-- {dir_basename + '/':38s}{desc}"

    if stage_name == STAGE_DOCKING:
        body_lines = _build_docking_tree(base_path, config)
    else:
        body_lines = list(_STAGE_TREE_TEMPLATES.get(stage_name, []))

    lines = [header] + body_lines + ["|"]
    return "\n".join(lines) + "\n"


def _build_pipeline_flow_table(
    base_path: Path,
    stages: list[PipelineStage],
    initial_count: int,
    stage_timings: dict[str, float] | None,
) -> str:
    """Build a markdown table showing molecule counts and timings per stage."""
    stage_timings = stage_timings or {}

    rows: list[tuple[str, str, str]] = []
    rows.append(("Input", str(initial_count), ""))

    for stage in stages:
        if not stage.enabled:
            continue
        label = _STAGE_DESCRIPTIONS.get(stage.name, stage.name)
        count = _count_stage_molecules(base_path, stage.directory)
        count_str = str(count) if count is not None else "-"
        timing = stage_timings.get(stage.name)
        time_str = _format_duration(timing) if timing is not None else ""
        rows.append((label, count_str, time_str))

    # Calculate column widths
    col0_w = max(len(r[0]) for r in rows)
    col1_w = max(len(r[1]) for r in rows)
    col2_w = max((len(r[2]) for r in rows), default=0)

    header = f"| {'Stage':<{col0_w}} | {'Molecules':>{col1_w}} | {'Time':>{col2_w}} |"
    sep = f"|{'-' * (col0_w + 2)}|{'-' * (col1_w + 2)}|{'-' * (col2_w + 2)}|"
    body = "\n".join(
        f"| {r[0]:<{col0_w}} | {r[1]:>{col1_w}} | {r[2]:>{col2_w}} |" for r in rows
    )

    return f"{header}\n{sep}\n{body}"


# -- README template and data preparation (Task 4) --

_STRUCTURE_README_TEMPLATE = """\
# HEDGEHOG Run Info

Generated: {generated_at}

## Summary

- Initial molecules: {initial_count}
- Final molecules: {final_count}
- Retention rate: {retention}
- Stages enabled: {enabled_count}
- Stages completed: {completed_count}
{total_time_line}

## Pipeline Flow

{flow_table}

## Directory Structure

```
{dir_name}/
+-- input/
|   +-- sampled_molecules.csv           Input molecules
|
+-- stages/                             Pipeline stages (numbered by execution order)
{stages_content}+-- output/                             Final results
|   +-- final_molecules.csv            Final filtered molecules
|
+-- report.html                         Interactive HTML report
+-- report_data.json                    Report data
|
+-- configs/                            Configuration snapshots & runtime config
|   +-- master_config_resolved.yml
|   +-- config_*.yml
|   +-- run_models_mapping.csv
|
+-- run_YYYYMMDD_HHMMSS.log             Pipeline log file
```

## File Naming Conventions

- `filtered_molecules.csv` - Molecules that passed filters
- `failed_molecules.csv` - Molecules that failed filters
- `descriptors_all.csv` - All computed descriptors
- `metrics.csv` - Summary statistics
- `extended.csv` - Detailed results with all columns

## Stage Execution Summary

{stage_summary}
"""


def _prepare_readme_data(
    base_path: Path,
    stages: list[PipelineStage],
    initial_count: int,
    final_count: int,
    stage_timings: dict[str, float] | None = None,
    config: dict | None = None,
) -> dict:
    """Prepare template variables for _STRUCTURE_README_TEMPLATE."""
    stage_timings = stage_timings or {}
    enabled_stages = {s.name for s in stages if s.enabled}

    retention = (
        f"{100 * final_count / initial_count:.2f}%" if initial_count > 0 else "N/A"
    )

    # Build stage tree sections dynamically
    stage_sections = []
    for stage in stages:
        if stage.name in enabled_stages:
            stage_sections.append(
                _build_stage_tree_section(
                    stage.name, stage.directory, base_path, config
                )
            )
    stages_content = "".join(stage_sections)

    # Build pipeline flow table
    flow_table = _build_pipeline_flow_table(
        base_path, stages, initial_count, stage_timings
    )

    # Build stage execution summary with timings
    stage_summary_lines = []
    for stage in stages:
        if stage.completed:
            status = "COMPLETED"
        elif stage.enabled:
            status = "FAILED"
        else:
            status = "DISABLED"
        timing = stage_timings.get(stage.name)
        if timing is not None:
            stage_summary_lines.append(
                f"- **{stage.name}**: {status} ({_format_duration(timing)})"
            )
        else:
            stage_summary_lines.append(f"- **{stage.name}**: {status}")

    total_time = sum(stage_timings.values()) if stage_timings else None

    return {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "initial_count": initial_count,
        "final_count": final_count,
        "retention": retention,
        "enabled_count": len(enabled_stages),
        "completed_count": sum(1 for s in stages if s.completed),
        "total_time_line": (
            f"- Total pipeline time: {_format_duration(total_time)}"
            if total_time
            else ""
        ),
        "flow_table": flow_table,
        "dir_name": base_path.name,
        "stages_content": stages_content,
        "stage_summary": "\n".join(stage_summary_lines),
    }


def _generate_structure_readme(
    base_path: Path,
    stages: list[PipelineStage],
    initial_count: int,
    final_count: int,
    stage_timings: dict[str, float] | None = None,
    config: dict | None = None,
) -> None:
    """Generate RUN_INFO.md documenting the output structure for this run."""
    try:
        readme_path = base_path / "RUN_INFO.md"
        data = _prepare_readme_data(
            base_path, stages, initial_count, final_count, stage_timings, config
        )
        content = _STRUCTURE_README_TEMPLATE.format(**data)

        with open(readme_path, "w") as f:
            f.write(content)

        logger.info("Generated run info: %s", readme_path)
    except Exception as e:
        logger.warning("Failed to generate RUN_INFO.md: %s", e)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def calculate_metrics(data, config: dict, progress_callback=None) -> bool:
    """Calculate metrics for molecular data using the configured pipeline.

    Args:
        data: Input molecular data (pandas DataFrame)
        config: Pipeline configuration dictionary

    Returns:
        True if all enabled stages completed successfully, False otherwise
    """
    folder = Path(config[CONFIG_FOLDER_TO_SAVE])
    incomplete_marker = folder / FILE_RUN_INCOMPLETE
    success = False

    try:
        # Drop a marker so that interrupted runs are detectable.
        folder.mkdir(parents=True, exist_ok=True)
        incomplete_marker.write_text(
            f"Pipeline started: {datetime.now().isoformat()}\n"
            "This file exists while the pipeline is running or if it was interrupted/failed.\n"
            "It will be removed automatically on successful completion.\n"
            f"If the pipeline is still running please check progress in: {folder}/_workdir/\n"
            "\n"
            f"If the pipeline failed please review the run log: {folder}/run_*.log.\n"
        )

        _save_config_snapshot(config)
        pipeline = MolecularAnalysisPipeline(config, progress_callback)
        success = pipeline.run_pipeline(data)

        # Pipeline finished normally -- remove the marker.
        incomplete_marker.unlink(missing_ok=True)

        return success
    except InterruptedError:
        logger.info("Pipeline execution cancelled")
        raise
    except Exception:
        logger.exception("Pipeline execution failed")
        success = False
    finally:
        _cleanup_lingering_processes()

    return success


def _cleanup_lingering_processes(timeout_seconds: float = 2.0) -> None:
    """Best-effort cleanup of leftover workers without long shutdown stalls."""
    # Some third-party libraries (e.g., joblib/loky) can leave worker processes
    # alive past pipeline completion. Use a global timeout budget so cleanup time
    # does not scale linearly with worker count.
    try:
        children = multiprocessing.active_children()
    except Exception:
        children = []

    deadline = time.monotonic() + max(0.0, timeout_seconds)

    for child in children:
        try:
            if child.is_alive():
                child.terminate()
        except Exception:
            pass

    for child in children:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break
        try:
            child.join(timeout=min(0.2, remaining))
        except Exception:
            pass

    # If processes are still alive after the grace period, hard-kill them.
    for child in children:
        try:
            if child.is_alive():
                child.kill()
        except Exception:
            pass

    for child in children:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break
        try:
            if child.is_alive():
                child.join(timeout=min(0.1, remaining))
        except Exception:
            pass

    try:
        from joblib.externals.loky import reusable_executor as loky_reusable_executor

        loky_executor = getattr(loky_reusable_executor, "_executor", None)
        if loky_executor is not None:
            loky_executor.shutdown(wait=True, kill_workers=True)
    except Exception:
        pass

    # A second pass helps when child processes were created by external runtimes
    # (e.g. loky) and were not visible during the first active_children() call.
    try:
        tail_children = multiprocessing.active_children()
    except Exception:
        tail_children = []

    for child in tail_children:
        try:
            if child.is_alive():
                child.terminate()
        except Exception:
            pass

    for child in tail_children:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break
        try:
            if child.is_alive():
                child.join(timeout=min(0.1, remaining))
        except Exception:
            pass

    # Prevent Python's atexit Pool finalizers from re-entering pool shutdown and
    # blocking interpreter exit when workers were already force-terminated above.
    try:
        from multiprocessing import util as mp_util

        finalizers = getattr(mp_util, "_finalizer_registry", {})
        for finalizer in list(finalizers.values()):
            callback = getattr(finalizer, "_callback", None)
            cb_name = getattr(
                callback,
                "__qualname__",
                getattr(callback, "__name__", ""),
            )
            if cb_name.endswith("Pool._terminate_pool"):
                finalizer.cancel()
    except Exception:
        pass
