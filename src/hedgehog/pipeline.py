import shutil
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from hedgehog.configs.logger import load_config, logger
from hedgehog.reporting import ReportGenerator
from hedgehog.stages.descriptors.main import main as descriptors_main
from hedgehog.stages.docking.utils import run_docking as docking_main
from hedgehog.stages.dockingFilters.main import docking_filters_main
from hedgehog.stages.structFilters.main import main as structural_filters_main
from hedgehog.stages.synthesis.main import main as synthesis_main
from hedgehog.utils.input_paths import find_latest_input_source as _find_input

# Directory names
DIR_INPUT = "input"
DIR_STAGES = "stages"
DIR_OUTPUT = "output"
DIR_CONFIGS = "configs"

# Stage subdirectories
DIR_DESCRIPTORS_INITIAL = "stages/01_descriptors_initial"
DIR_STRUCT_FILTERS_PRE = "stages/02_structural_filters_pre"
DIR_STRUCT_FILTERS_POST = "stages/03_structural_filters_post"
DIR_SYNTHESIS = "stages/04_synthesis"
DIR_DOCKING = "stages/05_docking"
DIR_DOCKING_FILTERS = "stages/06_docking_filters"
DIR_DESCRIPTORS_FINAL = "stages/07_descriptors_final"

# Legacy names for backwards compatibility
DIR_DESCRIPTORS = DIR_DESCRIPTORS_INITIAL
DIR_STRUCT_FILTERS = DIR_STRUCT_FILTERS_POST
DIR_BEFORE_DESCRIPTORS = DIR_STRUCT_FILTERS_PRE
DIR_FINAL_DESCRIPTORS = DIR_DESCRIPTORS_FINAL
DIR_RUN_CONFIGS = DIR_CONFIGS

# File names
FILE_SAMPLED_MOLECULES = "sampled_molecules.csv"
FILE_FINAL_MOLECULES = "final_molecules.csv"
FILE_FILTERED_MOLECULES = "filtered_molecules.csv"
FILE_PASS_SMILES_TEMPLATE = "filtered_molecules.csv"
FILE_MASTER_CONFIG = "master_config_resolved.yml"
FILE_GNINA_OUTPUT = "gnina_out.sdf"

# Stage names
STAGE_STRUCT_INI_FILTERS = "struct_ini_filters"
STAGE_DESCRIPTORS = "descriptors"
STAGE_STRUCT_FILTERS = "struct_filters"
STAGE_SYNTHESIS = "synthesis"
STAGE_DOCKING = "docking"
STAGE_DOCKING_FILTERS = "docking_filters"
STAGE_FINAL_DESCRIPTORS = "final_descriptors"

# Config keys
CONFIG_DESCRIPTORS = "config_descriptors"
CONFIG_STRUCT_FILTERS = "config_structFilters"
CONFIG_SYNTHESIS = "config_synthesis"
CONFIG_DOCKING = "config_docking"
CONFIG_DOCKING_FILTERS = "config_docking_filters"
CONFIG_RUN_KEY = "run"
CONFIG_RUN_BEFORE_DESCRIPTORS = "run_before_descriptors"
CONFIG_TOOLS = "tools"
CONFIG_FOLDER_TO_SAVE = "folder_to_save"

# Command-line override keys
OVERRIDE_STRUCT_FILTERS_BEFORE = "_run_struct_filters_before_descriptors_override"
OVERRIDE_SINGLE_STAGE = "_run_single_stage_override"

# Docking tools
DOCKING_TOOL_SMINA = "smina"
DOCKING_TOOL_GNINA = "gnina"
DOCKING_TOOL_BOTH = "both"
DOCKING_RESULTS_DIR_TEMPLATE = {
    DOCKING_TOOL_SMINA: f"{DIR_DOCKING}/smina",
    DOCKING_TOOL_GNINA: f"{DIR_DOCKING}/gnina",
}


FILE_RUN_INCOMPLETE = ".RUN_INCOMPLETE"


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

    def start(self, message: str | None = None) -> None:
        self._emit_event(
            {
                "type": "stage_start",
                "stage": self.stage,
                "stage_index": self.stage_index,
                "total_stages": self.total_stages,
                "message": message,
            }
        )

    def progress(self, current: int, total: int, message: str | None = None) -> None:
        self._emit_event(
            {
                "type": "stage_progress",
                "stage": self.stage,
                "stage_index": self.stage_index,
                "total_stages": self.total_stages,
                "current": int(current),
                "total": int(total),
                "message": message,
            }
        )

    def complete(self, ok: bool = True, message: str | None = None) -> None:
        self._emit_event(
            {
                "type": "stage_complete",
                "stage": self.stage,
                "stage_index": self.stage_index,
                "total_stages": self.total_stages,
                "ok": bool(ok),
                "message": message,
            }
        )


def _log_stage_header(stage_label: str) -> None:
    """Log a formatted stage header."""
    separator = "[#B29EEE]" + "\u2500" * 59 + "[/#B29EEE]"
    logger.info("")
    logger.info(separator)
    logger.info("[#B29EEE]  %s[/#B29EEE]", stage_label)
    logger.info(separator)
    logger.info("")


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
        DIR_DESCRIPTORS_INITIAL: Path(DIR_DESCRIPTORS_INITIAL)
        / "filtered"
        / FILE_FILTERED_MOLECULES,
        DIR_STRUCT_FILTERS_POST: Path(DIR_STRUCT_FILTERS_POST)
        / FILE_FILTERED_MOLECULES,
        DIR_STRUCT_FILTERS_PRE: Path(DIR_STRUCT_FILTERS_PRE) / FILE_FILTERED_MOLECULES,
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
        DIR_STRUCT_FILTERS_PRE,
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
        except Exception as e:
            logger.error("Error running descriptors: %s", e)
            return False

    def run_structural_filters(
        self, stage_dir: str, reporter: StageProgressReporter | None = None
    ) -> bool:
        """Run structural filters on molecules."""
        try:
            is_post_descriptors = stage_dir != DIR_STRUCT_FILTERS_PRE
            has_descriptor_data = self.data_checker.check_stage_data(
                DIR_DESCRIPTORS_INITIAL
            )

            if is_post_descriptors and not has_descriptor_data:
                if self.config.get(OVERRIDE_SINGLE_STAGE) == STAGE_STRUCT_FILTERS:
                    logger.info(
                        "No previous stage data found, will use molecules from config"
                    )
                else:
                    logger.warning(
                        "No data available for structural filters in %s", stage_dir
                    )
                    return False

            config_struct_filters = load_config(self.config[CONFIG_STRUCT_FILTERS])
            if not config_struct_filters.get(CONFIG_RUN_KEY, False):
                logger.info("Structural filters disabled in config")
                return False

            structural_filters_main(self.config, stage_dir, reporter=reporter)
            return True
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

        if DOCKING_TOOL_BOTH in tools_list or not tools_list:
            return [DOCKING_TOOL_SMINA, DOCKING_TOOL_GNINA]
        return tools_list

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
                smina_dir = (
                    base_folder / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_SMINA]
                )
                if _directory_has_files(smina_dir):
                    return True
                logger.debug("SMINA results not found in %s", smina_dir)
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
                logger.error(
                    "Docking finished but no results detected in output directories"
                )
                return False
            return True
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
        except Exception as e:
            logger.error("Error running docking filters: %s", e)
            return False


class MolecularAnalysisPipeline:
    """Orchestrates the execution of the molecular analysis pipeline."""

    # Stage definitions: (name, config_key, directory)
    _STAGE_DEFINITIONS = [
        (STAGE_STRUCT_INI_FILTERS, CONFIG_STRUCT_FILTERS, DIR_BEFORE_DESCRIPTORS),
        (STAGE_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_DESCRIPTORS),
        (STAGE_STRUCT_FILTERS, CONFIG_STRUCT_FILTERS, DIR_STRUCT_FILTERS),
        (STAGE_SYNTHESIS, CONFIG_SYNTHESIS, DIR_SYNTHESIS),
        (STAGE_DOCKING, CONFIG_DOCKING, DIR_DOCKING),
        (STAGE_DOCKING_FILTERS, CONFIG_DOCKING_FILTERS, DIR_DOCKING_FILTERS),
        (STAGE_FINAL_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_FINAL_DESCRIPTORS),
    ]

    # Stage display labels for logging
    _STAGE_LABELS = {
        STAGE_STRUCT_INI_FILTERS: "Stage 1': Pre-descriptors Structural Filters",
        STAGE_DESCRIPTORS: "Stage 1: Molecular Descriptors",
        STAGE_STRUCT_FILTERS: "Stage 2: Post-descriptors Structural Filters",
        STAGE_SYNTHESIS: "Stage 3: Synthesis Analysis",
        STAGE_DOCKING: "Stage 4: Molecular Docking",
        STAGE_DOCKING_FILTERS: "Stage 5: Docking Filters",
        STAGE_FINAL_DESCRIPTORS: "Stage 6': Final Descriptors Calculation",
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
        self._initialize_stages()

    def _initialize_stages(self) -> None:
        """Initialize stages by loading their configs and determining enabled status."""
        single_stage_override = self.config.get(OVERRIDE_SINGLE_STAGE)

        for stage in self.stages:
            try:
                stage_config = load_config(self.config[stage.config_key])
                stage.enabled = stage_config.get(CONFIG_RUN_KEY, False)

                if stage.name == STAGE_STRUCT_INI_FILTERS:
                    run_before = self.config.get(
                        OVERRIDE_STRUCT_FILTERS_BEFORE,
                        stage_config.get(CONFIG_RUN_BEFORE_DESCRIPTORS, False),
                    )
                    stage.enabled = stage.enabled and run_before

                if single_stage_override:
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

    def run_pipeline(self, data) -> bool:
        """Execute the full molecular analysis pipeline."""
        self.current_data = data
        success_count = 0
        total_enabled = sum(1 for s in self.stages if s.enabled)

        # Execute each stage
        stage_results = [
            self._run_pre_descriptors_filters(),
            self._run_descriptors(data),
            self._run_post_descriptors_filters(),
            self._run_synthesis(),
            self._run_docking(),
            self._run_docking_filters(),
            self._run_final_descriptors(),
        ]

        # Check for early exit conditions
        stage_names = [defn[0] for defn in self._STAGE_DEFINITIONS]
        for name, (completed, early_exit) in zip(stage_names, stage_results):
            if completed:
                self._stage_by_name[name].completed = True
                success_count += 1
            if early_exit:
                return self._finalize_pipeline(data, success_count, total_enabled)

        logger.info(
            "Pipeline completed: %d/%d stages successful", success_count, total_enabled
        )
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
        current_idx = enabled_stages.index(stage.name) if stage.name in enabled_stages else 0
        reporter = StageProgressReporter(
            emit_event=self._emit_progress_event,
            stage=stage.name,
            stage_index=current_idx + 1,
            total_stages=len(enabled_stages),
        )
        reporter.start()

        _log_stage_header(self._STAGE_LABELS[stage.name])
        start_time = time.perf_counter()
        try:
            completed = runner_func(*args, reporter=reporter)
        except TypeError:
            completed = runner_func(*args)
        elapsed = time.perf_counter() - start_time
        self.stage_timings[stage.name] = elapsed
        logger.info("Stage %s completed in %.1f seconds", stage.name, elapsed)
        reporter.complete(ok=bool(completed))

        if completed:
            return True, False
        if on_failure:
            return on_failure()
        return False, False

    def _run_pre_descriptors_filters(self) -> tuple[bool, bool]:
        """Run pre-descriptors structural filters stage."""
        return self._run_stage(
            STAGE_STRUCT_INI_FILTERS,
            self.stage_runner.run_structural_filters,
            DIR_STRUCT_FILTERS_PRE,
        )

    def _run_descriptors(self, data) -> tuple[bool, bool]:
        """Run descriptors calculation stage."""
        descriptors_input = data

        # If pre-descriptors structural filters are enabled, use their output as the
        # input to descriptors (otherwise stage 1' becomes analysis-only and the
        # molecule counts become inconsistent across stages).
        pre_stage = self._stage_by_name[STAGE_STRUCT_INI_FILTERS]
        if pre_stage.enabled:
            pre_path = self._build_data_path(DIR_STRUCT_FILTERS_PRE)
            if pre_path.exists():
                try:
                    pre_df = pd.read_csv(pre_path)
                    if len(pre_df) > 0:
                        descriptors_input = pre_df
                        self.current_data = pre_df
                except Exception as e:
                    logger.warning(
                        "Could not load pre-descriptors structural filters output (%s): %s",
                        pre_path,
                        e,
                    )

        return self._run_stage(
            STAGE_DESCRIPTORS, self.stage_runner.run_descriptors, descriptors_input
        )

    def _run_post_descriptors_filters(self) -> tuple[bool, bool]:
        """Run post-descriptors structural filters stage."""

        def _on_failure():
            logger.info("No molecules left after descriptors; ending pipeline early.")
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
        return self._run_stage(STAGE_DOCKING, self.stage_runner.run_docking)

    def _run_docking_filters(self) -> tuple[bool, bool]:
        """Run docking filters stage."""
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

    def _finalize_pipeline(self, data, success_count: int, total_enabled: int) -> bool:
        """Finalize pipeline execution with summary and output."""
        self._log_pipeline_summary()

        initial_count = len(data)
        final_data = self.get_latest_data(
            skip_descriptors=True, fallback_on_empty=False
        )
        final_count = len(final_data) if final_data is not None else 0

        self._log_molecule_summary(initial_count, final_count)
        self._save_final_output(final_data, final_count)
        _generate_structure_readme(
            self.data_checker.base_path,
            self.stages,
            initial_count,
            final_count,
            stage_timings=self.stage_timings,
            config=self.config,
        )
        self._generate_html_report(initial_count, final_count)

        return not any(self._stage_is_failed(stage) for stage in self.stages)

    def _stage_is_failed(self, stage: PipelineStage) -> bool:
        """Return True if an enabled stage is considered a failure.

        Stages are not treated as failures when they are skipped due to having
        no molecules available at their required input boundary.
        """
        if not stage.enabled:
            return False
        if stage.completed:
            return False

        if stage.name == STAGE_FINAL_DESCRIPTORS:
            # Final descriptors depends on the latest non-descriptor stage output.
            # If that latest output is empty (0 molecules), we treat it as skipped,
            # not as a pipeline failure.
            sources = [
                DIR_DOCKING_FILTERS,
                DIR_SYNTHESIS,
                DIR_STRUCT_FILTERS_POST,
                DIR_STRUCT_FILTERS_PRE,
            ]
            latest = self._find_data_source(sources)
            if not latest:
                return False
            try:
                latest_df = pd.read_csv(self._build_data_path(latest))
            except Exception:
                # If the latest output exists but cannot be read, treat as failure.
                return True
            return len(latest_df) > 0

        # Other stages are skippable when their required upstream output is absent/empty.
        skip_conditions = {
            STAGE_STRUCT_FILTERS: DIR_DESCRIPTORS,
            STAGE_SYNTHESIS: DIR_STRUCT_FILTERS,
        }
        required_data = skip_conditions.get(stage.name)
        if required_data and not self.data_checker.stage_has_molecules(required_data):
            return False

        return True

    def _generate_html_report(self, initial_count: int, final_count: int) -> None:
        """Generate HTML report for the pipeline run."""
        try:
            report_generator = ReportGenerator(
                base_path=self.data_checker.base_path,
                stages=self.stages,
                config=self.config,
                initial_count=initial_count,
                final_count=final_count,
            )
            report_path = report_generator.generate()
            logger.info("HTML report generated: %s", report_path)
        except Exception as e:
            logger.warning("Failed to generate HTML report: %s", e)

    def _log_molecule_summary(self, initial_count: int, final_count: int) -> None:
        """Log molecule count summary."""
        logger.info("---------> Molecule Count Summary")
        logger.info("Molecules before filtration: %d", initial_count)
        logger.info("Molecules remaining after all stages: %d", final_count)
        if initial_count > 0:
            retention = 100 * final_count / initial_count
            logger.info("Retention rate: %.2f%%", retention)

    def _save_final_output(self, final_data, final_count: int) -> None:
        """Save final molecules to output directory."""
        final_output_path = (
            self.data_checker.base_path / DIR_OUTPUT / FILE_FINAL_MOLECULES
        )
        final_output_path.parent.mkdir(parents=True, exist_ok=True)

        id_cols = ["smiles", "model_name", "mol_idx"]

        if final_data is None or len(final_data) == 0:
            # Always create the file so downstream tooling can rely on its presence.
            pd.DataFrame(columns=id_cols).to_csv(final_output_path, index=False)
            logger.info("Saved 0 final molecules to %s", final_output_path)
            return

        cols = [c for c in id_cols if c in final_data.columns]
        final_data[cols].to_csv(final_output_path, index=False)
        logger.info("Saved %d final molecules to %s", final_count, final_output_path)

    def _log_pipeline_summary(self) -> None:
        """Log a summary of pipeline execution status and timings."""
        _log_stage_header("Pipeline Execution Summary")

        for stage in self.stages:
            status = self._get_stage_status(stage)
            timing = self.stage_timings.get(stage.name)
            if timing is not None:
                logger.info("%s: %s (%.1fs)", stage.name, status, timing)
            else:
                logger.info("%s: %s", stage.name, status)

        # Log total time for completed stages
        if self.stage_timings:
            total_time = sum(self.stage_timings.values())
            logger.info("Total pipeline time: %.1f seconds", total_time)

    def _get_stage_status(self, stage: PipelineStage) -> str:
        """Get display status for a stage."""
        if not stage.enabled:
            return "DISABLED"
        if stage.completed:
            return "[#B29EEE]\u2713 COMPLETED[/#B29EEE]"

        if stage.name == STAGE_FINAL_DESCRIPTORS:
            # Final descriptors can run on the latest available non-descriptor output
            # (docking filters, synthesis, structural filters, etc.). Mark as skipped
            # only if that latest output has no molecules.
            sources = [
                DIR_DOCKING_FILTERS,
                DIR_SYNTHESIS,
                DIR_STRUCT_FILTERS_POST,
                DIR_STRUCT_FILTERS_PRE,
            ]
            latest = self._find_data_source(sources)
            if not latest:
                return "[yellow]\u27c2 SKIPPED (no molecules)[/yellow]"
            try:
                latest_df = pd.read_csv(self._build_data_path(latest))
            except Exception:
                return "[red]\u2717 FAILED[/red]"
            if len(latest_df) == 0:
                return "[yellow]\u27c2 SKIPPED (no molecules)[/yellow]"

        # Check for skipped conditions
        skip_conditions = {
            STAGE_STRUCT_FILTERS: DIR_DESCRIPTORS,
            STAGE_SYNTHESIS: DIR_STRUCT_FILTERS,
        }

        required_data = skip_conditions.get(stage.name)
        if required_data and not self.data_checker.stage_has_molecules(required_data):
            return "[yellow]\u27c2 SKIPPED (no molecules)[/yellow]"

        return "[red]\u2717 FAILED[/red]"


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


# Stage descriptions for README generation
_STAGE_DESCRIPTIONS = {
    STAGE_STRUCT_INI_FILTERS: "Pre-descriptors structural filters",
    STAGE_DESCRIPTORS: "Physicochemical descriptors",
    STAGE_STRUCT_FILTERS: "Post-descriptors structural filters",
    STAGE_SYNTHESIS: "Retrosynthesis analysis",
    STAGE_DOCKING: "Molecular docking",
    STAGE_DOCKING_FILTERS: "Post-docking pose & interaction filters",
    STAGE_FINAL_DESCRIPTORS: "Final descriptor calculation",
}

# Directory-level tree templates per stage (static parts)
_STAGE_TREE_TEMPLATES: dict[str, list[str]] = {
    STAGE_STRUCT_INI_FILTERS: [
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

    # Fall back to config if directories don't exist yet
    if not has_smina and not has_gnina and config:
        docking_cfg_path = config.get(CONFIG_DOCKING)
        if docking_cfg_path:
            try:
                docking_cfg = load_config(docking_cfg_path)
                tools = docking_cfg.get(CONFIG_TOOLS, "")
                has_smina = tools in (DOCKING_TOOL_SMINA, DOCKING_TOOL_BOTH)
                has_gnina = tools in (DOCKING_TOOL_GNINA, DOCKING_TOOL_BOTH)
            except Exception:
                has_smina = has_gnina = True  # Show both as fallback

    lines = [
        "|   +-- ligands.csv                Prepared ligands",
        "|   +-- job_meta.json              Job metadata",
    ]
    if has_smina:
        lines.append("|   +-- smina/")
        lines.append("|   |   +-- smina_out.sdf          Aggregated SMINA results (1 pose/molecule)")
    if has_gnina:
        lines.append("|   +-- gnina/")
        lines.append("|   |   +-- gnina_out.sdf          Aggregated GNINA results (1 pose/molecule)")

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
        enabled_stages = {s.name for s in stages if s.enabled}
        enabled_count = len(enabled_stages)
        completed_count = sum(1 for s in stages if s.completed)

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
        stage_timings = stage_timings or {}
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
        stage_summary = "\n".join(stage_summary_lines)

        total_time = sum(stage_timings.values()) if stage_timings else None

        content = f"""\
# HEDGEHOG Run Info

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Summary

- Initial molecules: {initial_count}
- Final molecules: {final_count}
- Retention rate: {retention}
- Stages enabled: {enabled_count}
- Stages completed: {completed_count}
{f"- Total pipeline time: {_format_duration(total_time)}" if total_time else ""}

## Pipeline Flow

{flow_table}

## Directory Structure

```
{base_path.name}/
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

        with open(readme_path, "w") as f:
            f.write(content)

        logger.info("Generated run info: %s", readme_path)
    except Exception as e:
        logger.warning("Failed to generate RUN_INFO.md: %s", e)


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

    try:
        # Drop a marker so that interrupted runs are detectable.
        folder.mkdir(parents=True, exist_ok=True)
        incomplete_marker.write_text(
            f"Pipeline started: {datetime.now().isoformat()}\n"
            "This file is removed automatically on successful completion.\n"
            "If you see it, the run was interrupted or failed.\n"
        )

        _save_config_snapshot(config)
        pipeline = MolecularAnalysisPipeline(config, progress_callback)
        success = pipeline.run_pipeline(data)

        # Pipeline finished normally — remove the marker.
        incomplete_marker.unlink(missing_ok=True)
        return success
    except Exception as e:
        logger.error("Pipeline execution failed: %s", e)
        return False
