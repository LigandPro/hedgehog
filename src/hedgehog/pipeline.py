import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.descriptors.main import main as descriptors_main
from hedgehog.stages.docking.utils import run_docking as docking_main
from hedgehog.stages.structFilters.main import main as structural_filters_main
from hedgehog.stages.synthesis.main import main as synthesis_main

# Directory names
DIR_INPUT = "input"
DIR_STAGES = "stages"
DIR_OUTPUT = "output"
DIR_CONFIGS = "configs"
DIR_LOGS = "logs"

# Stage subdirectories
DIR_DESCRIPTORS_INITIAL = "stages/01_descriptors_initial"
DIR_STRUCT_FILTERS_PRE = "stages/02_structural_filters_pre"
DIR_STRUCT_FILTERS_POST = "stages/03_structural_filters_post"
DIR_SYNTHESIS = "stages/04_synthesis"
DIR_DOCKING = "stages/05_docking"
DIR_DESCRIPTORS_FINAL = "stages/06_descriptors_final"

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
STAGE_FINAL_DESCRIPTORS = "final_descriptors"

# Config keys
CONFIG_DESCRIPTORS = "config_descriptors"
CONFIG_STRUCT_FILTERS = "config_structFilters"
CONFIG_SYNTHESIS = "config_synthesis"
CONFIG_DOCKING = "config_docking"
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


def _log_stage_header(stage_label: str) -> None:
    """Log a formatted stage header."""
    logger.info("")
    logger.info("[#B29EEE]" + "\u2500" * 59 + "[/#B29EEE]")
    logger.info(f"[#B29EEE]  {stage_label}[/#B29EEE]")
    logger.info("[#B29EEE]" + "\u2500" * 59 + "[/#B29EEE]")
    logger.info("")


def _file_exists_and_not_empty(file_path: Path) -> bool:
    """Check if a file exists and is not empty."""
    try:
        return file_path.exists() and file_path.stat().st_size > 0
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
        # Legacy paths
        "Descriptors": Path("Descriptors") / "passDescriptorsSMILES.csv",
        "StructFilters": Path("StructFilters") / "passStructFiltersSMILES.csv",
        "beforeDescriptors": Path("beforeDescriptors_StructFilters")
        / "passStructFiltersSMILES.csv",
        "Synthesis": Path("Synthesis") / "passSynthesisSMILES.csv",
    }

    def __init__(self, config: dict):
        self.config = config
        self.base_path = Path(config[CONFIG_FOLDER_TO_SAVE])

    def check_stage_data(self, stage_name: str) -> bool:
        """Check if data exists for a given stage."""
        path = self._get_stage_output_path(stage_name.strip())
        return _file_exists_and_not_empty(path) if path else False

    def _get_stage_output_path(self, stage_name: str) -> Path | None:
        """Get the expected output file path for a stage."""
        relative_path = self._STAGE_OUTPUT_PATHS.get(stage_name)
        return self.base_path / relative_path if relative_path else None


class PipelineStageRunner:
    """Executes individual pipeline stages and manages stage data flow."""

    DATA_SOURCE_PRIORITY = [
        DIR_SYNTHESIS,
        DIR_STRUCT_FILTERS_POST,
        DIR_DESCRIPTORS_INITIAL,
        DIR_STRUCT_FILTERS_PRE,
    ]

    def __init__(self, config: dict, data_checker: DataChecker):
        self.config = config
        self.data_checker = data_checker

    def find_latest_data_source(self) -> str | None:
        """Find the most recent stage with available output data."""
        for source in self.DATA_SOURCE_PRIORITY:
            if self.data_checker.check_stage_data(source):
                logger.debug(f"Found data from stage: {source}")
                return source
        logger.debug("No processed data found from any stage")
        return None

    def run_descriptors(self, data, subfolder: str | None = None) -> bool:
        """Run molecular descriptors calculation."""
        try:
            config_descriptors = load_config(self.config[CONFIG_DESCRIPTORS])
            if not config_descriptors.get(CONFIG_RUN_KEY, False):
                logger.info("Descriptors calculation disabled in config")
                return False
            descriptors_main(data, self.config, subfolder=subfolder)
            return True
        except Exception as e:
            logger.error(f"Error running descriptors: {e}")
            return False

    def run_structural_filters(self, stage_dir: str) -> bool:
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
                        f"No data available for structural filters in {stage_dir}"
                    )
                    return False

            config_struct_filters = load_config(self.config[CONFIG_STRUCT_FILTERS])
            if not config_struct_filters.get(CONFIG_RUN_KEY, False):
                logger.info("Structural filters disabled in config")
                return False

            structural_filters_main(self.config, stage_dir)
            return True
        except Exception as e:
            logger.error(f"Error running structural filters: {e}")
            return False

    def run_synthesis(self) -> bool:
        """Run synthesis analysis."""
        try:
            if not self._validate_synthesis_input():
                return False

            config_synthesis = load_config(self.config[CONFIG_SYNTHESIS])
            if not config_synthesis.get(CONFIG_RUN_KEY, False):
                logger.info("Synthesis disabled in config")
                return False

            synthesis_main(self.config)

            output_path = (
                self.data_checker.base_path / DIR_SYNTHESIS / FILE_FILTERED_MOLECULES
            )
            if not output_path.exists():
                logger.error("Synthesis finished but no output file detected")
                return False
            return True
        except Exception as e:
            logger.error(f"Error running synthesis: {e}")
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
            logger.info(f"Using {FILE_SAMPLED_MOLECULES}")
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
                logger.debug(f"GNINA results not found at {gnina_sdf}")
            elif tool == DOCKING_TOOL_SMINA:
                smina_dir = (
                    base_folder / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_SMINA]
                )
                if _directory_has_files(smina_dir):
                    return True
                logger.debug(f"SMINA results not found in {smina_dir}")
        return False

    def _get_gnina_output_dir(self, cfg: dict, base_folder: Path) -> Path:
        """Get the GNINA output directory from config."""
        gnina_config = cfg.get("gnina_config", {})
        cfg_out_dir = gnina_config.get("output_dir") or cfg.get("gnina_output_dir")
        if cfg_out_dir:
            out_dir = Path(cfg_out_dir)
            return out_dir if out_dir.is_absolute() else base_folder / out_dir
        return base_folder / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_GNINA]

    def run_docking(self) -> bool:
        """Run molecular docking."""
        try:
            config_docking = load_config(self.config[CONFIG_DOCKING])
            if not config_docking.get(CONFIG_RUN_KEY, False):
                logger.info("Docking disabled in config")
                return False

            if not docking_main(self.config):
                return False

            if not self.docking_results_present():
                logger.error(
                    "Docking finished but no results detected in output directories"
                )
                return False
            return True
        except Exception as e:
            logger.error(f"Error running docking: {e}")
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
        (STAGE_FINAL_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_FINAL_DESCRIPTORS),
    ]

    # Stage display labels for logging
    _STAGE_LABELS = {
        STAGE_STRUCT_INI_FILTERS: "Stage 1': Pre-descriptors Structural Filters",
        STAGE_DESCRIPTORS: "Stage 1: Molecular Descriptors",
        STAGE_STRUCT_FILTERS: "Stage 2: Post-descriptors Structural Filters",
        STAGE_SYNTHESIS: "Stage 3: Synthesis Analysis",
        STAGE_DOCKING: "Stage 4: Molecular Docking",
        STAGE_FINAL_DESCRIPTORS: "Stage 5': Final Descriptors Calculation",
    }

    def __init__(self, config: dict):
        self.config = config
        self.data_checker = DataChecker(config)
        self.stage_runner = PipelineStageRunner(config, self.data_checker)
        self.current_data = None
        self.stages = [PipelineStage(*defn) for defn in self._STAGE_DEFINITIONS]
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
                        logger.info(f"Single stage mode: enabling only {stage.name}")

                logger.debug(
                    f"Stage {stage.name}: {'Enabled' if stage.enabled else 'Disabled'}"
                )
            except Exception as e:
                logger.warning(f"Could not load config for {stage.name}: {e}")
                stage.enabled = False

    def get_latest_data(self, skip_descriptors: bool = False):
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

            if len(data) == 0:
                return self._try_fallback_sources(priority, latest_source, data)
            return data
        except Exception as e:
            logger.warning(f"Could not load data from {latest_source}: {e}")
            return self.current_data

    def _find_data_source(self, priority: list[str]) -> str | None:
        """Find first available data source from priority list."""
        for source in priority:
            if self.data_checker.check_stage_data(source):
                logger.debug(f"Found data from stage: {source}")
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
        logger.warning(f"No molecules in {current_source}, trying previous step...")
        current_idx = priority.index(current_source)

        for next_source in priority[current_idx + 1 :]:
            if self.data_checker.check_stage_data(next_source):
                path = self._build_data_path(next_source)
                try:
                    data = pd.read_csv(path)
                    if len(data) > 0:
                        logger.info(
                            f"Loaded data from {next_source}: {len(data)} molecules (previous step)"
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
            self._run_final_descriptors(),
        ]

        # Check for early exit conditions
        for i, (completed, early_exit) in enumerate(stage_results):
            if completed:
                self.stages[i].completed = True
                success_count += 1
            if early_exit:
                return self._finalize_pipeline(data, success_count, total_enabled)

        logger.info(
            f"Pipeline completed: {success_count}/{total_enabled} stages successful"
        )
        return self._finalize_pipeline(data, success_count, total_enabled)

    def _run_stage(self, stage_idx: int, runner_func, *args) -> tuple[bool, bool]:
        """Run a pipeline stage with standard logging. Returns (completed, early_exit)."""
        stage = self.stages[stage_idx]
        if not stage.enabled:
            return False, False

        _log_stage_header(self._STAGE_LABELS[stage.name])
        completed = runner_func(*args)
        return completed, False

    def _run_pre_descriptors_filters(self) -> tuple[bool, bool]:
        """Run pre-descriptors structural filters stage."""
        return self._run_stage(
            0, self.stage_runner.run_structural_filters, DIR_STRUCT_FILTERS_PRE
        )

    def _run_descriptors(self, data) -> tuple[bool, bool]:
        """Run descriptors calculation stage."""
        return self._run_stage(1, self.stage_runner.run_descriptors, data)

    def _run_post_descriptors_filters(self) -> tuple[bool, bool]:
        """Run post-descriptors structural filters stage."""
        stage = self.stages[2]
        if not stage.enabled:
            return False, False

        _log_stage_header(self._STAGE_LABELS[stage.name])
        if self.stage_runner.run_structural_filters(DIR_STRUCT_FILTERS_POST):
            return True, False

        logger.info("No molecules left after descriptors; ending pipeline early.")
        return False, True

    def _run_synthesis(self) -> tuple[bool, bool]:
        """Run synthesis analysis stage."""
        stage = self.stages[3]
        if not stage.enabled:
            return False, False

        _log_stage_header(self._STAGE_LABELS[stage.name])
        if self.stage_runner.run_synthesis():
            return True, False

        return self._handle_synthesis_failure()

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
        return self._run_stage(4, self.stage_runner.run_docking)

    def _run_final_descriptors(self) -> tuple[bool, bool]:
        """Run final descriptors calculation stage."""
        stage = self.stages[5]
        if not stage.enabled:
            return False, False

        _log_stage_header(self._STAGE_LABELS[stage.name])
        final_data = self.get_latest_data(skip_descriptors=True)

        if final_data is None or len(final_data) == 0:
            logger.info(
                "No molecules available for final descriptors calculation (skipping descriptors stage as source)"
            )
            if final_data is not None and len(final_data) == 0:
                logger.info(
                    "Ending pipeline: no molecules from previous steps (excluding descriptors)"
                )
            return False, False

        completed = self.stage_runner.run_descriptors(
            final_data, subfolder=DIR_FINAL_DESCRIPTORS
        )
        return completed, False

    def _finalize_pipeline(self, data, success_count: int, total_enabled: int) -> bool:
        """Finalize pipeline execution with summary and output."""
        self._log_pipeline_summary()

        initial_count = len(data)
        final_data = self.get_latest_data(skip_descriptors=True)
        final_count = len(final_data) if final_data is not None else 0

        self._log_molecule_summary(initial_count, final_count)
        self._save_final_output(final_data, final_count)
        _generate_structure_readme(
            self.data_checker.base_path, self.stages, initial_count, final_count
        )

        return success_count == total_enabled

    def _log_molecule_summary(self, initial_count: int, final_count: int) -> None:
        """Log molecule count summary."""
        logger.info("---------> Molecule Count Summary")
        logger.info(f"Molecules before filtration: {initial_count}")
        logger.info(f"Molecules remaining after all stages: {final_count}")
        if initial_count > 0:
            retention = 100 * final_count / initial_count
            logger.info(f"Retention rate: {retention:.2f}%")

    def _save_final_output(self, final_data, final_count: int) -> None:
        """Save final molecules to output directory."""
        if final_data is None or len(final_data) == 0:
            return

        final_output_path = (
            self.data_checker.base_path / DIR_OUTPUT / FILE_FINAL_MOLECULES
        )
        final_output_path.parent.mkdir(parents=True, exist_ok=True)
        final_data.to_csv(final_output_path, index=False)
        logger.info(f"Saved {final_count} final molecules to {final_output_path}")

    def _log_pipeline_summary(self) -> None:
        """Log a summary of pipeline execution status."""
        _log_stage_header("Pipeline Execution Summary")

        for stage in self.stages:
            status = self._get_stage_status(stage)
            logger.info(f"{stage.name}: {status}")

    def _get_stage_status(self, stage: PipelineStage) -> str:
        """Get display status for a stage."""
        if not stage.enabled:
            return "DISABLED"
        if stage.completed:
            return "[#B29EEE]\u2713 COMPLETED[/#B29EEE]"

        # Check for skipped conditions
        skip_conditions = {
            STAGE_STRUCT_FILTERS: DIR_DESCRIPTORS,
            STAGE_FINAL_DESCRIPTORS: DIR_STRUCT_FILTERS,
            STAGE_SYNTHESIS: DIR_STRUCT_FILTERS,
        }

        required_data = skip_conditions.get(stage.name)
        if required_data and not self.data_checker.check_stage_data(required_data):
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
                logger.warning(f"Could not copy config file for {key}: {copy_err}")

        logger.info(f"Saved run config snapshot to: {dest_dir}")
    except Exception as snapshot_err:
        logger.warning(f"Config snapshot failed: {snapshot_err}")


# Stage directory structure templates for README generation
_README_STAGE_SECTIONS = {
    STAGE_STRUCT_INI_FILTERS: """\
|   +-- 02_structural_filters_pre/      Pre-descriptors structural filters
|   |   +-- {filter_name}/             Per-filter results
|   |   |   +-- metrics.csv
|   |   |   +-- extended.csv
|   |   |   +-- filtered_molecules.csv
|   |   +-- filtered_molecules.csv     Combined passed molecules
|   |   +-- failed_molecules.csv       Combined failed molecules
|   |   +-- molecule_counts_comparison.png
|   |   +-- restriction_ratios_comparison.png
|
""",
    STAGE_DESCRIPTORS: """\
|   +-- 01_descriptors_initial/         Physicochemical descriptors
|   |   +-- metrics/
|   |   |   +-- descriptors_all.csv    All computed descriptors
|   |   |   +-- skipped_molecules.csv  Failed to parse
|   |   +-- filtered/
|   |   |   +-- filtered_molecules.csv Passed descriptor thresholds
|   |   |   +-- failed_molecules.csv   Failed descriptor thresholds
|   |   |   +-- descriptors_passed.csv Detailed metrics (passed)
|   |   |   +-- descriptors_failed.csv Detailed metrics (failed)
|   |   |   +-- pass_flags.csv         Pass/fail flags per descriptor
|   |   +-- plots/
|   |       +-- descriptors_distribution.png
|
""",
    STAGE_STRUCT_FILTERS: """\
|   +-- 03_structural_filters_post/     Post-descriptors structural filters
|   |   +-- {filter_name}/             Per-filter results (pains, brenk, nih, etc.)
|   |   |   +-- metrics.csv
|   |   |   +-- extended.csv
|   |   |   +-- filtered_molecules.csv
|   |   +-- filtered_molecules.csv     Combined passed molecules
|   |   +-- failed_molecules.csv       Combined failed molecules
|   |   +-- molecule_counts_comparison.png
|   |   +-- restriction_ratios_comparison.png
|
""",
    STAGE_SYNTHESIS: """\
|   +-- 04_synthesis/                   Retrosynthesis analysis
|   |   +-- synthesis_scores.csv       RAScore, SAScore, SCScore, SYBA
|   |   +-- synthesis_extended.csv     With retrosynthesis results
|   |   +-- filtered_molecules.csv     Synthesizable molecules
|   |   +-- input_smiles.smi           Input for AiZynthFinder
|   |   +-- retrosynthesis_results.json
|
""",
    STAGE_DOCKING: """\
|   +-- 05_docking/                     Molecular docking
|   |   +-- ligands.csv                Prepared ligands
|   |   +-- smina/
|   |   |   +-- poses/                 Docking poses (.pdbqt)
|   |   |   +-- scores.csv
|   |   +-- gnina/
|   |       +-- output.sdf
|   |       +-- scores.csv
|
""",
    STAGE_FINAL_DESCRIPTORS: """\
|   +-- 06_descriptors_final/           Final descriptor calculation
|       +-- metrics/
|       +-- filtered/
|       +-- plots/
|
""",
}


def _generate_structure_readme(
    base_path: Path, stages: list[PipelineStage], initial_count: int, final_count: int
) -> None:
    """Generate README.md documenting the output structure for this run."""
    try:
        readme_path = base_path / "README.md"
        enabled_stages = {s.name for s in stages if s.enabled}
        enabled_count = len(enabled_stages)
        completed_count = sum(1 for s in stages if s.completed)

        retention = (
            f"{100 * final_count / initial_count:.2f}%" if initial_count > 0 else "N/A"
        )

        # Build stage sections
        stage_sections = []
        for stage_name, section in _README_STAGE_SECTIONS.items():
            if stage_name in enabled_stages:
                stage_sections.append(section)

        stages_content = "".join(stage_sections)

        # Build stage execution summary
        stage_summary_lines = []
        for stage in stages:
            status = (
                "COMPLETED"
                if stage.completed
                else ("FAILED" if stage.enabled else "DISABLED")
            )
            stage_summary_lines.append(f"- {stage.name}: {status}")
        stage_summary = "\n".join(stage_summary_lines)

        content = f"""\
# HEDGEHOG Pipeline Output

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Summary

- Initial molecules: {initial_count}
- Final molecules: {final_count}
- Retention rate: {retention}
- Stages enabled: {enabled_count}
- Stages completed: {completed_count}

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
+-- configs/                            Configuration snapshots
|   +-- master_config_resolved.yml
|   +-- config_*.yml
|
+-- logs/                               Pipeline logs
    +-- pipeline_*.log
```

## File Naming Conventions

### Standard Output Files

- `filtered_molecules.csv` - Molecules that passed filters
- `failed_molecules.csv` - Molecules that failed filters
- `descriptors_all.csv` - All computed descriptors
- `metrics.csv` - Summary statistics
- `extended.csv` - Detailed results with all columns

### Column Structure

All CSV files use consistent column ordering:
1. `smiles` - SMILES string
2. `model_name` - Model/source identifier
3. `mol_idx` - Molecule index
4. Additional columns (stage-specific)

## Stage Execution Summary

{stage_summary}

## Notes

- This structure follows the new HEDGEHOG hierarchical organization
- Legacy flat structure is no longer used for new runs
- All paths use snake_case naming convention
- Stage numbering (01, 02, 03...) indicates execution order

## Related Documentation

See project repository for full documentation on:
- Configuration options
- Stage-specific parameters
- File format specifications
- Pipeline API reference
"""

        with open(readme_path, "w") as f:
            f.write(content)

        logger.info(f"Generated structure documentation: {readme_path}")
    except Exception as e:
        logger.warning(f"Failed to generate structure README: {e}")


def calculate_metrics(data, config: dict) -> bool:
    """Calculate metrics for molecular data using the configured pipeline.

    Args:
        data: Input molecular data (pandas DataFrame)
        config: Pipeline configuration dictionary

    Returns:
        True if all enabled stages completed successfully, False otherwise
    """
    try:
        _save_config_snapshot(config)
        pipeline = MolecularAnalysisPipeline(config)
        return pipeline.run_pipeline(data)
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}")
        return False
