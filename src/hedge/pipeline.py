import shutil
from pathlib import Path

import pandas as pd
import yaml

from hedge.configs.logger import load_config, logger
from hedge.stages.descriptors.main import main as descriptors_main
from hedge.stages.docking.utils import run_docking as docking_main
from hedge.stages.structFilters.main import main as structural_filters_main
from hedge.stages.synthesis.main import main as synthesis_main

# Constants - Directory names
DIR_DESCRIPTORS = "Descriptors"
DIR_STRUCT_FILTERS = "StructFilters"
DIR_BEFORE_DESCRIPTORS = "beforeDescriptors"
DIR_SYNTHESIS = "Synthesis"
DIR_DOCKING = "Docking"
DIR_FINAL_DESCRIPTORS = "finalDescriptors"
DIR_RUN_CONFIGS = "run_configs"

# Constants - File names
FILE_PASS_SMILES_TEMPLATE = "pass{stage}SMILES.csv"
FILE_MASTER_CONFIG = "master_config_resolved.yml"
FILE_GNINA_OUTPUT = "gnina_out.sdf"

# Constants - Stage names
STAGE_STRUCT_INI_FILTERS = "struct_ini_filters"
STAGE_DESCRIPTORS = "descriptors"
STAGE_STRUCT_FILTERS = "struct_filters"
STAGE_SYNTHESIS = "synthesis"
STAGE_DOCKING = "docking"
STAGE_FINAL_DESCRIPTORS = "final_descriptors"

# Constants - Config keys
CONFIG_DESCRIPTORS = "config_descriptors"
CONFIG_STRUCT_FILTERS = "config_structFilters"
CONFIG_SYNTHESIS = "config_synthesis"
CONFIG_DOCKING = "config_docking"
CONFIG_RUN_KEY = "run"
CONFIG_RUN_BEFORE_DESCRIPTORS = "run_before_descriptors"
CONFIG_TOOLS = "tools"
CONFIG_FOLDER_TO_SAVE = "folder_to_save"

# Constants - Command-line override keys
OVERRIDE_STRUCT_FILTERS_BEFORE = "_run_struct_filters_before_descriptors_override"
OVERRIDE_SINGLE_STAGE = "_run_single_stage_override"

# Constants - Docking tools
DOCKING_TOOL_SMINA = "smina"
DOCKING_TOOL_GNINA = "gnina"
DOCKING_TOOL_BOTH = "both"
DOCKING_RESULTS_DIR_TEMPLATE = {
    DOCKING_TOOL_SMINA: DIR_DOCKING + "/smina_results",
    DOCKING_TOOL_GNINA: DIR_DOCKING + "/gnina_results"
}


class PipelineStage:
    """Represents a single stage in the molecular analysis pipeline."""

    def __init__(self, name, config_key, directory):
        self.name = name
        self.config_key = config_key
        self.directory = directory
        self.enabled = False
        self.completed = False


class DataChecker:
    """Checks for existence of stage output data files."""

    def __init__(self, config):
        self.config = config
        self.base_path = Path(config[CONFIG_FOLDER_TO_SAVE])

    def check_stage_data(self, stage_name):
        """Check if data exists for a given stage."""
        try:
            stage_name = stage_name.strip()
            path = self._get_stage_output_path(stage_name)
            return path is not None and path.exists() and path.stat().st_size > 0
        except Exception:
            return False

    def _get_stage_output_path(self, stage_name):
        """Get the expected output file path for a stage."""
        if stage_name == DIR_DESCRIPTORS:
            return self.base_path / DIR_DESCRIPTORS / FILE_PASS_SMILES_TEMPLATE.format(stage=DIR_DESCRIPTORS)
        if stage_name == DIR_STRUCT_FILTERS:
            return self.base_path / DIR_STRUCT_FILTERS / FILE_PASS_SMILES_TEMPLATE.format(stage=DIR_STRUCT_FILTERS)
        if stage_name == DIR_BEFORE_DESCRIPTORS:
            return self.base_path / f"{DIR_BEFORE_DESCRIPTORS}_StructFilters" / FILE_PASS_SMILES_TEMPLATE.format(stage=DIR_STRUCT_FILTERS)
        if stage_name == DIR_SYNTHESIS:
            return self.base_path / DIR_SYNTHESIS / FILE_PASS_SMILES_TEMPLATE.format(stage=DIR_SYNTHESIS)
        return None


class PipelineStageRunner:
    """Executes individual pipeline stages and manages stage data flow."""

    DATA_SOURCE_PRIORITY = [DIR_SYNTHESIS, DIR_STRUCT_FILTERS, DIR_DESCRIPTORS, DIR_BEFORE_DESCRIPTORS]

    def __init__(self, config, data_checker):
        self.config = config
        self.data_checker = data_checker

    def find_latest_data_source(self):
        """Find the most recent stage with available output data."""
        for source in self.DATA_SOURCE_PRIORITY:
            if self.data_checker.check_stage_data(source):
                logger.debug(f"Found data from stage: {source}")
                return source

        logger.debug("No processed data found from any stage")
        return None

    def run_descriptors(self, data, subfolder=None):
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


    def run_structural_filters(self, prefix):
        """Run structural filters on molecules."""
        try:
            if prefix != DIR_BEFORE_DESCRIPTORS and not self.data_checker.check_stage_data(prefix):
                if self.config.get(OVERRIDE_SINGLE_STAGE) == STAGE_STRUCT_FILTERS:
                    logger.info("No previous stage data found, will use molecules from config")
                else:
                    logger.warning(f"No data available for structural filters with prefix: {prefix}")
                    return False

            config_struct_filters = load_config(self.config[CONFIG_STRUCT_FILTERS])
            if not config_struct_filters.get(CONFIG_RUN_KEY, False):
                logger.info("Structural filters disabled in config")
                return False

            structural_filters_main(self.config, prefix)
            return True
        except Exception as e:
            logger.error(f"Error running structural filters: {e}")
            return False


    def run_synthesis(self):
        """Run synthesis analysis."""
        try:
            latest_data_source = self.find_latest_data_source()
            if not latest_data_source:
                if self.config.get(OVERRIDE_SINGLE_STAGE) == STAGE_SYNTHESIS:
                    sampled_path = self.data_checker.base_path / "sampledMols.csv"
                    if sampled_path.exists():
                        logger.info("Using sampledMols.csv")
                    else:
                        logger.warning("No data available for synthesis, check `config.yml` file")
                        return False
                else:
                    logger.warning("No data available for synthesis, check provided path")
                    return False

            config_synthesis = load_config(self.config[CONFIG_SYNTHESIS])
            if not config_synthesis.get(CONFIG_RUN_KEY, False):
                logger.info("Synthesis disabled in config")
                return False

            synthesis_main(self.config)

            output_path = self.data_checker.base_path / DIR_SYNTHESIS / FILE_PASS_SMILES_TEMPLATE.format(stage=DIR_SYNTHESIS)
            if not output_path.exists():
                logger.error("Synthesis finished but no output file detected")
                return False
            return True
        except Exception as e:
            logger.error(f"Error running synthesis: {e}")
            return False

    def _parse_docking_tools(self, tools_cfg):
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

    def docking_results_present(self):
        """Check if docking results exist for all configured tools."""
        try:
            cfg = load_config(self.config[CONFIG_DOCKING])
        except Exception:
            return False

        base_folder = self.data_checker.base_path.resolve()
        tools_cfg = cfg.get(CONFIG_TOOLS, DOCKING_TOOL_BOTH)
        tools_list = self._parse_docking_tools(tools_cfg)

        present_map = {}

        if DOCKING_TOOL_GNINA in tools_list:
            gnina_dir = self._get_gnina_output_dir(cfg, base_folder)
            gnina_sdf = gnina_dir / FILE_GNINA_OUTPUT
            gnina_present = self._file_exists_and_not_empty(gnina_sdf)
            present_map[DOCKING_TOOL_GNINA] = gnina_present
            if not gnina_present:
                logger.debug(f"GNINA results not found at {gnina_sdf}")

        if DOCKING_TOOL_SMINA in tools_list:
            smina_dir = base_folder / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_SMINA]
            smina_present = self._directory_has_files(smina_dir)
            present_map[DOCKING_TOOL_SMINA] = smina_present
            if not smina_present:
                logger.debug(f"SMINA results not found in {smina_dir}")

        return any(present_map.get(tool, False) for tool in tools_list)

    def _get_gnina_output_dir(self, cfg, base_folder):
        """Get the GNINA output directory from config."""
        cfg_out_dir = cfg.get("gnina_output_dir")
        if cfg_out_dir:
            out_dir = Path(cfg_out_dir)
            return out_dir if out_dir.is_absolute() else base_folder / out_dir
        return base_folder / DOCKING_RESULTS_DIR_TEMPLATE[DOCKING_TOOL_GNINA]

    def _file_exists_and_not_empty(self, file_path):
        """Check if a file exists and is not empty."""
        try:
            return file_path.exists() and file_path.stat().st_size > 0
        except Exception:
            return False

    def _directory_has_files(self, dir_path):
        """Check if a directory exists and contains files."""
        try:
            return dir_path.exists() and any(p.is_file() for p in dir_path.iterdir())
        except Exception:
            return False

    def run_docking(self):
        """Run molecular docking."""
        try:
            config_docking = load_config(self.config[CONFIG_DOCKING])
            if not config_docking.get(CONFIG_RUN_KEY, False):
                logger.info("Docking disabled in config")
                return False

            success = bool(docking_main(self.config))
            if not success:
                return False

            if not self.docking_results_present():
                logger.error("Docking finished but no results detected in output directories")
                return False
            return True
        except Exception as e:
            logger.error(f"Error running docking: {e}")
            return False


class MolecularAnalysisPipeline:
    """Orchestrates the execution of the molecular analysis pipeline."""

    def __init__(self, config):
        self.config = config
        self.data_checker = DataChecker(config)
        self.stage_runner = PipelineStageRunner(config, self.data_checker)
        self.current_data = None

        self.stages = [PipelineStage(STAGE_STRUCT_INI_FILTERS, CONFIG_STRUCT_FILTERS, DIR_BEFORE_DESCRIPTORS),
                       PipelineStage(STAGE_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_DESCRIPTORS),
                       PipelineStage(STAGE_STRUCT_FILTERS, CONFIG_STRUCT_FILTERS, DIR_STRUCT_FILTERS),
                       PipelineStage(STAGE_SYNTHESIS, CONFIG_SYNTHESIS, DIR_SYNTHESIS),
                       PipelineStage(STAGE_DOCKING, CONFIG_DOCKING, DIR_DOCKING),
                       PipelineStage(STAGE_FINAL_DESCRIPTORS, CONFIG_DESCRIPTORS, DIR_FINAL_DESCRIPTORS)
                      ]

        self._initialize_stages()

    def _initialize_stages(self):
        """Initialize stages by loading their configs and determining enabled status."""
        single_stage_override = self.config.get(OVERRIDE_SINGLE_STAGE)
        for stage in self.stages:
            try:
                stage_config = load_config(self.config[stage.config_key])
                stage.enabled = stage_config.get(CONFIG_RUN_KEY, False)

                if stage.name == STAGE_STRUCT_INI_FILTERS:
                    if OVERRIDE_STRUCT_FILTERS_BEFORE in self.config:
                        run_before = self.config[OVERRIDE_STRUCT_FILTERS_BEFORE]
                    else:
                        run_before = stage_config.get(CONFIG_RUN_BEFORE_DESCRIPTORS, False)

                    stage.enabled = stage.enabled and run_before

                if single_stage_override:
                    if stage.name == single_stage_override:
                        stage.enabled = True
                        logger.info(f"Single stage mode: enabling only {stage.name}")
                    else:
                        stage.enabled = False

                logger.debug(f"Stage {stage.name}: {'Enabled' if stage.enabled else 'Disabled'}")
            except Exception as e:
                logger.warning(f"Could not load config for {stage.name}: {e}")
                stage.enabled = False

    def get_latest_data(self, skip_descriptors=False):
        """Load the most recent stage output data.
        
        Args:
            skip_descriptors: If True, skip Descriptors stage when looking for data source
        """
        priority = self.stage_runner.DATA_SOURCE_PRIORITY.copy()
        if skip_descriptors:
            priority = [p for p in priority if p != DIR_DESCRIPTORS]

        latest_source = None
        for source in priority:
            if self.data_checker.check_stage_data(source):
                latest_source = source
                logger.debug(f"Found data from stage: {source}")
                break

        if latest_source:
            try:
                filename = FILE_PASS_SMILES_TEMPLATE.format(stage=latest_source)
                path = self.data_checker.base_path / latest_source / filename
                data = pd.read_csv(path)
                logger.info(f"Loaded latest data from {latest_source}: {len(data)} molecules")

                if len(data) == 0:
                    logger.warning(f"No molecules in {latest_source}, trying previous step...")
                    current_idx = priority.index(latest_source)
                    for next_source in priority[current_idx + 1:]:
                        if self.data_checker.check_stage_data(next_source):
                            filename = FILE_PASS_SMILES_TEMPLATE.format(stage=next_source)
                            path = self.data_checker.base_path / next_source / filename
                            data = pd.read_csv(path)
                            if len(data) > 0:
                                logger.info(f"Loaded data from {next_source}: {len(data)} molecules (previous step)")
                                return data
                    logger.warning("All checked data sources are empty")
                    return data

                return data
            except Exception as e:
                logger.warning(f"Could not load data from {latest_source}: {e}")

        return self.current_data

    def run_pipeline(self, data):
        """Execute the full molecular analysis pipeline."""
        self.current_data = data
        success_count = 0
        total_enabled_stages = sum(1 for stage in self.stages if stage.enabled)

        # Stage 1': Pre-descriptors structural filters
        if self.stages[0].enabled:
            logger.info("---------> [#B29EEE]Stage 1': Pre-descriptors Structural Filters[/#B29EEE]")
            if self.stage_runner.run_structural_filters(DIR_BEFORE_DESCRIPTORS):
                self.stages[0].completed = True
                success_count += 1

        # Stage 1: Descriptors calculation
        if self.stages[1].enabled:
            logger.info("---------> [#B29EEE]Stage 1: Molecular Descriptors[/#B29EEE]")
            if self.stage_runner.run_descriptors(data):
                self.stages[1].completed = True
                success_count += 1

        # Stage 2: Post-descriptors structural filters
        if self.stages[2].enabled:
            logger.info("---------> [#B29EEE]Stage 2: Post-descriptors Structural Filters[/#B29EEE]")
            if self.stage_runner.run_structural_filters(DIR_DESCRIPTORS):
                self.stages[2].completed = True
                success_count += 1
            else:
                logger.info("No molecules left after descriptors; ending pipeline early.")
                self._log_pipeline_summary()
                initial_count = len(data)
                final_data = self.get_latest_data(skip_descriptors=True)
                final_count = len(final_data) if final_data is not None else 0
                logger.info("---------> Molecule Count Summary")
                logger.info(f"Molecules before filtration: {initial_count}")
                logger.info(f"Molecules remaining after all stages: {final_count}")
                if initial_count > 0:
                    logger.info(f"  Retention rate: {100*final_count/initial_count:.2f}%")
                if final_data is not None and len(final_data) > 0:
                    final_output_path = self.data_checker.base_path / "finalMolecules.csv"
                    final_data.to_csv(final_output_path, index=False)
                    logger.info(f"Saved {final_count} final molecules to {final_output_path}")
                return success_count == total_enabled_stages

        # Stage 3: Synthesis
        if self.stages[3].enabled:
            logger.info("---------> [#B29EEE]Stage 3: Synthesis Analysis[/#B29EEE]")
            if self.stage_runner.run_synthesis():
                self.stages[3].completed = True
                success_count += 1
            else:
                output_path = self.data_checker.base_path / DIR_SYNTHESIS / FILE_PASS_SMILES_TEMPLATE.format(stage=DIR_SYNTHESIS)
                if output_path.exists():
                    try:
                        df_check = pd.read_csv(output_path)
                        if len(df_check) == 0:
                            logger.info("No molecules left after synthesis; ending pipeline early.")
                            self._log_pipeline_summary()
                            initial_count = len(data)
                            final_data = self.get_latest_data(skip_descriptors=True)
                            final_count = len(final_data) if final_data is not None else 0
                            logger.info("---------> Molecule Count Summary")
                            logger.info(f"Molecules before filtration: {initial_count}")
                            logger.info(f"Molecules remaining after all stages: {final_count}")
                            if initial_count > 0:
                                logger.info(f"Retention rate: {100*final_count/initial_count:.2f}%")
                            if final_data is not None and len(final_data) > 0:
                                final_output_path = self.data_checker.base_path / "finalMolecules.csv"
                                final_data.to_csv(final_output_path, index=False)
                                logger.info(f"Saved {final_count} final molecules to {final_output_path}")
                            return success_count == total_enabled_stages
                        logger.warning("Synthesis stage failed but output file exists. Continuing with available molecules.")
                    except Exception:
                        logger.warning("Synthesis stage failed. Check logs for details.")
                else:
                    logger.error("Synthesis stage failed (no output file created). Check logs above for error details.")
                    logger.info("Continuing pipeline without synthesis results...")

        # Stage 4: Docking
        if self.stages[4].enabled:
            logger.info("---------> [#B29EEE]Stage 4: Molecular Docking[/#B29EEE]")
            if self.stage_runner.run_docking():
                self.stages[4].completed = True
                success_count += 1

        # Stage 5': Final descriptors calculation
        if self.stages[5].enabled:
            logger.info("---------> [#B29EEE]Stage 5': Final Descriptors Calculation[/#B29EEE]")
            final_data = self.get_latest_data(skip_descriptors=True)
            if final_data is not None and len(final_data) > 0:
                if self.stage_runner.run_descriptors(final_data, subfolder=DIR_FINAL_DESCRIPTORS):
                    self.stages[5].completed = True
                    success_count += 1
            else:
                logger.info("No molecules available for final descriptors calculation (skipping descriptors stage as source)")
                if final_data is not None and len(final_data) == 0:
                    logger.info("Ending pipeline: no molecules from previous steps (excluding descriptors)")

        logger.info(f"Pipeline completed: {success_count}/{total_enabled_stages} stages successful")
        self._log_pipeline_summary()

        initial_count = len(data)
        final_data = self.get_latest_data(skip_descriptors=True)
        final_count = len(final_data) if final_data is not None else 0
        logger.info("---------> Molecule Count Summary")
        logger.info(f"Molecules before filtration: {initial_count}")
        logger.info(f"Molecules remaining after all stages: {final_count}")
        if initial_count > 0:
            logger.info(f"Retention rate: {100*final_count/initial_count:.2f}%")

        if final_data is not None and len(final_data) > 0:
            final_output_path = self.data_checker.base_path / "finalMolecules.csv"
            final_data.to_csv(final_output_path, index=False)
            logger.info(f"Saved {final_count} final molecules to {final_output_path}")

        return success_count == total_enabled_stages

    def _log_pipeline_summary(self):
        """Log a summary of pipeline execution status."""
        logger.info("---------> [#B29EEE]Pipeline Execution Summary[/#B29EEE]")

        for stage in self.stages:
            if stage.enabled:
                if stage.completed:
                    status = "[#B29EEE]✓ COMPLETED[/#B29EEE]"
                elif (stage.name == STAGE_STRUCT_FILTERS and not self.data_checker.check_stage_data(DIR_DESCRIPTORS)) or (stage.name == STAGE_FINAL_DESCRIPTORS and not self.data_checker.check_stage_data(DIR_STRUCT_FILTERS)) or (stage.name == STAGE_SYNTHESIS and not self.data_checker.check_stage_data(DIR_STRUCT_FILTERS)):
                    status = "[yellow]⟂ SKIPPED (no molecules)[/yellow]"
                else:
                    status = "[red]✗ FAILED[/red]"

                logger.info(f"{stage.name}: {status}")
            else:
                logger.info(f"{stage.name}: DISABLED")


def _save_config_snapshot(config):
    """Save a snapshot of configuration files for provenance."""
    try:
        base_path = Path(config[CONFIG_FOLDER_TO_SAVE])
        dest_dir = base_path / DIR_RUN_CONFIGS
        dest_dir.mkdir(parents=True, exist_ok=True)

        master_config_path = dest_dir / FILE_MASTER_CONFIG
        with open(master_config_path, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)

        config_keys = [k for k in config.keys() if k.startswith("config_")]
        for key in config_keys:
            path_str = config.get(key)
            if path_str:
                try:
                    src_path = Path(path_str)
                    if src_path.exists():
                        shutil.copyfile(src_path, dest_dir / src_path.name)
                except Exception as copy_err:
                    logger.warning(f"Could not copy config file for {key}: {copy_err}")

        logger.info(f"Saved run config snapshot to: {dest_dir}")
    except Exception as snapshot_err:
        logger.warning(f"Config snapshot failed: {snapshot_err}")


def calculate_metrics(data, config):
    """
    Calculate metrics for molecular data using the configured pipeline.
    
    Args:
        data: Input molecular data
        config: Pipeline configuration dictionary
    
    Returns
    -------
        True if all enabled stages completed successfully, False otherwise
    """
    try:
        _save_config_snapshot(config)

        pipeline = MolecularAnalysisPipeline(config)
        return pipeline.run_pipeline(data)
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}")
        return False
