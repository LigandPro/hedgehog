"""Report generator for HEDGEHOG pipeline results."""

import base64
import json
import logging
import statistics
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, PackageLoader

from hedgehog.reporting import plots

logger = logging.getLogger(__name__)

# Stage directory names
STAGE_DIRS = {
    "descriptors_initial": "stages/01_descriptors_initial",
    "struct_filters_pre": "stages/02_structural_filters_pre",
    "struct_filters_post": "stages/03_structural_filters_post",
    "synthesis": "stages/04_synthesis",
    "docking": "stages/05_docking",
    "docking_filters": "stages/06_docking_filters",
    "descriptors_final": "stages/07_descriptors_final",
}

# Stage display names
STAGE_DISPLAY_NAMES = {
    "descriptors_initial": "Initial Descriptors",
    "struct_filters_pre": "Pre-Descriptors Filters",
    "struct_filters_post": "Post-Descriptors Filters",
    "synthesis": "Synthesis Analysis",
    "docking": "Molecular Docking",
    "docking_filters": "Docking Filters",
    "descriptors_final": "Final Descriptors",
}

# Key descriptors to show in report
KEY_DESCRIPTORS = [
    "MolWt",
    "LogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
]

# Mapping for case-insensitive descriptor lookup (lowercase -> display name)
DESCRIPTOR_ALIASES = {
    "molwt": "MolWt",
    "logp": "LogP",
    "clogp": "cLogP",
    "tpsa": "TPSA",
    "numhdonors": "NumHDonors",
    "numhacceptors": "NumHAcceptors",
    "numrotatablebonds": "NumRotatableBonds",
    "hbd": "NumHDonors",
    "hba": "NumHAcceptors",
    "qed": "QED",
    "fsp3": "Fsp3",
    "n_atoms": "n_atoms",
    "n_heavy_atoms": "n_heavy_atoms",
    "n_het_atoms": "n_het_atoms",
    "n_n_atoms": "n_N_atoms",
    "fn_atoms": "fN_atoms",
    "charged_mol": "charged_mol",
    "sw": "SW",
    "ring_size": "ring_size",
    "n_rings": "n_rings",
    "n_aroma_rings": "n_aroma_rings",
    "n_fused_aromatic_rings": "n_fused_aromatic_rings",
    "n_rigid_bonds": "n_rigid_bonds",
    "n_rot_bonds": "n_rot_bonds",
}

# Columns to exclude from descriptor analysis (not numeric descriptors)
DESCRIPTOR_EXCLUDE_COLS = {"smiles", "model_name", "mol_idx", "chars", "name", "id"}


class ReportGenerator:
    """Generates comprehensive HTML reports for HEDGEHOG pipeline runs."""

    def __init__(
        self,
        base_path: Path,
        stages: list[Any],
        config: dict[str, Any],
        initial_count: int,
        final_count: int,
    ):
        """Initialize the report generator.

        Args:
            base_path: Base path of the pipeline run output
            stages: List of PipelineStage objects
            config: Pipeline configuration dictionary
            initial_count: Number of molecules at pipeline start
            final_count: Number of molecules at pipeline end
        """
        self.base_path = Path(base_path)
        self.stages = stages
        self.config = config
        self.initial_count = initial_count
        self.final_count = final_count
        self.output_dir = self.base_path
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate(self) -> Path:
        """Generate the full HTML report.

        Returns:
            Path to the generated report file
        """
        logger.info("Generating pipeline report...")

        # Collect all data
        data = self._collect_data()

        # Generate plots
        plot_htmls = self._generate_plots(data)

        # Render HTML
        html_content = self._render_template(data, plot_htmls)

        # Save report
        report_path = self._save_report(html_content)

        # Save JSON data
        self._save_json_data(data)

        logger.info("Report generated: %s", report_path)
        return report_path

    def _collect_data(self) -> dict[str, Any]:
        """Collect all metrics from stage outputs.

        Returns:
            Dictionary with all collected data
        """
        funnel_data = self._get_all_funnel_data()
        return {
            "metadata": self._get_metadata(),
            "summary": self._get_summary(),
            "funnel": funnel_data["all"],
            "funnel_by_model": funnel_data["by_model"],
            "available_models": funnel_data["models"],
            "stages": self._get_stage_stats(),
            "models": self._get_model_stats(),
            "descriptors": self._get_descriptor_stats(),
            "descriptors_detailed": self._get_descriptors_detailed(),
            "filters": self._get_filter_stats(),
            "filters_detailed": self._get_filters_detailed(),
            "synthesis": self._get_synthesis_stats(),
            "synthesis_detailed": self._get_synthesis_detailed(),
            "retrosynthesis": self._get_retrosynthesis_detailed(),
            "docking": self._get_docking_stats(),
            "docking_detailed": self._get_docking_detailed(),
            "docking_filters_detailed": self._get_docking_filters_detailed(),
            "descriptors_final": self._get_descriptor_stats("descriptors_final"),
            "descriptors_final_detailed": self._get_descriptors_detailed(
                "descriptors_final"
            ),
            "existing_plots": self._get_existing_plots(),
            "config": self._get_config_summary(),
        }

    def _get_metadata(self) -> dict[str, Any]:
        """Get report metadata."""
        return {
            "generated_at": datetime.now().isoformat(),
            "hedgehog_version": "1.0.0",
            "run_path": str(self.base_path),
        }

    def _get_summary(self) -> dict[str, Any]:
        """Get executive summary statistics."""
        retention_rate = (
            self.final_count / self.initial_count if self.initial_count > 0 else 0
        )

        # If stages not provided, detect from directory structure
        if self.stages:
            enabled_stages = [s for s in self.stages if s.enabled]
            completed_stages = [s for s in self.stages if s.completed]
            stage_statuses = [
                {
                    "name": s.name,
                    "enabled": s.enabled,
                    "completed": s.completed,
                    "status": "completed"
                    if s.completed
                    else ("failed" if s.enabled else "disabled"),
                }
                for s in self.stages
            ]
        else:
            # Auto-detect stages from existing directories
            stage_order = [
                ("descriptors_initial", "Initial Descriptors"),
                ("struct_filters_pre", "Structural Filters"),
                ("synthesis", "Synthesis Analysis"),
                ("docking", "Molecular Docking"),
                ("docking_filters", "Docking Filters"),
                ("descriptors_final", "Final Descriptors"),
            ]
            enabled_stages = []
            completed_stages = []
            stage_statuses = []
            for stage_key, display_name in stage_order:
                stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
                if stage_dir.exists():
                    enabled_stages.append(stage_key)
                    # Check if stage has output files (completed)
                    has_output = any(stage_dir.glob("*.csv")) or any(
                        stage_dir.glob("*/*.csv")
                    )
                    if has_output:
                        completed_stages.append(stage_key)
                    stage_statuses.append(
                        {
                            "name": display_name,
                            "enabled": True,
                            "completed": has_output,
                            "status": "completed" if has_output else "failed",
                        }
                    )

        return {
            "initial_molecules": self.initial_count,
            "final_molecules": self.final_count,
            "retention_rate": retention_rate,
            "retention_percent": f"{retention_rate * 100:.2f}%",
            "stages_enabled": len(enabled_stages),
            "stages_completed": len(completed_stages),
            "stage_statuses": stage_statuses,
        }

    def _get_funnel_data(self) -> list[dict[str, Any]]:
        """Get molecule funnel data through pipeline stages."""
        funnel = [{"stage": "Initial", "count": self.initial_count}]

        # Pipeline order: Pre-Filters → Descriptors → Post-Filters → Synthesis → Docking → Docking Filters → Final Descriptors
        stage_order = [
            ("struct_filters_pre", "Pre-Filters"),
            ("descriptors_initial", "Descriptors"),
            ("struct_filters_post", "Post-Filters"),
            ("synthesis", "Synthesis"),
            ("docking", "Docking"),
            ("docking_filters", "Docking Filters"),
            ("descriptors_final", "Final Descriptors"),
        ]

        for stage_key, display_name in stage_order:
            count = self._get_stage_output_count(stage_key)
            if count is not None:
                funnel.append({"stage": display_name, "count": count})

        # Add final count if different from last stage
        if funnel and funnel[-1]["count"] != self.final_count:
            funnel.append({"stage": "Final", "count": self.final_count})

        return funnel

    def _get_stage_output_count(self, stage_key: str) -> int | None:
        """Get molecule count from stage output file."""
        stage_dir = STAGE_DIRS.get(stage_key)
        if not stage_dir:
            return None

        # Try common output file locations
        paths_to_try = [
            self.base_path / stage_dir / "filtered_molecules.csv",
            self.base_path / stage_dir / "filtered" / "filtered_molecules.csv",
            self.base_path / stage_dir / "ligands.csv",  # For docking stage
        ]

        for path in paths_to_try:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    return len(df)
                except Exception as e:
                    logger.debug("Could not read %s: %s", path, e)
                    continue
        return None

    def _get_stage_output_count_by_model(
        self, stage_key: str, model_name: str | None = None
    ) -> int | None:
        """Get molecule count from stage output file, optionally filtered by model.

        Args:
            stage_key: Stage key from STAGE_DIRS
            model_name: Optional model name to filter by

        Returns:
            Count of molecules, or None if data unavailable
        """
        stage_dir = STAGE_DIRS.get(stage_key)
        if not stage_dir:
            return None

        paths_to_try = [
            self.base_path / stage_dir / "filtered_molecules.csv",
            self.base_path / stage_dir / "filtered" / "filtered_molecules.csv",
            self.base_path / stage_dir / "ligands.csv",
        ]

        for path in paths_to_try:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    if model_name and "model_name" in df.columns:
                        df = df[df["model_name"] == model_name]
                    return len(df)
                except Exception as e:
                    logger.debug("Could not read %s: %s", path, e)
                    continue
        return None

    def _get_available_models(self) -> list[str]:
        """Get list of available model names from input or output files.

        Returns:
            List of unique model names
        """
        models = set()

        # Try input file first
        input_path = self.base_path / "input" / "sampled_molecules.csv"
        if input_path.exists():
            try:
                df = pd.read_csv(input_path)
                if "model_name" in df.columns:
                    models.update(df["model_name"].dropna().unique())
            except Exception as e:
                logger.debug("Could not read %s: %s", input_path, e)

        # Try output file
        output_path = self.output_dir / "final_molecules.csv"
        if output_path.exists():
            try:
                df = pd.read_csv(output_path)
                if "model_name" in df.columns:
                    models.update(df["model_name"].dropna().unique())
            except Exception as e:
                logger.debug("Could not read %s: %s", output_path, e)

        return sorted(models)

    def _get_initial_count_by_model(self, model_name: str | None = None) -> int:
        """Get initial molecule count, optionally filtered by model.

        Args:
            model_name: Optional model name to filter by

        Returns:
            Count of initial molecules
        """
        if not model_name:
            return self.initial_count

        input_path = self.base_path / "input" / "sampled_molecules.csv"
        if input_path.exists():
            try:
                df = pd.read_csv(input_path)
                if "model_name" in df.columns:
                    return len(df[df["model_name"] == model_name])
            except Exception as e:
                logger.debug("Could not read %s: %s", input_path, e)
        return 0

    def _get_funnel_data_by_model(
        self, model_name: str | None = None
    ) -> list[dict[str, Any]]:
        """Get molecule funnel data for a specific model or all models.

        Args:
            model_name: Model name to filter by, or None for all models

        Returns:
            List of funnel stage data
        """
        initial_count = self._get_initial_count_by_model(model_name)
        funnel = [{"stage": "Initial", "count": initial_count}]

        stage_order = [
            ("struct_filters_pre", "Pre-Filters"),
            ("descriptors_initial", "Descriptors"),
            ("struct_filters_post", "Post-Filters"),
            ("synthesis", "Synthesis"),
            ("docking", "Docking"),
            ("docking_filters", "Docking Filters"),
            ("descriptors_final", "Final Descriptors"),
        ]

        for stage_key, display_name in stage_order:
            count = self._get_stage_output_count_by_model(stage_key, model_name)
            if count is not None:
                funnel.append({"stage": display_name, "count": count})

        return funnel

    def _get_all_funnel_data(self) -> dict[str, Any]:
        """Get funnel data for all models and per-model breakdown.

        Returns:
            Dictionary with 'all' funnel and 'by_model' funnel data
        """
        result = {
            "all": self._get_funnel_data(),
            "by_model": {},
            "models": self._get_available_models(),
        }

        for model in result["models"]:
            result["by_model"][model] = self._get_funnel_data_by_model(model)

        return result

    def _get_stage_stats(self) -> list[dict[str, Any]]:
        """Get pass/fail statistics for each stage."""
        stats = []

        for stage_key, display_name in STAGE_DISPLAY_NAMES.items():
            stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
            if not stage_dir.exists():
                continue

            passed = 0
            failed = 0

            # Try to read filtered and failed molecules
            passed_path = stage_dir / "filtered_molecules.csv"
            if not passed_path.exists():
                passed_path = stage_dir / "filtered" / "filtered_molecules.csv"

            failed_path = stage_dir / "failed_molecules.csv"
            if not failed_path.exists():
                failed_path = stage_dir / "filtered" / "failed_molecules.csv"

            if passed_path.exists():
                try:
                    passed = len(pd.read_csv(passed_path))
                except Exception:
                    pass

            if failed_path.exists():
                try:
                    failed = len(pd.read_csv(failed_path))
                except Exception:
                    pass

            if passed > 0 or failed > 0:
                stats.append(
                    {
                        "stage": display_name,
                        "passed": passed,
                        "failed": failed,
                        "total": passed + failed,
                    }
                )

        return stats

    def _get_model_stats(self) -> list[dict[str, Any]]:
        """Get per-model statistics."""
        # Try to load final molecules to get model breakdown
        final_path = self.output_dir / "final_molecules.csv"
        if not final_path.exists():
            return []

        try:
            df = pd.read_csv(final_path)
        except Exception:
            return []

        if "model_name" not in df.columns:
            return []

        model_stats = []
        for model in df["model_name"].unique():
            final_count = len(df[df["model_name"] == model])

            # Try to get initial count from input
            initial_count = self._get_initial_model_count(model)

            model_stats.append(
                {
                    "model_name": model,
                    "initial": initial_count or final_count,
                    "final": final_count,
                    "losses": self._get_model_losses(model),
                }
            )

        return model_stats

    def _get_initial_model_count(self, model: str) -> int | None:
        """Get initial molecule count for a model."""
        input_path = self.base_path / "input" / "sampled_molecules.csv"
        if not input_path.exists():
            return None

        try:
            df = pd.read_csv(input_path)
            if "model_name" in df.columns:
                return len(df[df["model_name"] == model])
        except Exception:
            pass
        return None

    def _get_model_losses(self, model: str) -> dict[str, int]:
        """Get molecule losses per stage for a model."""
        losses = {}

        stage_pairs = [
            ("descriptors", "descriptors_initial"),
            ("struct_filters", "struct_filters_post"),
            ("synthesis", "synthesis"),
            ("docking", "docking"),
            ("docking_filters", "docking_filters"),
        ]

        for loss_key, stage_key in stage_pairs:
            failed_path = (
                self.base_path / STAGE_DIRS.get(stage_key, "") / "failed_molecules.csv"
            )
            if not failed_path.exists():
                failed_path = (
                    self.base_path
                    / STAGE_DIRS.get(stage_key, "")
                    / "filtered"
                    / "failed_molecules.csv"
                )

            if failed_path.exists():
                try:
                    df = pd.read_csv(failed_path)
                    if "model_name" in df.columns:
                        losses[loss_key] = len(df[df["model_name"] == model])
                    else:
                        losses[loss_key] = 0
                except Exception:
                    losses[loss_key] = 0
            else:
                losses[loss_key] = 0

        return losses

    def _get_descriptor_stats(
        self, stage_key: str = "descriptors_initial"
    ) -> dict[str, Any]:
        """Get descriptor statistics.

        Args:
            stage_key: Stage key from STAGE_DIRS (default: descriptors_initial)
        """
        desc_dir = self.base_path / STAGE_DIRS[stage_key]

        # Try multiple possible locations
        paths_to_try = [
            desc_dir / "metrics" / "descriptors_all.csv",
            desc_dir / "filtered" / "descriptors_passed.csv",
        ]

        df = None
        for path in paths_to_try:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    break
                except Exception:
                    continue

        if df is None:
            return {}

        stats = {"distributions": {}, "summary": {}}

        # Build case-insensitive column mapping
        col_map = self._build_descriptor_column_map(df.columns)

        for desc in KEY_DESCRIPTORS:
            col = col_map.get(desc)
            if col and col in df.columns:
                values = df[col].dropna().tolist()
                if values:
                    stats["distributions"][desc] = values
                    stats["summary"][desc] = {
                        "mean": float(df[col].mean()),
                        "std": float(df[col].std()),
                        "min": float(df[col].min()),
                        "max": float(df[col].max()),
                    }

        return stats

    def _build_descriptor_column_map(self, columns: list[str]) -> dict[str, str]:
        """Build mapping from standard descriptor names to actual column names.

        Args:
            columns: List of column names in the dataframe

        Returns:
            Dict mapping standard name -> actual column name
        """
        col_map = {}
        col_lower = {c.lower(): c for c in columns}

        for desc in KEY_DESCRIPTORS + ["QED", "Fsp3"]:
            # Direct match
            if desc in columns:
                col_map[desc] = desc
            # Case-insensitive match
            elif desc.lower() in col_lower:
                col_map[desc] = col_lower[desc.lower()]

        # Check aliases
        for alias, standard in DESCRIPTOR_ALIASES.items():
            if alias in col_lower and standard not in col_map:
                col_map[standard] = col_lower[alias]

        return col_map

    def _get_filter_stats(self) -> dict[str, Any]:
        """Get filter statistics."""
        filter_stats = {"by_filter": {}, "totals": {}, "by_model": {}}

        # Check both pre and post filter directories
        for stage_key in ["struct_filters_pre", "struct_filters_post"]:
            stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
            if not stage_dir.exists():
                continue

            # Look for per-filter subdirectories
            for subdir in stage_dir.iterdir():
                if not subdir.is_dir():
                    continue
                if subdir.name in ["plots", "logs"]:
                    continue

                metrics_path = subdir / "metrics.csv"
                if metrics_path.exists():
                    try:
                        df = pd.read_csv(metrics_path)
                        if "failed_count" in df.columns:
                            total_failed = df["failed_count"].sum()
                            filter_stats["totals"][subdir.name] = int(total_failed)

                            # Get per-model breakdown if available
                            if "model_name" in df.columns:
                                model_counts = (
                                    df.groupby("model_name")["failed_count"]
                                    .sum()
                                    .to_dict()
                                )
                                filter_stats["by_filter"][subdir.name] = {
                                    str(k): int(v) for k, v in model_counts.items()
                                }
                    except Exception:
                        continue

        return filter_stats

    def _get_synthesis_stats(self) -> dict[str, Any]:
        """Get synthesis analysis statistics."""
        synth_dir = self.base_path / STAGE_DIRS["synthesis"]
        scores_path = synth_dir / "synthesis_scores.csv"

        if not scores_path.exists():
            return {}

        try:
            df = pd.read_csv(scores_path)
        except Exception:
            return {}

        stats = {"distributions": {}, "scatter_data": {}}

        # Get score distributions
        score_columns = ["sa_score", "syba_score", "ra_score", "sc_score"]
        for col in score_columns:
            if col in df.columns:
                values = df[col].dropna().tolist()
                if values:
                    stats["distributions"][col] = values

        # Get scatter plot data
        if "sa_score" in df.columns and "syba_score" in df.columns:
            stats["scatter_data"] = {
                "sa_scores": df["sa_score"].dropna().tolist(),
                "syba_scores": df["syba_score"].dropna().tolist(),
                "model_names": df["model_name"].tolist()
                if "model_name" in df.columns
                else [],
            }

        return stats

    def _get_model_name_lookup(self) -> dict[str, str]:
        """Build a lookup table from molecule name to model name.

        Returns:
            Dictionary mapping molecule name/ID to model name
        """
        lookup = {}
        ligands_csv = self.base_path / STAGE_DIRS["docking"] / "ligands.csv"
        if ligands_csv.exists():
            try:
                df = pd.read_csv(ligands_csv)
                if "model_name" in df.columns:
                    for id_col in ["name", "mol_idx", "molecule_id"]:
                        if id_col in df.columns:
                            for _, row in df.iterrows():
                                if pd.notna(row.get(id_col)) and pd.notna(
                                    row.get("model_name")
                                ):
                                    lookup[str(row[id_col])] = row["model_name"]
                            break
            except Exception:
                pass
        return lookup

    def _detect_docking_tools(self) -> list[str]:
        """Detect available docking tools from directory structure.

        Scans the docking directory for subdirectories and SDF files to
        identify which docking tools were used.

        Returns:
            List of tool names (e.g., ["gnina", "smina"])
        """
        docking_dir = self.base_path / STAGE_DIRS["docking"]
        if not docking_dir.exists():
            return []

        tools = set()

        # Check for subdirectories with SDF results
        for subdir in docking_dir.iterdir():
            if subdir.is_dir():
                # Check for *_out.sdf or output.sdf in subdirectory
                for sdf_file in subdir.glob("*.sdf"):
                    if "out" in sdf_file.name.lower():
                        tools.add(subdir.name)
                        break

        # Check for *_out.sdf files in main docking directory
        for sdf_file in docking_dir.glob("*_out.sdf"):
            tool_name = sdf_file.stem.replace("_out", "")
            tools.add(tool_name)

        # Check for config files (gnina_config.ini, smina_config.ini)
        for config_file in docking_dir.glob("*_config.ini"):
            tool_name = config_file.stem.replace("_config", "")
            tools.add(tool_name)

        return sorted(tools)

    def _parse_docking_sdf(
        self, sdf_path: Path, model_lookup: dict[str, str] | None = None
    ) -> list[dict[str, Any]]:
        """Parse SDF file to extract docking scores and molecule info.

        Args:
            sdf_path: Path to SDF file with docking results
            model_lookup: Optional dict mapping molecule names to model names

        Returns:
            List of dicts with molecule_id, affinity, model_name, and any extra scores
        """
        try:
            from rdkit import Chem
        except ImportError:
            logger.warning("RDKit not available for SDF parsing")
            return []

        if not sdf_path.exists():
            return []

        results = []
        mol_best_scores: dict[str, dict] = {}  # Track best score per molecule

        try:
            supplier = Chem.SDMolSupplier(str(sdf_path))
            for mol in supplier:
                if mol is None:
                    continue

                # Get molecule name
                mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
                if not mol_name:
                    continue

                # Get affinity score
                affinity = None
                for prop_name in ["minimizedAffinity", "affinity", "score"]:
                    if mol.HasProp(prop_name):
                        try:
                            affinity = float(mol.GetProp(prop_name))
                            break
                        except (ValueError, TypeError):
                            pass

                if affinity is None:
                    continue

                # Get model name from SDF property or lookup
                model_name = "Unknown"
                if mol.HasProp("s_sm_model_name"):
                    model_name = mol.GetProp("s_sm_model_name")
                elif model_lookup and mol_name in model_lookup:
                    model_name = model_lookup[mol_name]

                # Collect additional scores (GNINA CNN scores)
                extra_scores = {}
                for prop_name in ["CNNscore", "CNNaffinity", "CNN_VS"]:
                    if mol.HasProp(prop_name):
                        try:
                            extra_scores[prop_name] = float(mol.GetProp(prop_name))
                        except (ValueError, TypeError):
                            pass

                record = {
                    "molecule_id": mol_name,
                    "affinity": affinity,
                    "model_name": model_name,
                    **extra_scores,
                }

                # Keep only best score per molecule (lowest affinity)
                if mol_name not in mol_best_scores:
                    mol_best_scores[mol_name] = record
                elif affinity < mol_best_scores[mol_name]["affinity"]:
                    mol_best_scores[mol_name] = record

            results = list(mol_best_scores.values())

        except Exception as e:
            logger.warning("Error parsing SDF file %s: %s", sdf_path, e)

        return results

    def _get_docking_stats(self) -> dict[str, Any]:
        """Get docking statistics."""
        docking_dir = self.base_path / STAGE_DIRS["docking"]

        # Detect available tools dynamically but always include common keys
        # expected by the report template.
        detected = self._detect_docking_tools()
        tools = sorted(set(detected) | {"gnina", "smina"})

        stats = {tool: {} for tool in tools}
        model_lookup = self._get_model_name_lookup()

        for tool in tools:
            # Try CSV first (legacy format)
            csv_paths = [
                docking_dir / tool / "scores.csv",
                docking_dir / f"{tool}_scores.csv",
                docking_dir / tool / "docking_results.csv",
            ]

            csv_found = False
            for csv_path in csv_paths:
                if csv_path.exists():
                    try:
                        df = pd.read_csv(csv_path)
                        for score_col in ["affinity", "score", "minimizedAffinity"]:
                            if score_col in df.columns:
                                stats[tool] = {
                                    "scores": df[score_col].dropna().tolist(),
                                    "count": len(df),
                                }
                                csv_found = True
                                break
                    except Exception:
                        pass
                    if csv_found:
                        break

            # Try SDF if CSV not found
            if not csv_found:
                sdf_paths = [
                    docking_dir / tool / f"{tool}_out.sdf",
                    docking_dir / tool / "output.sdf",
                    docking_dir / f"{tool}_out.sdf",
                ]

                for sdf_path in sdf_paths:
                    if sdf_path.exists():
                        records = self._parse_docking_sdf(sdf_path, model_lookup)
                        if records:
                            scores = [r["affinity"] for r in records]
                            stats[tool] = {
                                "scores": scores,
                                "count": len(records),
                            }
                            break

        return stats

    # =========================================================================
    # DETAILED DATA COLLECTION METHODS
    # =========================================================================

    def _get_descriptors_detailed(
        self, stage_key: str = "descriptors_initial"
    ) -> dict[str, Any]:
        """Get detailed descriptor data for enhanced visualization.

        Reads ALL numeric columns from the CSV file, not just predefined ones.

        Args:
            stage_key: Stage key from STAGE_DIRS (default: descriptors_initial)

        Returns:
            Dictionary with raw_data, summary_by_model, and key_descriptors
        """
        desc_dir = self.base_path / STAGE_DIRS[stage_key]

        paths_to_try = [
            desc_dir / "metrics" / "descriptors_all.csv",
            desc_dir / "filtered" / "descriptors_passed.csv",
            desc_dir / "filtered_molecules.csv",
        ]

        df = None
        for path in paths_to_try:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    break
                except Exception:
                    continue

        if df is None:
            return {}

        result = {
            "raw_data": [],
            "summary_by_model": {},
            "key_descriptors": ["MolWt", "LogP", "TPSA", "QED"],
        }

        # Build column name normalization map (actual_col -> display_name)
        # Use DESCRIPTOR_ALIASES for known mappings, otherwise keep original name
        col_normalize = {}
        for col in df.columns:
            col_lower = col.lower()
            if col_lower in DESCRIPTOR_EXCLUDE_COLS:
                continue
            if col_lower in DESCRIPTOR_ALIASES:
                col_normalize[col] = DESCRIPTOR_ALIASES[col_lower]
            else:
                # Keep original name for unknown columns
                col_normalize[col] = col

        # Identify numeric columns
        numeric_cols = []
        for col in col_normalize.keys():
            if df[col].dtype in ["int64", "float64", "int32", "float32"]:
                numeric_cols.append(col)

        # Collect raw data records with normalized column names
        for _, row in df.iterrows():
            record = {}
            # Add model_name
            if "model_name" in df.columns and pd.notna(row.get("model_name")):
                record["model_name"] = row.get("model_name")

            # Add ALL numeric descriptors with normalized names
            for actual_col in numeric_cols:
                display_name = col_normalize[actual_col]
                if pd.notna(row.get(actual_col)):
                    val = row.get(actual_col)
                    # Convert to Python native types
                    if isinstance(val, (int, float)):
                        record[display_name] = float(val)

            if record:
                result["raw_data"].append(record)

        # Calculate summary by model
        if "model_name" in df.columns:
            for model in df["model_name"].dropna().unique():
                model_df = df[df["model_name"] == model]
                model_summary = {}
                for actual_col in numeric_cols:
                    display_name = col_normalize[actual_col]
                    values = model_df[actual_col].dropna()
                    if len(values) > 0:
                        model_summary[display_name] = float(values.mean())
                if model_summary:
                    result["summary_by_model"][model] = model_summary

        return result

    def _build_descriptors_js_data(
        self, desc_detailed: dict[str, Any], available_models: list[str]
    ) -> dict[str, Any]:
        """Build JSON data for JavaScript descriptors visualization.

        Args:
            desc_detailed: Detailed descriptor data from _get_descriptors_detailed
            available_models: List of model names

        Returns:
            Dictionary with structure for JavaScript plotting
        """
        raw_data = desc_detailed.get("raw_data", [])
        if not raw_data:
            return {}

        # Thresholds for Lipinski/drug-likeness visualization
        thresholds = {
            "MolWt": {"min": 100, "max": 500},
            "LogP": {"min": -2, "max": 5},
            "cLogP": {"min": -2, "max": 5},
            "TPSA": {"min": 20, "max": 140},
            "NumHDonors": {"min": 0, "max": 5},
            "NumHAcceptors": {"min": 0, "max": 10},
            "NumRotatableBonds": {"min": 0, "max": 10},
            "n_rot_bonds": {"min": 0, "max": 10},
            "QED": {"min": 0.3, "max": 1},
            "Fsp3": {"min": 0, "max": 1},
            "n_atoms": {"min": 10, "max": 70},
            "n_heavy_atoms": {"min": 10, "max": 50},
            "n_rings": {"min": 1, "max": 6},
            "n_aroma_rings": {"min": 0, "max": 4},
            "n_rigid_bonds": {"min": 0, "max": 30},
        }

        # Get all descriptor names from raw_data
        all_descriptors = set()
        for record in raw_data:
            for key in record.keys():
                if key != "model_name":
                    all_descriptors.add(key)
        descriptor_names = sorted(list(all_descriptors))

        # Build "all" data (aggregated across all models)
        all_data = {}
        for desc_name in descriptor_names:
            values = [
                r.get(desc_name) for r in raw_data if r.get(desc_name) is not None
            ]
            if values:
                all_data[desc_name] = {
                    "values": values,
                    "mean": statistics.mean(values),
                    "median": statistics.median(values),
                    "std": statistics.stdev(values) if len(values) > 1 else 0,
                    "min": min(values),
                    "max": max(values),
                }

        # Build per-model data
        by_model = {}
        for model in available_models:
            model_records = [r for r in raw_data if r.get("model_name") == model]
            if not model_records:
                continue
            model_data = {}
            for desc_name in descriptor_names:
                values = [
                    r.get(desc_name)
                    for r in model_records
                    if r.get(desc_name) is not None
                ]
                if values:
                    model_data[desc_name] = {
                        "values": values,
                        "mean": statistics.mean(values),
                        "median": statistics.median(values),
                        "std": statistics.stdev(values) if len(values) > 1 else 0,
                        "min": min(values),
                        "max": max(values),
                    }
            if model_data:
                by_model[model] = model_data

        # Build compare data
        compare_data = {
            "is_comparison": True,
            "models": available_models,
            "model_colors": plots.COMPARE_PALETTE[: len(available_models)],
            "data": by_model,
        }

        return {
            "all": all_data,
            "by_model": by_model,
            "__compare__": compare_data,
            "models": available_models,
            "descriptors": descriptor_names,
            "thresholds": thresholds,
        }

    def _build_filters_js_data(self, available_models: list[str]) -> dict[str, Any]:
        """Build JSON data for JavaScript filters visualization.

        Args:
            available_models: List of model names

        Returns:
            Dictionary with structure for JavaScript plotting.
            Sub-metrics are converted from ratios to absolute molecule counts.
        """
        result = {
            "filters": [],
            "models": available_models,
            "filter_data": {},
        }

        # Read filter data from each filter subdirectory
        for stage_key in ["struct_filters_pre", "struct_filters_post"]:
            stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
            if not stage_dir.exists():
                continue

            for subdir in sorted(stage_dir.iterdir()):
                if not subdir.is_dir() or subdir.name in ["plots", "logs"]:
                    continue

                filter_name = subdir.name
                if filter_name in result["filter_data"]:
                    continue  # Already processed

                metrics_path = subdir / "metrics.csv"
                if not metrics_path.exists():
                    continue

                try:
                    df = pd.read_csv(metrics_path)
                except Exception:
                    continue

                if "model_name" not in df.columns:
                    continue

                filter_info = {
                    "models": available_models,
                    "num_mol": {},
                    "pass_rate": {},  # Pass rate percentage per model
                    "sub_metrics": {},
                    "sub_metric_names": [],
                }

                # Get number of molecules per model
                if "num_mol" in df.columns:
                    for model in available_models:
                        model_df = df[df["model_name"] == model]
                        if not model_df.empty:
                            filter_info["num_mol"][model] = int(
                                model_df["num_mol"].iloc[0]
                            )

                # Find ratio columns and convert to absolute counts
                # Look for columns with banned_ratio pattern
                exclude_cols = {"model_name", "num_mol", "smiles", "mol_idx"}
                all_cols = [c for c in df.columns if c not in exclude_cols]

                # Identify ratio columns (those that represent molecule fractions)
                ratio_cols = [
                    c
                    for c in all_cols
                    if c.endswith("_banned_ratio")
                    or c == "banned_ratio"
                    or c.startswith("banned_ratio_")  # e.g. banned_ratio_s_1
                    or c.endswith("_ratio")
                    or c == "all_banned_ratio"
                    or c == "any_banned_ratio"
                ]

                # If we have ratio columns, use them for sub-metrics
                if ratio_cols:
                    # Build clean names: remove banned_ratio patterns
                    sub_metric_names = []
                    col_to_name = {}
                    for col in ratio_cols:
                        if col == "banned_ratio":
                            name = "failed"
                        elif col == "all_banned_ratio":
                            name = "all"
                        elif col == "any_banned_ratio":
                            continue  # Skip any_banned_ratio (usually 0)
                        elif col.startswith("banned_ratio_"):
                            # e.g. banned_ratio_s_1 -> s_1
                            name = col.replace("banned_ratio_", "")
                        elif col.endswith("_banned_ratio"):
                            name = col.replace("_banned_ratio", "")
                        elif col.endswith("_ratio"):
                            name = col.replace("_ratio", "")
                        else:
                            name = col
                        sub_metric_names.append(name)
                        col_to_name[col] = name

                    filter_info["sub_metric_names"] = sub_metric_names

                    # Determine which column represents the main banned ratio
                    main_ratio_col = None
                    for candidate in ["all_banned_ratio", "banned_ratio"]:
                        if candidate in df.columns:
                            main_ratio_col = candidate
                            break

                    # Check if we need to add passed/failed pair
                    has_only_failed = sub_metric_names == ["failed"] or (
                        len(sub_metric_names) == 1 and "failed" in sub_metric_names
                    )
                    if has_only_failed:
                        filter_info["sub_metric_names"] = ["passed", "failed"]

                    for model in available_models:
                        model_df = df[df["model_name"] == model]
                        if not model_df.empty:
                            num_mol = int(model_df["num_mol"].iloc[0])
                            model_metrics = {}
                            for col, name in col_to_name.items():
                                if col in model_df.columns:
                                    ratio = model_df[col].iloc[0]
                                    if pd.notna(ratio):
                                        # Convert ratio to absolute count
                                        model_metrics[name] = int(
                                            round(num_mol * ratio)
                                        )
                                    else:
                                        model_metrics[name] = 0

                            # Add passed count if we only have failed
                            if has_only_failed and "failed" in model_metrics:
                                model_metrics["passed"] = (
                                    num_mol - model_metrics["failed"]
                                )

                            filter_info["sub_metrics"][model] = model_metrics

                            # Calculate pass_rate from main ratio column
                            if main_ratio_col and main_ratio_col in model_df.columns:
                                main_ratio = model_df[main_ratio_col].iloc[0]
                                if pd.notna(main_ratio):
                                    filter_info["pass_rate"][model] = round(
                                        (1 - main_ratio) * 100, 1
                                    )
                                else:
                                    filter_info["pass_rate"][model] = 100.0
                            else:
                                filter_info["pass_rate"][model] = 100.0
                else:
                    # No ratio columns - create simple passed/failed
                    if "banned_ratio" not in all_cols:
                        # Skip filters without meaningful data
                        continue
                    filter_info["sub_metric_names"] = ["passed", "failed"]
                    for model in available_models:
                        model_df = df[df["model_name"] == model]
                        if not model_df.empty:
                            num_mol = int(model_df["num_mol"].iloc[0])
                            ratio = model_df["banned_ratio"].iloc[0]
                            failed = (
                                int(round(num_mol * ratio)) if pd.notna(ratio) else 0
                            )
                            passed = num_mol - failed
                            filter_info["sub_metrics"][model] = {
                                "passed": passed,
                                "failed": failed,
                            }
                            # Calculate pass_rate
                            if pd.notna(ratio):
                                filter_info["pass_rate"][model] = round(
                                    (1 - ratio) * 100, 1
                                )
                            else:
                                filter_info["pass_rate"][model] = 100.0

                if filter_info["sub_metrics"]:
                    result["filters"].append(filter_name)
                    result["filter_data"][filter_name] = filter_info

        return result

    def _get_filters_detailed(self) -> dict[str, Any]:
        """Get detailed filter data for enhanced visualization.

        Returns:
            Dictionary with by_filter, banned_ratios, common_alerts_reasons
        """
        result = {
            "by_filter": {},
            "banned_ratios": {},
            "common_alerts_reasons": {},
            "filter_metrics": {},
        }

        for stage_key in ["struct_filters_pre", "struct_filters_post"]:
            stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
            if not stage_dir.exists():
                continue

            for subdir in stage_dir.iterdir():
                if not subdir.is_dir() or subdir.name in ["plots", "logs"]:
                    continue

                filter_name = subdir.name
                metrics_path = subdir / "metrics.csv"

                if metrics_path.exists():
                    try:
                        df = pd.read_csv(metrics_path)

                        # Get per-model breakdown
                        if "model_name" in df.columns:
                            if "failed_count" in df.columns:
                                model_counts = (
                                    df.groupby("model_name")["failed_count"]
                                    .sum()
                                    .to_dict()
                                )
                                result["by_filter"][filter_name] = {
                                    str(k): int(v) for k, v in model_counts.items()
                                }

                            # Calculate banned_ratio per model
                            if "banned_ratio" in df.columns:
                                ratios = (
                                    df.groupby("model_name")["banned_ratio"]
                                    .mean()
                                    .to_dict()
                                )
                                result["banned_ratios"][filter_name] = {
                                    str(k): float(v) for k, v in ratios.items()
                                }
                            elif (
                                "total_count" in df.columns
                                and "failed_count" in df.columns
                            ):
                                # Calculate from counts
                                model_ratios = {}
                                for model in df["model_name"].unique():
                                    model_df = df[df["model_name"] == model]
                                    total = model_df["total_count"].sum()
                                    failed = model_df["failed_count"].sum()
                                    if total > 0:
                                        model_ratios[str(model)] = failed / total
                                if model_ratios:
                                    result["banned_ratios"][filter_name] = model_ratios

                        # Collect filter-specific metrics
                        filter_metrics = {}

                        # Lilly filter metrics
                        if filter_name == "lilly" and "demerit_score" in df.columns:
                            filter_metrics["avg_demerit_score"] = float(
                                df["demerit_score"].mean()
                            )

                        # Molcomplexity metrics
                        if filter_name == "molcomplexity":
                            if "bertz" in df.columns:
                                filter_metrics["avg_bertz"] = float(df["bertz"].mean())
                            if "sas" in df.columns:
                                filter_metrics["avg_sas"] = float(df["sas"].mean())

                        if filter_metrics:
                            result["filter_metrics"][filter_name] = filter_metrics

                    except Exception:
                        continue

                # Check for common_alerts specific data
                if filter_name == "common_alerts":
                    alerts_path = subdir / "alerts_summary.csv"
                    if not alerts_path.exists():
                        alerts_path = subdir / "failed_molecules.csv"

                    if alerts_path.exists():
                        try:
                            alerts_df = pd.read_csv(alerts_path)
                            # Look for alert type columns
                            alert_cols = [
                                "alert_type",
                                "reason",
                                "alert_name",
                                "filter_reason",
                            ]
                            for col in alert_cols:
                                if col in alerts_df.columns:
                                    reasons = alerts_df[col].value_counts().to_dict()
                                    result["common_alerts_reasons"] = {
                                        str(k): int(v) for k, v in reasons.items()
                                    }
                                    break
                        except Exception:
                            pass

        return result

    def _get_synthesis_detailed(self) -> dict[str, Any]:
        """Get detailed synthesis data for enhanced visualization.

        Returns:
            Dictionary with raw_data, solved/unsolved counts, time stats, by_model
        """
        synth_dir = self.base_path / STAGE_DIRS["synthesis"]
        # Prefer synthesis_extended.csv which has solved/search_time columns
        scores_path = synth_dir / "synthesis_extended.csv"
        if not scores_path.exists():
            scores_path = synth_dir / "synthesis_scores.csv"

        if not scores_path.exists():
            return {}

        try:
            df = pd.read_csv(scores_path)
        except Exception:
            return {}

        result = {
            "raw_data": [],
            "sa_scores": [],
            "syba_scores": [],
            "ra_scores": [],
            "solved_count": 0,
            "unsolved_count": 0,
            "summary": {},
            "by_model": {},
        }

        # Extract SA, SYBA, and RA scores
        if "sa_score" in df.columns:
            result["sa_scores"] = df["sa_score"].dropna().tolist()
            if result["sa_scores"]:
                result["summary"]["avg_sa_score"] = float(df["sa_score"].mean())

        if "syba_score" in df.columns:
            result["syba_scores"] = df["syba_score"].dropna().tolist()
            if result["syba_scores"]:
                result["summary"]["avg_syba_score"] = float(df["syba_score"].mean())

        if "ra_score" in df.columns:
            result["ra_scores"] = df["ra_score"].dropna().tolist()
            if result["ra_scores"]:
                result["summary"]["avg_ra_score"] = float(
                    df["ra_score"].dropna().mean()
                )

        # Count solved/unsolved (look for solved, route_found, etc.)
        solved_col = None
        solved_cols = ["solved", "route_found", "synthesis_solved"]
        for col in solved_cols:
            if col in df.columns:
                solved_col = col
                solved = (
                    df[col].sum()
                    if df[col].dtype == bool
                    else df[col].astype(bool).sum()
                )
                result["solved_count"] = int(solved)
                result["unsolved_count"] = len(df) - int(solved)
                result["summary"]["pct_solved"] = (
                    100 * solved / len(df) if len(df) > 0 else 0
                )
                break

        # Collect time data if available
        time_col = None
        time_cols = ["search_time", "route_time", "synthesis_time"]
        for col in time_cols:
            if col in df.columns:
                times = df[col].dropna().tolist()
                if times:
                    time_col = col
                    result["summary"]["avg_search_time"] = float(df[col].mean())

                    # Collect raw data with model info for box plot
                    if "model_name" in df.columns:
                        for _, row in df.iterrows():
                            if pd.notna(row.get(col)):
                                result["raw_data"].append(
                                    {
                                        "model_name": row.get("model_name", "Unknown"),
                                        "search_time": row.get(col),
                                    }
                                )
                break

        # Group data by model for comparison views
        if "model_name" in df.columns:
            for model in df["model_name"].dropna().unique():
                model_df = df[df["model_name"] == model]
                model_data = {
                    "sa_scores": model_df["sa_score"].dropna().tolist()
                    if "sa_score" in df.columns
                    else [],
                    "syba_scores": model_df["syba_score"].dropna().tolist()
                    if "syba_score" in df.columns
                    else [],
                    "ra_scores": model_df["ra_score"].dropna().tolist()
                    if "ra_score" in df.columns
                    else [],
                    "solved_count": 0,
                    "unsolved_count": 0,
                    "summary": {},
                }

                # Calculate per-model summary
                if "sa_score" in model_df.columns and len(model_data["sa_scores"]) > 0:
                    model_data["summary"]["avg_sa_score"] = float(
                        model_df["sa_score"].mean()
                    )

                if (
                    "syba_score" in model_df.columns
                    and len(model_data["syba_scores"]) > 0
                ):
                    model_data["summary"]["avg_syba_score"] = float(
                        model_df["syba_score"].mean()
                    )

                if "ra_score" in model_df.columns and len(model_data["ra_scores"]) > 0:
                    model_data["summary"]["avg_ra_score"] = float(
                        model_df["ra_score"].dropna().mean()
                    )

                # Count solved/unsolved per model
                if solved_col and solved_col in model_df.columns:
                    model_solved = (
                        model_df[solved_col].sum()
                        if model_df[solved_col].dtype == bool
                        else (model_df[solved_col] == True).sum()  # noqa: E712
                    )
                    model_data["solved_count"] = int(model_solved)
                    model_data["unsolved_count"] = len(model_df) - int(model_solved)
                    model_data["summary"]["pct_solved"] = (
                        100 * model_solved / len(model_df) if len(model_df) > 0 else 0
                    )

                # Time data per model
                if time_col and time_col in model_df.columns:
                    model_times = model_df[time_col].dropna()
                    if len(model_times) > 0:
                        model_data["summary"]["avg_search_time"] = float(
                            model_times.mean()
                        )

                result["by_model"][model] = model_data

        return result

    def _get_retrosynthesis_detailed(self) -> dict[str, Any]:
        """Get detailed retrosynthesis data from AiZynthFinder results.

        Returns:
            Dictionary with route_scores, steps, precursors, solve_rate, summary
        """
        synth_dir = self.base_path / STAGE_DIRS["synthesis"]
        json_path = synth_dir / "retrosynthesis_results.json"

        if not json_path.exists():
            return {}

        try:
            with open(json_path) as f:
                data = json.load(f)
        except Exception:
            return {}

        if "data" not in data:
            return {}

        result = {
            "route_scores": [],
            "steps": [],
            "precursors": [],
            "solved_count": 0,
            "total_count": 0,
            "summary": {},
        }

        for item in data["data"]:
            result["total_count"] += 1
            if item.get("is_solved", False):
                result["solved_count"] += 1
                if item.get("top_score") is not None:
                    result["route_scores"].append(item["top_score"])
                if item.get("number_of_steps") is not None:
                    result["steps"].append(item["number_of_steps"])
                if item.get("number_of_precursors") is not None:
                    result["precursors"].append(item["number_of_precursors"])

        # Calculate summary statistics
        if result["total_count"] > 0:
            result["summary"]["solve_rate"] = (
                100.0 * result["solved_count"] / result["total_count"]
            )
        if result["route_scores"]:
            result["summary"]["avg_route_score"] = statistics.mean(
                result["route_scores"]
            )
        if result["steps"]:
            result["summary"]["avg_steps"] = statistics.mean(result["steps"])
        if result["precursors"]:
            result["summary"]["avg_precursors"] = statistics.mean(result["precursors"])

        return result

    def _get_existing_plots(self) -> dict[str, str]:
        """Collect existing plot images from stage directories.

        Returns:
            Dictionary mapping plot name -> base64-encoded image data URI
        """
        plots_found = {}

        # Define expected plot locations
        plot_locations = {
            # Descriptors
            "descriptors_initial_distribution": (
                STAGE_DIRS["descriptors_initial"]
                + "/plots/descriptors_distribution.png"
            ),
            "descriptors_final_distribution": (
                STAGE_DIRS["descriptors_final"] + "/plots/descriptors_distribution.png"
            ),
            # Structural filters
            "filters_pre_counts": (
                STAGE_DIRS["struct_filters_pre"]
                + "/plots/molecule_counts_comparison.png"
            ),
            "filters_pre_ratios": (
                STAGE_DIRS["struct_filters_pre"]
                + "/plots/restriction_ratios_comparison.png"
            ),
            "filters_post_counts": (
                STAGE_DIRS["struct_filters_post"]
                + "/plots/molecule_counts_comparison.png"
            ),
            "filters_post_ratios": (
                STAGE_DIRS["struct_filters_post"]
                + "/plots/restriction_ratios_comparison.png"
            ),
        }

        for name, rel_path in plot_locations.items():
            full_path = self.base_path / rel_path
            if full_path.exists():
                try:
                    with open(full_path, "rb") as f:
                        img_data = base64.b64encode(f.read()).decode("utf-8")
                    plots_found[name] = f"data:image/png;base64,{img_data}"
                    logger.debug("Found plot: %s", name)
                except Exception as e:
                    logger.warning("Failed to load plot %s: %s", name, e)

        return plots_found

    def _get_docking_detailed(self) -> dict[str, Any]:
        """Get detailed docking data for enhanced visualization.

        Supports both CSV and SDF file formats, with automatic tool detection.

        Returns:
            Dictionary with raw_data, top_molecules, summary stats, by_model
        """
        docking_dir = self.base_path / STAGE_DIRS["docking"]

        # Detect available tools dynamically but always include common keys
        # expected by the report template.
        detected = self._detect_docking_tools()
        tools = sorted(set(detected) | {"gnina", "smina"})

        # Initialize result with detected tools
        result = {}
        for tool in tools:
            result[tool] = {
                "raw_data": [],
                "top_molecules": [],
                "summary": {},
                "by_model": {},
            }

        model_lookup = self._get_model_name_lookup()

        for tool in tools:
            records = []

            # Try CSV files first
            csv_paths = [
                docking_dir / tool / "scores.csv",
                docking_dir / tool / "docking_results.csv",
                docking_dir / f"{tool}_scores.csv",
            ]

            df = None
            for path in csv_paths:
                if path.exists():
                    try:
                        df = pd.read_csv(path)
                        break
                    except Exception:
                        continue

            if df is not None:
                # Process CSV data
                score_col = None
                for col in ["affinity", "score", "minimizedAffinity", "docking_score"]:
                    if col in df.columns:
                        score_col = col
                        break

                if score_col:
                    id_col = None
                    for col in ["molecule_id", "mol_id", "name", "smiles"]:
                        if col in df.columns:
                            id_col = col
                            break

                    for _, row in df.iterrows():
                        if pd.notna(row.get(score_col)):
                            record = {
                                "molecule_id": row.get(id_col, "Unknown")
                                if id_col
                                else "Unknown",
                                "affinity": float(row[score_col]),
                                "model_name": row.get("model_name", "Unknown")
                                if "model_name" in df.columns
                                else model_lookup.get(
                                    str(row.get(id_col, "")), "Unknown"
                                ),
                            }
                            records.append(record)
            else:
                # Try SDF files
                sdf_paths = [
                    docking_dir / tool / f"{tool}_out.sdf",
                    docking_dir / tool / "output.sdf",
                    docking_dir / f"{tool}_out.sdf",
                ]

                for sdf_path in sdf_paths:
                    if sdf_path.exists():
                        records = self._parse_docking_sdf(sdf_path, model_lookup)
                        break

            if not records:
                continue

            # Store raw data
            result[tool]["raw_data"] = records

            # Calculate summary stats
            scores = [r["affinity"] for r in records]
            result[tool]["summary"] = {
                "avg_affinity": sum(scores) / len(scores),
                "best_affinity": min(scores),
                "count": len(records),
            }

            # Top 10 molecules by affinity (lowest is best)
            sorted_records = sorted(records, key=lambda x: x["affinity"])
            result[tool]["top_molecules"] = sorted_records[:10]

            # Group data by model
            models = set(
                r["model_name"] for r in records if r["model_name"] != "Unknown"
            )
            for model in models:
                model_records = [r for r in records if r["model_name"] == model]
                if not model_records:
                    continue

                model_scores = [r["affinity"] for r in model_records]
                model_sorted = sorted(model_records, key=lambda x: x["affinity"])

                result[tool]["by_model"][model] = {
                    "scores": model_scores,
                    "summary": {
                        "avg_affinity": sum(model_scores) / len(model_scores),
                        "best_affinity": min(model_scores),
                        "count": len(model_records),
                    },
                    "top_molecules": model_sorted[:10],
                }

        return result

    def _build_descriptor_comparison_data(
        self,
        initial_detailed: dict[str, Any],
        final_detailed: dict[str, Any],
    ) -> dict[str, Any]:
        """Build comparison data between initial and final descriptors.

        Args:
            initial_detailed: Detailed data from initial descriptors stage
            final_detailed: Detailed data from final descriptors stage

        Returns:
            Dictionary with per-descriptor initial/final value arrays and stats
        """
        initial_raw = initial_detailed.get("raw_data", [])
        final_raw = final_detailed.get("raw_data", [])
        if not initial_raw or not final_raw:
            return {}

        # Get all descriptor names present in both datasets
        initial_descs = set()
        for record in initial_raw:
            initial_descs.update(k for k in record if k != "model_name")

        final_descs = set()
        for record in final_raw:
            final_descs.update(k for k in record if k != "model_name")

        common_descs = sorted(initial_descs & final_descs)
        if not common_descs:
            return {}

        comparison = {}
        for desc in common_descs:
            init_vals = [r[desc] for r in initial_raw if r.get(desc) is not None]
            final_vals = [r[desc] for r in final_raw if r.get(desc) is not None]
            if init_vals and final_vals:
                comparison[desc] = {
                    "initial_values": init_vals,
                    "final_values": final_vals,
                    "mean_initial": statistics.mean(init_vals),
                    "mean_final": statistics.mean(final_vals),
                }

        return {"descriptors": list(comparison.keys()), "data": comparison}

    def _get_docking_filters_detailed(self) -> dict[str, Any]:
        """Get detailed docking filters data for report visualization.

        Reads metrics.csv and filtered_molecules.csv from the docking filters
        stage directory. Dynamically detects enabled filters from pass_* columns.

        Returns:
            Dictionary with total/passed/failed counts, per-filter stats,
            numeric metric distributions, and per-model breakdown.
        """
        df_dir = self.base_path / STAGE_DIRS["docking_filters"]
        metrics_path = df_dir / "metrics.csv"

        if not metrics_path.exists():
            return {}

        try:
            df = pd.read_csv(metrics_path)
        except Exception:
            return {}

        if df.empty:
            return {}

        # Detect pass_* columns (individual filters), excluding aggregate 'pass'
        pass_cols = [c for c in df.columns if c.startswith("pass_") and c != "pass"]

        total_poses = len(df)
        passed_poses = int(df["pass"].sum()) if "pass" in df.columns else 0
        pass_rate = (
            round(100.0 * passed_poses / total_poses, 1) if total_poses > 0 else 0.0
        )

        # Unique molecules that passed
        filtered_path = df_dir / "filtered_molecules.csv"
        unique_molecules_passed = 0
        if filtered_path.exists():
            try:
                fdf = pd.read_csv(filtered_path)
                unique_molecules_passed = len(fdf)
            except Exception:
                pass

        # Read aggregation mode from pipeline config
        aggregation_mode = "all"
        dock_filt_config = self.config.get("docking_filters", {})
        if isinstance(dock_filt_config, dict):
            agg = dock_filt_config.get("aggregation", {})
            if isinstance(agg, dict):
                aggregation_mode = agg.get("mode", "all")

        # Per-filter pass/fail stats
        per_filter = {}
        for col in pass_cols:
            filter_name = col.replace("pass_", "")
            if col in df.columns:
                total = int(df[col].notna().sum())
                passed = int(df[col].sum())
                per_filter[filter_name] = {
                    "passed": passed,
                    "total": total,
                    "pass_rate": round(100.0 * passed / total, 1) if total > 0 else 0.0,
                }

        # Numeric metric distributions for histograms
        metric_columns = [
            "clashes",
            "strain_energy",
            "min_conformer_rmsd",
            "shape_score",
            "n_hbonds",
            "frac_atoms_outside_box",
        ]
        numeric_metrics = {}
        for col in metric_columns:
            if col in df.columns:
                values = df[col].dropna().tolist()
                if values:
                    numeric_metrics[col] = values

        # Per-model breakdown
        by_model = {}
        if "model_name" in df.columns:
            for model in df["model_name"].dropna().unique():
                model_df = df[df["model_name"] == model]
                m_total = len(model_df)
                m_passed = (
                    int(model_df["pass"].sum()) if "pass" in model_df.columns else 0
                )
                by_model[str(model)] = {
                    "total": m_total,
                    "passed": m_passed,
                    "pass_rate": round(100.0 * m_passed / m_total, 1)
                    if m_total > 0
                    else 0.0,
                }

        # Extract thresholds from config for display
        thresholds = {}
        if isinstance(dock_filt_config, dict):
            pq = dock_filt_config.get("pose_quality", {})
            if isinstance(pq, dict):
                if pq.get("max_clashes") is not None:
                    thresholds["clashes"] = {"max": pq["max_clashes"]}
                if pq.get("max_strain_energy") is not None:
                    thresholds["strain_energy"] = {"max": pq["max_strain_energy"]}
            cd = dock_filt_config.get("conformer_deviation", {})
            if isinstance(cd, dict):
                if cd.get("max_rmsd_to_conformer") is not None:
                    thresholds["min_conformer_rmsd"] = {
                        "max": cd["max_rmsd_to_conformer"]
                    }
            sb = dock_filt_config.get("search_box", {})
            if isinstance(sb, dict):
                if sb.get("max_outside_fraction") is not None:
                    thresholds["frac_atoms_outside_box"] = {
                        "max": sb["max_outside_fraction"]
                    }
            ss = dock_filt_config.get("shepherd_score", {})
            if isinstance(ss, dict):
                if ss.get("min_shape_score") is not None:
                    thresholds["shape_score"] = {"min": ss["min_shape_score"]}

        return {
            "total_poses": total_poses,
            "passed_poses": passed_poses,
            "pass_rate": pass_rate,
            "unique_molecules_passed": unique_molecules_passed,
            "aggregation_mode": aggregation_mode,
            "per_filter": per_filter,
            "numeric_metrics": numeric_metrics,
            "thresholds": thresholds,
            "by_model": by_model,
        }

    def _get_config_summary(self) -> dict[str, Any]:
        """Get configuration summary."""
        # Get folder from config or use base_path
        folder = self.config.get("folder_to_save", "")
        if not folder:
            folder = str(self.base_path)

        # Get stages from self.stages or detect from directories
        if self.stages:
            stages_enabled = [s.name for s in self.stages if s.enabled]
        else:
            # Auto-detect from existing stage directories
            stages_enabled = []
            stages_dir = self.base_path / "stages"
            if stages_dir.exists():
                stage_names = {
                    "01_descriptors_initial": "Descriptors (Initial)",
                    "02_structural_filters_pre": "Structural Filters (Pre)",
                    "03_structural_filters_post": "Structural Filters (Post)",
                    "04_synthesis": "Synthesis Analysis",
                    "05_docking": "Molecular Docking",
                    "06_docking_filters": "Docking Filters",
                    "07_descriptors_final": "Descriptors (Final)",
                }
                for dir_name, display_name in stage_names.items():
                    if (stages_dir / dir_name).exists():
                        stages_enabled.append(display_name)

        return {
            "folder_to_save": folder,
            "stages_enabled": stages_enabled,
        }

    def _generate_plots(self, data: dict[str, Any]) -> dict[str, str]:
        """Generate all visualization plots.

        Args:
            data: Collected pipeline data

        Returns:
            Dictionary mapping plot names to HTML strings
        """
        plot_htmls = {}

        # Funnel chart (classic)
        plot_htmls["funnel"] = plots.plot_funnel(data.get("funnel", []))

        # Sankey diagram (improved flow visualization)
        plot_htmls["sankey"] = plots.plot_sankey(data.get("funnel", []))

        # Sankey data for each model (for JavaScript filtering)
        funnel_by_model = data.get("funnel_by_model", {})
        available_models = data.get("available_models", [])

        sankey_by_model = {"all": plots.plot_sankey_json(data.get("funnel", []))}
        for model, model_funnel in funnel_by_model.items():
            sankey_by_model[model] = plots.plot_sankey_json(model_funnel)

        # Add comparison data for "Compare All" option
        if available_models and funnel_by_model:
            sankey_by_model["__compare__"] = plots.plot_sankey_compare_json(
                funnel_by_model, available_models
            )

        plot_htmls["sankey_data"] = json.dumps(sankey_by_model)

        # Stage summary
        plot_htmls["stage_summary"] = plots.plot_stage_summary(data.get("stages", []))

        # Model comparison
        model_stats = data.get("models", [])
        if model_stats:
            plot_htmls["model_comparison"] = plots.plot_model_comparison(model_stats)
            plot_htmls["model_losses"] = plots.plot_model_stacked_losses(model_stats)

        # Descriptor distributions (original)
        desc_data = data.get("descriptors", {})
        if desc_data.get("distributions"):
            plot_htmls["descriptors"] = plots.plot_descriptor_distributions(
                desc_data["distributions"]
            )

        # Descriptors detailed (enhanced)
        desc_detailed = data.get("descriptors_detailed", {})
        if desc_detailed.get("raw_data"):
            plot_htmls["descriptors_violin"] = plots.plot_descriptors_violin_by_model(
                desc_detailed["raw_data"],
                desc_detailed.get("key_descriptors", ["MolWt", "LogP", "TPSA", "QED"]),
            )
            plot_htmls["descriptors_hbd_hba"] = plots.plot_descriptors_hbd_hba_box(
                desc_detailed["raw_data"]
            )
        if desc_detailed.get("summary_by_model"):
            plot_htmls["descriptors_table"] = plots.plot_descriptors_summary_table(
                desc_detailed["summary_by_model"]
            )

        # Filter analysis (original)
        filter_data = data.get("filters", {})
        if filter_data.get("by_filter"):
            plot_htmls["filter_heatmap"] = plots.plot_filter_heatmap(
                filter_data["by_filter"]
            )
        if filter_data.get("totals"):
            plot_htmls["filter_failures"] = plots.plot_top_filter_failures(
                filter_data["totals"]
            )

        # Filters detailed (enhanced)
        filters_detailed = data.get("filters_detailed", {})
        if filters_detailed.get("by_filter"):
            plot_htmls["filter_stacked_bar"] = plots.plot_filter_stacked_bar(
                filters_detailed["by_filter"]
            )
        if filters_detailed.get("banned_ratios"):
            plot_htmls["filter_banned_heatmap"] = (
                plots.plot_filter_banned_ratio_heatmap(
                    filters_detailed["banned_ratios"]
                )
            )
        if filters_detailed.get("common_alerts_reasons"):
            plot_htmls["filter_top_reasons"] = plots.plot_filter_top_reasons_bar(
                filters_detailed["common_alerts_reasons"]
            )

        # Synthesis scores (original)
        synth_data = data.get("synthesis", {})
        if synth_data.get("distributions"):
            plot_htmls["synthesis_dist"] = plots.plot_synthesis_distributions(
                synth_data["distributions"]
            )
        scatter = synth_data.get("scatter_data", {})
        if scatter.get("sa_scores") and scatter.get("syba_scores"):
            plot_htmls["synthesis_scatter"] = plots.plot_synthesis_scatter(
                scatter["sa_scores"],
                scatter["syba_scores"],
                scatter.get("model_names", []),
            )

        # Synthesis detailed (enhanced)
        synth_detailed = data.get("synthesis_detailed", {})
        if synth_detailed.get("sa_scores"):
            plot_htmls["synthesis_sa_hist"] = plots.plot_synthesis_sa_histogram(
                synth_detailed["sa_scores"]
            )
        if synth_detailed.get("syba_scores"):
            plot_htmls["synthesis_syba_hist"] = plots.plot_synthesis_syba_histogram(
                synth_detailed["syba_scores"]
            )
        if synth_detailed.get("ra_scores"):
            plot_htmls["synthesis_ra_hist"] = plots.plot_synthesis_ra_histogram(
                synth_detailed["ra_scores"]
            )
        if (
            synth_detailed.get("solved_count", 0) > 0
            or synth_detailed.get("unsolved_count", 0) > 0
        ):
            plot_htmls["synthesis_solved_pie"] = plots.plot_synthesis_solved_pie(
                synth_detailed.get("solved_count", 0),
                synth_detailed.get("unsolved_count", 0),
            )
        if synth_detailed.get("raw_data"):
            plot_htmls["synthesis_time_box"] = plots.plot_synthesis_time_box(
                synth_detailed["raw_data"]
            )

        # Retrosynthesis (AiZynthFinder) plots
        retrosynth = data.get("retrosynthesis", {})
        if retrosynth.get("route_scores"):
            plot_htmls["retrosynthesis_route_score_hist"] = (
                plots.plot_retrosynthesis_route_score_histogram(
                    retrosynth["route_scores"]
                )
            )
        if retrosynth.get("steps"):
            plot_htmls["retrosynthesis_steps_hist"] = (
                plots.plot_retrosynthesis_steps_histogram(retrosynth["steps"])
            )

        # Docking scores (original)
        docking_data = data.get("docking", {})
        if docking_data.get("gnina", {}).get("scores"):
            plot_htmls["docking_gnina"] = plots.plot_docking_distribution(
                docking_data["gnina"]["scores"], "gnina"
            )
        if docking_data.get("smina", {}).get("scores"):
            plot_htmls["docking_smina"] = plots.plot_docking_distribution(
                docking_data["smina"]["scores"], "smina"
            )

        # Docking detailed (enhanced)
        docking_detailed = data.get("docking_detailed", {})
        for tool in ["gnina", "smina"]:
            tool_data = docking_detailed.get(tool, {})
            if tool_data.get("raw_data"):
                scores = [
                    d.get("affinity")
                    for d in tool_data["raw_data"]
                    if d.get("affinity") is not None
                ]
                if scores:
                    plot_htmls[f"docking_{tool}_affinity_hist"] = (
                        plots.plot_docking_affinity_histogram(scores, tool)
                    )
                plot_htmls[f"docking_{tool}_affinity_box"] = (
                    plots.plot_docking_affinity_box(tool_data["raw_data"])
                )
            if tool_data.get("top_molecules"):
                plot_htmls[f"docking_{tool}_top_molecules"] = (
                    plots.plot_docking_top_molecules(tool_data["top_molecules"])
                )

        # Generate JSON data for JavaScript model filtering (Synthesis)
        synth_detailed = data.get("synthesis_detailed", {})
        if synth_detailed:
            synthesis_data = {
                "all": {
                    "sa_scores": synth_detailed.get("sa_scores", []),
                    "syba_scores": synth_detailed.get("syba_scores", []),
                    "ra_scores": synth_detailed.get("ra_scores", []),
                    "solved_count": synth_detailed.get("solved_count", 0),
                    "unsolved_count": synth_detailed.get("unsolved_count", 0),
                    "summary": synth_detailed.get("summary", {}),
                }
            }
            # Add per-model data
            for model, model_data in synth_detailed.get("by_model", {}).items():
                synthesis_data[model] = model_data

            # Add compare mode data
            by_model = synth_detailed.get("by_model", {})
            if by_model:
                synthesis_data["__compare__"] = {
                    "is_comparison": True,
                    "models": list(by_model.keys()),
                    "model_colors": plots.COMPARE_PALETTE[: len(by_model)],
                    "data": by_model,
                }
            plot_htmls["synthesis_data"] = json.dumps(synthesis_data)

        # Generate JSON data for JavaScript descriptors visualization
        desc_detailed = data.get("descriptors_detailed", {})
        available_models = data.get("available_models", [])
        if desc_detailed.get("raw_data"):
            descriptors_data = self._build_descriptors_js_data(
                desc_detailed, available_models
            )
            plot_htmls["descriptors_data"] = json.dumps(descriptors_data)

        # Generate JSON data for JavaScript filters visualization
        filters_detailed = data.get("filters_detailed", {})
        if available_models:
            filters_js_data = self._build_filters_js_data(available_models)
            if filters_js_data.get("filter_data"):
                plot_htmls["filters_data"] = json.dumps(filters_js_data)

        # Generate JSON data for JavaScript model filtering (Docking)
        for tool in ["gnina", "smina"]:
            tool_data = docking_detailed.get(tool, {})
            if tool_data.get("raw_data") or tool_data.get("by_model"):
                # Collect all scores for "all" view
                all_scores = [
                    d.get("affinity")
                    for d in tool_data.get("raw_data", [])
                    if d.get("affinity") is not None
                ]
                docking_tool_data = {
                    "all": {
                        "scores": all_scores,
                        "summary": tool_data.get("summary", {}),
                        "top_molecules": tool_data.get("top_molecules", []),
                    }
                }
                # Add per-model data
                for model, model_data in tool_data.get("by_model", {}).items():
                    docking_tool_data[model] = model_data

                # Add compare mode data
                by_model = tool_data.get("by_model", {})
                if by_model:
                    docking_tool_data["__compare__"] = {
                        "is_comparison": True,
                        "models": list(by_model.keys()),
                        "model_colors": plots.COMPARE_PALETTE[: len(by_model)],
                        "data": by_model,
                    }
                plot_htmls[f"docking_{tool}_data"] = json.dumps(docking_tool_data)

        # =====================================================================
        # Docking Filters (Stage 06)
        # =====================================================================
        df_detailed = data.get("docking_filters_detailed", {})
        if df_detailed.get("per_filter"):
            plot_htmls["docking_filters_pass_fail"] = (
                plots.plot_docking_filters_pass_fail_bar(df_detailed["per_filter"])
            )
        if df_detailed.get("numeric_metrics"):
            plot_htmls["docking_filters_metric_hists"] = (
                plots.plot_docking_filters_metric_histograms(
                    df_detailed["numeric_metrics"],
                    df_detailed.get("thresholds", {}),
                )
            )
        if df_detailed.get("by_model"):
            plot_htmls["docking_filters_by_model"] = (
                plots.plot_docking_filters_by_model_bar(df_detailed["by_model"])
            )

            # Build JSON data for JS model toggle
            df_by_model_js = {"all": df_detailed}
            for model, model_data in df_detailed["by_model"].items():
                df_by_model_js[model] = model_data
            plot_htmls["docking_filters_data"] = json.dumps(df_by_model_js)

        # =====================================================================
        # Final Descriptors (Stage 07)
        # =====================================================================
        final_desc_detailed = data.get("descriptors_final_detailed", {})
        if final_desc_detailed.get("raw_data"):
            final_desc_js = self._build_descriptors_js_data(
                final_desc_detailed, available_models
            )
            plot_htmls["final_descriptors_data"] = json.dumps(final_desc_js)

            if final_desc_detailed.get("summary_by_model"):
                plot_htmls["final_descriptors_table"] = (
                    plots.plot_descriptors_summary_table(
                        final_desc_detailed["summary_by_model"]
                    )
                )

        # Build initial vs final descriptor comparison data
        initial_desc = data.get("descriptors_detailed", {})
        final_desc = data.get("descriptors_final_detailed", {})
        if initial_desc.get("raw_data") and final_desc.get("raw_data"):
            comparison = self._build_descriptor_comparison_data(
                initial_desc, final_desc
            )
            if comparison:
                plot_htmls["descriptors_comparison_data"] = json.dumps(comparison)

        return plot_htmls

    def _render_template(self, data: dict[str, Any], plot_htmls: dict[str, str]) -> str:
        """Render the HTML template with data and plots.

        Args:
            data: Collected pipeline data
            plot_htmls: Dictionary of plot HTML strings

        Returns:
            Rendered HTML string
        """
        try:
            env = Environment(
                loader=PackageLoader("hedgehog.reporting", "templates"),
                autoescape=True,
            )
            template = env.get_template("report.html")
        except Exception as e:
            logger.warning("Could not load template, using inline template: %s", e)
            return self._render_inline_template(data, plot_htmls)

        return template.render(
            data=data,
            plots=plot_htmls,
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

    def _build_model_options(self, models: list[str]) -> str:
        """Build HTML option elements for model dropdown.

        Args:
            models: List of model names

        Returns:
            HTML string with option elements
        """
        return "".join(f'<option value="{m}">{m}</option>' for m in models)

    def _render_inline_template(
        self, data: dict[str, Any], plot_htmls: dict[str, str]
    ) -> str:
        """Render report using inline template (fallback).

        Args:
            data: Collected pipeline data
            plot_htmls: Dictionary of plot HTML strings

        Returns:
            Rendered HTML string
        """
        summary = data.get("summary", {})
        metadata = data.get("metadata", {})

        # Build stage status HTML
        stage_rows = ""
        for status in summary.get("stage_statuses", []):
            icon = "✓" if status["completed"] else ("✗" if status["enabled"] else "−")
            color = (
                "#2ecc71"
                if status["completed"]
                else ("#e74c3c" if status["enabled"] else "#95a5a6")
            )
            stage_rows += f"""
            <tr>
                <td>{status["name"]}</td>
                <td style="color: {color}; font-weight: bold;">{icon} {status["status"].upper()}</td>
            </tr>
            """

        # Build descriptor summary table
        desc_summary = data.get("descriptors", {}).get("summary", {})
        desc_rows = ""
        for desc, stats in desc_summary.items():
            desc_rows += f"""
            <tr>
                <td>{desc}</td>
                <td>{stats["mean"]:.2f}</td>
                <td>{stats["std"]:.2f}</td>
                <td>{stats["min"]:.2f}</td>
                <td>{stats["max"]:.2f}</td>
            </tr>
            """

        # Build model table
        model_rows = ""
        for model in data.get("models", []):
            retention = (
                100 * model["final"] / model["initial"] if model["initial"] > 0 else 0
            )
            model_rows += f"""
            <tr>
                <td>{model["model_name"]}</td>
                <td>{model["initial"]}</td>
                <td>{model["final"]}</td>
                <td>{retention:.1f}%</td>
            </tr>
            """

        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HEDGEHOG Pipeline Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --primary: #3498db;
            --success: #2ecc71;
            --danger: #e74c3c;
            --warning: #f39c12;
            --dark: #2c3e50;
            --light: #ecf0f1;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            color: var(--dark);
            background: #f5f6fa;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
        header {{
            background: linear-gradient(135deg, var(--primary), #2980b9);
            color: white;
            padding: 40px 20px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        header h1 {{ font-size: 2.5em; margin-bottom: 10px; }}
        header .meta {{ opacity: 0.9; font-size: 0.95em; }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .summary-card {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            text-align: center;
        }}
        .summary-card .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: var(--primary);
        }}
        .summary-card .label {{
            color: #7f8c8d;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        section {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        section h2 {{
            color: var(--dark);
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid var(--light);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--light);
        }}
        th {{
            background: var(--light);
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.85em;
            letter-spacing: 0.5px;
        }}
        tr:hover {{ background: #f8f9fa; }}
        .plot-container {{ margin: 20px 0; }}
        .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        @media (max-width: 768px) {{
            .two-col {{ grid-template-columns: 1fr; }}
            .summary-grid {{ grid-template-columns: repeat(2, 1fr); }}
        }}
        .retention-badge {{
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
        }}
        .retention-high {{ background: #d5f4e6; color: #27ae60; }}
        .retention-medium {{ background: #fef9e7; color: #f39c12; }}
        .retention-low {{ background: #fadbd8; color: #e74c3c; }}
        .model-filter {{
            display: flex;
            align-items: center;
            gap: 10px;
            margin-bottom: 15px;
        }}
        .model-filter label {{
            font-weight: 600;
            color: var(--dark);
        }}
        .model-filter select {{
            padding: 8px 12px;
            border: 1px solid var(--light);
            border-radius: 5px;
            font-size: 14px;
            background: white;
            cursor: pointer;
            min-width: 200px;
        }}
        .model-filter select:hover {{
            border-color: var(--primary);
        }}
        .view-toggle {{
            display: flex;
            gap: 10px;
            margin-bottom: 15px;
        }}
        .view-toggle button {{
            padding: 8px 16px;
            border: 1px solid var(--light);
            border-radius: 5px;
            background: white;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.2s;
        }}
        .view-toggle button.active {{
            background: var(--primary);
            color: white;
            border-color: var(--primary);
        }}
        .view-toggle button:hover:not(.active) {{
            border-color: var(--primary);
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>HEDGEHOG Pipeline Report</h1>
            <div class="meta">
                <p>Generated: {metadata.get("generated_at", "N/A")}</p>
                <p>Run path: {metadata.get("run_path", "N/A")}</p>
            </div>
        </header>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="value">{summary.get("initial_molecules", 0)}</div>
                <div class="label">Initial Molecules</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get("final_molecules", 0)}</div>
                <div class="label">Final Molecules</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get("retention_percent", "N/A")}</div>
                <div class="label">Retention Rate</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get("stages_completed", 0)}/{summary.get("stages_enabled", 0)}</div>
                <div class="label">Stages Completed</div>
            </div>
        </div>

        <section>
            <h2>Pipeline Flow</h2>
            <div class="model-filter">
                <label for="model-select">Filter by Model:</label>
                <select id="model-select" onchange="updateSankeyDiagram()">
                    <option value="all">All Models</option>
                    {self._build_model_options(data.get("available_models", []))}
                </select>
            </div>
            <div class="view-toggle">
                <button id="btn-sankey" class="active" onclick="showView('sankey')">Sankey Diagram</button>
                <button id="btn-funnel" onclick="showView('funnel')">Funnel Chart</button>
            </div>
            <div id="sankey-container" class="plot-container">{plot_htmls.get("sankey", "")}</div>
            <div id="funnel-container" class="plot-container" style="display: none;">{plot_htmls.get("funnel", "")}</div>
            <script>
                var sankeyDataByModel = {plot_htmls.get("sankey_data", "{{}}")};

                function showView(view) {{
                    document.getElementById('sankey-container').style.display = view === 'sankey' ? 'block' : 'none';
                    document.getElementById('funnel-container').style.display = view === 'funnel' ? 'block' : 'none';
                    document.getElementById('btn-sankey').className = view === 'sankey' ? 'active' : '';
                    document.getElementById('btn-funnel').className = view === 'funnel' ? 'active' : '';
                }}

                function updateSankeyDiagram() {{
                    var model = document.getElementById('model-select').value;
                    var data = sankeyDataByModel[model];
                    var container = document.getElementById('sankey-container');
                    if (!data || !data.labels || data.labels.length < 2) {{
                        container.textContent = 'No data available for this model';
                        container.style.textAlign = 'center';
                        container.style.color = 'gray';
                        container.style.padding = '40px';
                        return;
                    }}

                    var trace = {{
                        type: 'sankey',
                        arrangement: 'snap',
                        node: {{
                            pad: 20,
                            thickness: 20,
                            line: {{ color: 'black', width: 0.5 }},
                            label: data.labels,
                            color: data.node_colors,
                            x: data.x_positions,
                            y: data.y_positions
                        }},
                        link: {{
                            source: data.sources,
                            target: data.targets,
                            value: data.values,
                            color: data.link_colors
                        }}
                    }};

                    var layout = {{
                        font: {{ size: 12 }},
                        height: 450,
                        margin: {{ l: 20, r: 20, t: 20, b: 20 }}
                    }};

                    Plotly.react('sankey-container', [trace], layout);
                }}
            </script>
        </section>

        <section>
            <h2>Stage Summary</h2>
            <div class="two-col">
                <div>
                    <table>
                        <thead>
                            <tr><th>Stage</th><th>Status</th></tr>
                        </thead>
                        <tbody>{stage_rows}</tbody>
                    </table>
                </div>
                <div class="plot-container">{plot_htmls.get("stage_summary", "")}</div>
            </div>
        </section>

        {"<section><h2>Model Comparison</h2><div class='plot-container'>" + plot_htmls.get("model_comparison", "") + "</div><table><thead><tr><th>Model</th><th>Initial</th><th>Final</th><th>Retention</th></tr></thead><tbody>" + model_rows + "</tbody></table><div class='plot-container'>" + plot_htmls.get("model_losses", "") + "</div></section>" if model_rows else ""}

        {"<section><h2>Descriptor Analysis</h2><div class='plot-container'>" + plot_htmls.get("descriptors", "") + "</div><table><thead><tr><th>Descriptor</th><th>Mean</th><th>Std</th><th>Min</th><th>Max</th></tr></thead><tbody>" + desc_rows + "</tbody></table></section>" if desc_rows else ""}

        {"<section><h2>Filter Analysis</h2><div class='two-col'><div class='plot-container'>" + plot_htmls.get("filter_heatmap", "") + "</div><div class='plot-container'>" + plot_htmls.get("filter_failures", "") + "</div></div></section>" if plot_htmls.get("filter_heatmap") or plot_htmls.get("filter_failures") else ""}

        {"<section><h2>Synthesis Scores</h2><div class='plot-container'>" + plot_htmls.get("synthesis_dist", "") + "</div><div class='plot-container'>" + plot_htmls.get("synthesis_scatter", "") + "</div></section>" if plot_htmls.get("synthesis_dist") else ""}

        {"<section><h2>Docking Results</h2><div class='two-col'><div class='plot-container'>" + plot_htmls.get("docking_gnina", "") + "</div><div class='plot-container'>" + plot_htmls.get("docking_smina", "") + "</div></div></section>" if plot_htmls.get("docking_gnina") or plot_htmls.get("docking_smina") else ""}

        <footer style="text-align: center; padding: 20px; color: #7f8c8d; font-size: 0.9em;">
            <p>HEDGEHOG - Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration</p>
        </footer>
    </div>
</body>
</html>
        """
        return html

    def _save_report(self, html_content: str) -> Path:
        """Save the HTML report to file.

        Args:
            html_content: Rendered HTML content

        Returns:
            Path to the saved report
        """
        report_path = self.output_dir / "report.html"
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        return report_path

    def _save_json_data(self, data: dict[str, Any]) -> Path:
        """Save report data as JSON for programmatic access.

        Args:
            data: Collected pipeline data

        Returns:
            Path to the saved JSON file
        """
        json_path = self.output_dir / "report_data.json"

        # Clean data for JSON serialization
        clean_data = self._make_json_serializable(data)

        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(clean_data, f, indent=2)

        return json_path

    def _make_json_serializable(self, obj: Any) -> Any:
        """Convert object to JSON-serializable format.

        Args:
            obj: Object to convert

        Returns:
            JSON-serializable version of the object
        """
        if isinstance(obj, dict):
            return {k: self._make_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_json_serializable(item) for item in obj]
        elif isinstance(obj, (int, float, str, bool, type(None))):
            return obj
        elif hasattr(obj, "tolist"):  # numpy arrays
            return obj.tolist()
        else:
            return str(obj)
