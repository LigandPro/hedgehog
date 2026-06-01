"""Generate a stage-audit Jupyter notebook for HEDGEHOG report outputs."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from typing import Any

NOTEBOOK_FILENAME = "stage_filter_audit.ipynb"


def _source_lines(text: str) -> list[str]:
    normalized = textwrap.dedent(text).strip("\n")
    if not normalized:
        return []
    return [f"{line}\n" for line in normalized.splitlines()]


def _markdown_cell(text: str) -> dict[str, Any]:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": _source_lines(text),
    }


def _code_cell(code: str) -> dict[str, Any]:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": _source_lines(code),
    }


def build_stage_audit_notebook() -> dict[str, Any]:
    """Build the stage-audit notebook structure."""
    cells = [
        _markdown_cell(
            """
            # Experiment: HEDGEHOG Stage Audit Notebook

            Use this notebook to inspect which molecules passed, failed, or dropped at each
            pipeline stage and to compare them against the configured filter boundaries.

            Recommended workflow:

            1. Run the notebook top-to-bottom once.
            2. Change `STAGE_KEY`, `VIEW`, and `MODEL_NAME` in the controls cell.
            3. Re-run the display cells to inspect the relevant molecules with `mols2grid`.
            4. Use the threshold and boundary helpers to justify filter decisions in the final report.
            """
        ),
        _markdown_cell(
            """
            ## What this notebook reads

            The notebook works directly from the run directory and uses:

            - `report_data.json` for stage metadata
            - `configs/config_*.yml` for descriptor, synthesis, and docking-filter thresholds
            - `filtered_molecules.csv`, `failed_molecules.csv`, and detailed stage CSVs when available

            The grid viewer uses the `mols2grid` package.
            """
        ),
        _code_cell(
            """
            from __future__ import annotations

            import json
            from pathlib import Path

            import mols2grid
            import pandas as pd
            import yaml
            from IPython.display import Markdown, display

            RUN_DIR = next(
                (
                    candidate
                    for candidate in (Path.cwd(), *Path.cwd().parents)
                    if (candidate / "report_data.json").exists() and (candidate / "stages").exists()
                ),
                None,
            )
            if RUN_DIR is None:
                raise FileNotFoundError(
                    "Could not locate a HEDGEHOG run directory from the current working directory."
                )

            REPORT_DATA = json.loads((RUN_DIR / "report_data.json").read_text())
            CONFIG_DIR = RUN_DIR / "configs"


            def read_yaml(path: Path) -> dict:
                if not path.exists():
                    return {}
                with path.open() as handle:
                    data = yaml.safe_load(handle) or {}
                return data if isinstance(data, dict) else {}


            DESCRIPTORS_CONFIG = read_yaml(CONFIG_DIR / "config_descriptors.yml")
            SYNTHESIS_CONFIG = read_yaml(CONFIG_DIR / "config_synthesis.yml")
            DOCKING_FILTERS_CONFIG = read_yaml(CONFIG_DIR / "config_docking_filters.yml")
            STRUCT_FILTERS_CONFIG = read_yaml(CONFIG_DIR / "config_structFilters.yml")

            display(Markdown(f"**Run directory:** `{RUN_DIR}`"))
            """
        ),
        _code_cell(
            """
            IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]

            STAGE_SPECS = {
                "input": {
                    "label": "Input",
                    "detail": ["input/sampled_molecules.csv"],
                    "filtered": ["input/sampled_molecules.csv"],
                    "previous": None,
                },
                "mol_prep": {
                    "label": "Mol Prep",
                    "detail": ["stages/01_mol_prep/filtered_molecules.csv"],
                    "filtered": ["stages/01_mol_prep/filtered_molecules.csv"],
                    "failed": ["stages/01_mol_prep/failed_molecules.csv"],
                    "previous": "input",
                },
                "descriptors_initial": {
                    "label": "Initial Descriptors",
                    "detail": [
                        "stages/02_descriptors_initial/metrics/descriptors_all.csv",
                        "stages/02_descriptors_initial/filtered/descriptors_passed.csv",
                        "stages/02_descriptors_initial/filtered_molecules.csv",
                    ],
                    "filtered": [
                        "stages/02_descriptors_initial/filtered/filtered_molecules.csv",
                        "stages/02_descriptors_initial/filtered_molecules.csv",
                    ],
                    "passed": [
                        "stages/02_descriptors_initial/filtered/descriptors_passed.csv",
                    ],
                    "failed": [
                        "stages/02_descriptors_initial/filtered/descriptors_failed.csv",
                        "stages/02_descriptors_initial/filtered/failed_molecules.csv",
                        "stages/02_descriptors_initial/failed_molecules.csv",
                    ],
                    "previous": "mol_prep",
                },
                "struct_filters_post": {
                    "label": "Structural Filters",
                    "detail": [
                        "stages/03_structural_filters_post/failed_molecules.csv",
                        "stages/03_structural_filters_post/filtered_molecules.csv",
                    ],
                    "filtered": ["stages/03_structural_filters_post/filtered_molecules.csv"],
                    "failed": ["stages/03_structural_filters_post/failed_molecules.csv"],
                    "previous": "descriptors_initial",
                },
                "synthesis": {
                    "label": "Synthesis",
                    "detail": [
                        "stages/04_synthesis/synthesis_extended.csv",
                        "stages/04_synthesis/synthesis_scores.csv",
                    ],
                    "filtered": ["stages/04_synthesis/filtered_molecules.csv"],
                    "previous": "struct_filters_post",
                },
                "docking_filters": {
                    "label": "Docking Filters",
                    "detail": ["stages/06_docking_filters/metrics.csv"],
                    "filtered": [
                        "stages/06_docking_filters/filtered_poses.csv",
                        "stages/06_docking_filters/filtered_molecules.csv",
                    ],
                    "failed": ["stages/06_docking_filters/failed_molecules.csv"],
                    "pass_col": "pass",
                    "previous": "synthesis",
                },
                "descriptors_final": {
                    "label": "Final Descriptors",
                    "detail": [
                        "stages/07_descriptors_final/metrics/descriptors_all.csv",
                        "stages/07_descriptors_final/filtered/descriptors_passed.csv",
                        "stages/07_descriptors_final/filtered_molecules.csv",
                    ],
                    "filtered": [
                        "stages/07_descriptors_final/filtered/filtered_molecules.csv",
                        "stages/07_descriptors_final/filtered_molecules.csv",
                    ],
                    "passed": [
                        "stages/07_descriptors_final/filtered/descriptors_passed.csv",
                    ],
                    "failed": [
                        "stages/07_descriptors_final/filtered/descriptors_failed.csv",
                        "stages/07_descriptors_final/filtered/failed_molecules.csv",
                        "stages/07_descriptors_final/failed_molecules.csv",
                    ],
                    "previous": "docking_filters",
                },
                "final_output": {
                    "label": "Final Output",
                    "detail": ["output/final_molecules.csv"],
                    "filtered": ["output/final_molecules.csv"],
                    "previous": "descriptors_final",
                },
            }


            def candidate_paths(stage_key: str, field: str) -> list[Path]:
                spec = STAGE_SPECS[stage_key]
                return [RUN_DIR / rel_path for rel_path in spec.get(field, [])]


            def read_first_csv(stage_key: str, field: str) -> pd.DataFrame | None:
                for path in candidate_paths(stage_key, field):
                    if not path.exists():
                        continue
                    try:
                        df = pd.read_csv(path)
                    except Exception:
                        continue
                    df = df.copy()
                    df.attrs["source_path"] = str(path.relative_to(RUN_DIR))
                    return df
                return None


            def available_stage_keys() -> list[str]:
                keys = []
                for stage_key in STAGE_SPECS:
                    has_data = any(path.exists() for path in candidate_paths(stage_key, "detail"))
                    has_data = has_data or any(
                        path.exists() for path in candidate_paths(stage_key, "filtered")
                    )
                    if has_data:
                        keys.append(stage_key)
                return keys


            def merge_keys(left: pd.DataFrame, right: pd.DataFrame) -> list[str]:
                candidates = [
                    ["mol_idx", "model_name"],
                    ["mol_idx"],
                    ["smiles", "model_name"],
                    ["smiles"],
                ]
                for candidate in candidates:
                    if all(col in left.columns for col in candidate) and all(
                        col in right.columns for col in candidate
                    ):
                        return candidate
                shared = [
                    col for col in IDENTITY_COLUMNS if col in left.columns and col in right.columns
                ]
                if shared:
                    return shared
                raise KeyError("No shared identity columns were found for the merge.")


            def dedupe_rows(df: pd.DataFrame, by_molecule: bool = True) -> pd.DataFrame:
                if df is None or df.empty or not by_molecule:
                    return df
                if {"mol_idx", "model_name"} <= set(df.columns):
                    return df.drop_duplicates(subset=["mol_idx", "model_name"]).copy()
                if "mol_idx" in df.columns:
                    return df.drop_duplicates(subset=["mol_idx"]).copy()
                if {"smiles", "model_name"} <= set(df.columns):
                    return df.drop_duplicates(subset=["smiles", "model_name"]).copy()
                if "smiles" in df.columns:
                    return df.drop_duplicates(subset=["smiles"]).copy()
                return df.copy()


            def stage_passed_frame(stage_key: str) -> pd.DataFrame:
                explicit_passed = read_first_csv(stage_key, "passed")
                if explicit_passed is not None:
                    return explicit_passed.copy()

                detail_df = read_first_csv(stage_key, "detail")
                filtered_df = read_first_csv(stage_key, "filtered")
                pass_col = STAGE_SPECS[stage_key].get("pass_col")

                if detail_df is not None and pass_col and pass_col in detail_df.columns:
                    return detail_df[detail_df[pass_col].fillna(False)].copy()

                if detail_df is not None and filtered_df is not None:
                    keys = merge_keys(detail_df, filtered_df)
                    return detail_df.merge(
                        filtered_df[keys].drop_duplicates(),
                        on=keys,
                        how="inner",
                    )

                if filtered_df is not None:
                    return filtered_df.copy()
                if detail_df is not None:
                    return detail_df.copy()
                return pd.DataFrame(columns=IDENTITY_COLUMNS)


            def stage_failed_frame(stage_key: str) -> pd.DataFrame:
                explicit_failed = read_first_csv(stage_key, "failed")
                if explicit_failed is not None:
                    return explicit_failed.copy()

                detail_df = read_first_csv(stage_key, "detail")
                filtered_df = read_first_csv(stage_key, "filtered")
                pass_col = STAGE_SPECS[stage_key].get("pass_col")

                if detail_df is not None and pass_col and pass_col in detail_df.columns:
                    return detail_df[~detail_df[pass_col].fillna(False)].copy()

                previous_key = STAGE_SPECS[stage_key].get("previous")
                if previous_key is None:
                    return pd.DataFrame(columns=IDENTITY_COLUMNS)

                previous_df = stage_passed_frame(previous_key)
                current_df = filtered_df if filtered_df is not None else stage_passed_frame(stage_key)
                if previous_df.empty:
                    return pd.DataFrame(columns=IDENTITY_COLUMNS)
                if current_df.empty:
                    return previous_df.copy()

                keys = merge_keys(previous_df, current_df)
                merged = previous_df.merge(
                    current_df[keys].drop_duplicates(),
                    on=keys,
                    how="left",
                    indicator=True,
                )
                return merged.loc[merged["_merge"] == "left_only"].drop(columns=["_merge"])


            def stage_view(
                stage_key: str,
                view: str = "dropped",
                model_name: str | None = None,
                by_molecule: bool = True,
            ) -> pd.DataFrame:
                if view == "passed":
                    df = stage_passed_frame(stage_key)
                elif view == "failed":
                    df = stage_failed_frame(stage_key)
                elif view == "all":
                    detail_df = read_first_csv(stage_key, "detail")
                    df = detail_df if detail_df is not None else stage_passed_frame(stage_key)
                else:
                    df = stage_failed_frame(stage_key)

                if df is None or df.empty:
                    return pd.DataFrame(columns=IDENTITY_COLUMNS)

                out = df.copy()
                if model_name and "model_name" in out.columns:
                    out = out[out["model_name"] == model_name].copy()
                out = dedupe_rows(out, by_molecule=by_molecule)

                preferred = [col for col in IDENTITY_COLUMNS if col in out.columns]
                trailing = [col for col in out.columns if col not in preferred]
                return out[preferred + trailing]


            def descriptor_threshold_table() -> pd.DataFrame:
                borders = DESCRIPTORS_CONFIG.get("borders", {})
                rows = []
                grouped: dict[str, dict[str, object]] = {}
                for key, value in borders.items():
                    if key.endswith("_min") or key.endswith("_max"):
                        metric, bound = key.rsplit("_", 1)
                        grouped.setdefault(metric, {})[bound] = value
                for metric, bounds in sorted(grouped.items()):
                    for bound, value in bounds.items():
                        rows.append(
                            {
                                "metric": metric,
                                "bound": bound,
                                "value": value,
                                "config_key": f"borders.{metric}_{bound}",
                            }
                        )
                if "allowed_chars" in borders:
                    rows.append(
                        {
                            "metric": "allowed_chars",
                            "bound": "allowed",
                            "value": ",".join(map(str, borders["allowed_chars"])),
                            "config_key": "borders.allowed_chars",
                        }
                    )
                if "charged_mol_allowed" in borders:
                    rows.append(
                        {
                            "metric": "charged_mol_allowed",
                            "bound": "allowed",
                            "value": borders["charged_mol_allowed"],
                            "config_key": "borders.charged_mol_allowed",
                        }
                    )
                return pd.DataFrame(rows)


            def synthesis_threshold_table() -> pd.DataFrame:
                rows = []
                for key, value in sorted(SYNTHESIS_CONFIG.items()):
                    if key.endswith("_min") or key.endswith("_max"):
                        metric, bound = key.rsplit("_", 1)
                        rows.append(
                            {
                                "metric": metric,
                                "bound": bound,
                                "value": value,
                                "config_key": key,
                            }
                        )
                return pd.DataFrame(rows)


            def docking_filter_threshold_table() -> pd.DataFrame:
                rows = []

                search_box = DOCKING_FILTERS_CONFIG.get("search_box", {})
                if "max_outside_fraction" in search_box:
                    rows.append(
                        {
                            "metric": "frac_atoms_outside_box",
                            "bound": "max",
                            "value": search_box["max_outside_fraction"],
                            "config_key": "search_box.max_outside_fraction",
                        }
                    )

                pose_quality = DOCKING_FILTERS_CONFIG.get("pose_quality", {})
                if "max_clashes" in pose_quality:
                    rows.append(
                        {
                            "metric": "clashes",
                            "bound": "max",
                            "value": pose_quality["max_clashes"],
                            "config_key": "pose_quality.max_clashes",
                        }
                    )
                if "max_strain_energy" in pose_quality:
                    rows.append(
                        {
                            "metric": "strain_energy",
                            "bound": "max",
                            "value": pose_quality["max_strain_energy"],
                            "config_key": "pose_quality.max_strain_energy",
                        }
                    )

                interactions = DOCKING_FILTERS_CONFIG.get("interactions", {})
                if "min_hbonds" in interactions:
                    rows.append(
                        {
                            "metric": "hbonds",
                            "bound": "min",
                            "value": interactions["min_hbonds"],
                            "config_key": "interactions.min_hbonds",
                        }
                    )

                shepherd_score = DOCKING_FILTERS_CONFIG.get("shepherd_score", {})
                if "min_shape_score" in shepherd_score:
                    rows.append(
                        {
                            "metric": "shape_score",
                            "bound": "min",
                            "value": shepherd_score["min_shape_score"],
                            "config_key": "shepherd_score.min_shape_score",
                        }
                    )

                conformer = DOCKING_FILTERS_CONFIG.get("conformer_deviation", {})
                if "max_rmsd_to_conformer" in conformer:
                    rows.append(
                        {
                            "metric": "min_conformer_rmsd",
                            "bound": "max",
                            "value": conformer["max_rmsd_to_conformer"],
                            "config_key": "conformer_deviation.max_rmsd_to_conformer",
                        }
                    )
                return pd.DataFrame(rows)


            def stage_thresholds(stage_key: str) -> pd.DataFrame:
                if stage_key in {"descriptors_initial", "descriptors_final"}:
                    return descriptor_threshold_table()
                if stage_key == "synthesis":
                    return synthesis_threshold_table()
                if stage_key == "docking_filters":
                    return docking_filter_threshold_table()
                return pd.DataFrame(columns=["metric", "bound", "value", "config_key"])


            def stage_overview() -> pd.DataFrame:
                rows = []
                for stage_key in available_stage_keys():
                    passed_df = stage_passed_frame(stage_key)
                    dropped_df = stage_failed_frame(stage_key)
                    rows.append(
                        {
                            "stage_key": stage_key,
                            "stage_label": STAGE_SPECS[stage_key]["label"],
                            "passed_rows": len(passed_df),
                            "dropped_rows": len(dropped_df),
                            "threshold_rows": len(stage_thresholds(stage_key)),
                        }
                    )
                return pd.DataFrame(rows)


            def default_subset(df: pd.DataFrame) -> list[str]:
                lead = [col for col in ["model_name", "mol_idx", "reason", "step"] if col in df.columns]
                numeric = [
                    col
                    for col in [
                        "MolWt",
                        "LogP",
                        "TPSA",
                        "QED",
                        "sa_score",
                        "syba_score",
                        "ra_score",
                        "route_score",
                        "num_steps",
                        "gnina_minimizedAffinity",
                        "clashes",
                        "strain_energy",
                        "min_conformer_rmsd",
                        "frac_atoms_outside_box",
                        "shape_score",
                    ]
                    if col in df.columns
                ]
                pass_flags = [col for col in df.columns if col.startswith("pass_")][:6]
                subset = lead + numeric + pass_flags
                return subset if subset else [col for col in df.columns if col != "smiles"][:8]


            def render_stage_grid(
                stage_key: str,
                view: str = "dropped",
                model_name: str | None = None,
                limit: int = 100,
                by_molecule: bool = True,
                subset: list[str] | None = None,
                tooltip: list[str] | None = None,
            ) -> pd.DataFrame:
                df = stage_view(stage_key, view=view, model_name=model_name, by_molecule=by_molecule)
                if df.empty:
                    display(
                        Markdown(
                            f"No rows available for stage `{stage_key}` and view `{view}`."
                        )
                    )
                    return df

                subset = subset or default_subset(df)
                tooltip = tooltip or subset
                if "smiles" not in df.columns:
                    raise KeyError("The selected dataset does not contain a 'smiles' column.")

                display_df = df.head(limit).copy()
                mols2grid.display(
                    display_df,
                    smiles_col="smiles",
                    subset=subset,
                    tooltip=tooltip,
                )
                return display_df


            def molecules_near_threshold(
                stage_key: str,
                metric: str,
                window: float,
                view: str = "dropped",
                model_name: str | None = None,
                by_molecule: bool = True,
            ) -> pd.DataFrame:
                df = stage_view(stage_key, view=view, model_name=model_name, by_molecule=by_molecule)
                thresholds = stage_thresholds(stage_key)
                if df.empty:
                    return df
                if metric not in df.columns:
                    raise KeyError(f"Metric '{metric}' is not present in the selected dataset.")

                metric_thresholds = thresholds[thresholds["metric"] == metric]
                if metric_thresholds.empty:
                    raise KeyError(f"No thresholds were found for metric '{metric}' in stage '{stage_key}'.")

                numeric_values = pd.to_numeric(df[metric], errors="coerce")
                mask = pd.Series(False, index=df.index)
                for _, row in metric_thresholds.iterrows():
                    try:
                        threshold_value = float(row["value"])
                    except (TypeError, ValueError):
                        continue
                    if row["bound"] == "min":
                        mask |= numeric_values.between(threshold_value, threshold_value + window)
                    elif row["bound"] == "max":
                        mask |= numeric_values.between(threshold_value - window, threshold_value)
                return df.loc[mask].copy()
            """
        ),
        _code_cell(
            """
            overview_df = stage_overview()
            overview_df
            """
        ),
        _markdown_cell(
            """
            ## Controls

            Change the values below, then re-run the next cells.
            """
        ),
        _code_cell(
            """
            AVAILABLE_STAGE_KEYS = available_stage_keys()
            AVAILABLE_STAGE_KEYS

            STAGE_KEY = "descriptors_initial"
            VIEW = "dropped"  # passed, dropped, failed, all
            MODEL_NAME = None  # e.g. "ModelA"
            LIMIT = 120
            DEDUPE_BY_MOLECULE = True
            """
        ),
        _code_cell(
            """
            display(Markdown(f"### {STAGE_SPECS[STAGE_KEY]['label']}"))
            threshold_df = stage_thresholds(STAGE_KEY)
            if threshold_df.empty:
                display(Markdown("No numeric threshold table is defined for this stage."))
            else:
                display(threshold_df)

            stage_df = stage_view(
                STAGE_KEY,
                view=VIEW,
                model_name=MODEL_NAME,
                by_molecule=DEDUPE_BY_MOLECULE,
            )
            stage_df.head(10)
            """
        ),
        _code_cell(
            """
            _ = render_stage_grid(
                STAGE_KEY,
                view=VIEW,
                model_name=MODEL_NAME,
                limit=LIMIT,
                by_molecule=DEDUPE_BY_MOLECULE,
            )
            """
        ),
        _markdown_cell(
            """
            ## Boundary inspection

            Use this section for stages with numeric thresholds (descriptors, synthesis, docking filters).
            Set the metric and an absolute window, then re-run the next cells to inspect compounds sitting
            close to the configured boundary.
            """
        ),
        _code_cell(
            """
            METRIC = "MolWt"
            WINDOW = 15.0

            near_boundary_df = molecules_near_threshold(
                STAGE_KEY,
                metric=METRIC,
                window=WINDOW,
                view=VIEW,
                model_name=MODEL_NAME,
                by_molecule=DEDUPE_BY_MOLECULE,
            )
            near_boundary_df.head(10)
            """
        ),
        _code_cell(
            """
            if near_boundary_df.empty:
                display(Markdown("No molecules were found in the requested boundary window."))
            else:
                mols2grid.display(
                    near_boundary_df.head(LIMIT),
                    smiles_col="smiles",
                    subset=default_subset(near_boundary_df),
                    tooltip=default_subset(near_boundary_df),
                )
            """
        ),
    ]

    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "version": "3.12",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def write_stage_audit_notebook(base_path: Path) -> Path:
    """Write the stage-audit notebook into the run directory."""
    out_path = Path(base_path) / NOTEBOOK_FILENAME
    notebook = build_stage_audit_notebook()
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(notebook, handle, indent=2)
        handle.write("\n")
    return out_path
