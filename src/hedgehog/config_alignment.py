"""Target-molecule percentile alignment for HEDGEHOG stage configs."""

from __future__ import annotations

import ast
import copy
import math
import shutil
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from rdkit import Chem

from hedgehog._constants import (
    KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS,
    TOOL_GNINA,
    TOOL_MATCHA,
    TOOL_SMINA,
)
from hedgehog.configs.logger import load_config, logger
from hedgehog.docking.metadata import _parse_tools_config

ALIGNMENT_DIR_NAME = "target_alignment"
ALIGNED_CONFIG_NAME = "aligned_config.yml"
THRESHOLDS_NAME = "alignment_thresholds.yml"
PROBE_CONFIGS_DIR_NAME = "calibration_configs_unfiltered"
TARGET_CALIBRATION_RUN_DIR_NAME = "calibration_target_run"

_CONFIG_MOL_PREP = "config_mol_prep"
_CONFIG_DESCRIPTORS = "config_descriptors"
_CONFIG_STRUCT_FILTERS = "config_structFilters"
_CONFIG_SYNTHESIS = "config_synthesis"
_CONFIG_DOCKING = "config_docking"
_CONFIG_DOCKING_FILTERS = "config_docking_filters"

_SYNTHESIS_LEGACY_FILTERS = {
    "sa_score": ("sa_score_min", "sa_score_max"),
    "syba_score": ("syba_score_min", "syba_score_max"),
    "ra_score": ("ra_score_min", "ra_score_max"),
}

_DEFAULT_SYNTHESIS_SCORERS = ("sa", "syba", "rascore")
_SYNTHESIS_SCORER_COLUMNS = {
    "sa": "sa_score",
    "syba": "syba_score",
    "rascore": "ra_score",
    "ra": "ra_score",
    "sync": "sync_score",
    "scscore": "sc_score",
    "sc": "sc_score",
    "nonpher": "nonpher_complexity_score",
    "fsscore": "fs_score",
    "fs": "fs_score",
    "gasa": "gasa_score",
}


def _enabled_synthesis_score_columns(config: dict[str, Any]) -> set[str]:
    """Return score columns calculated by the source synthesis config."""
    raw_enabled = config.get("enabled_scores")
    if raw_enabled is None:
        raw_enabled = _DEFAULT_SYNTHESIS_SCORERS
    if isinstance(raw_enabled, str):
        names = raw_enabled.replace(",", " ").split()
    else:
        names = raw_enabled
    return {
        column
        for raw_name in names
        if (column := _SYNTHESIS_SCORER_COLUMNS.get(str(raw_name).strip().lower()))
        is not None
    }


_DOCKING_SCORE_PROPERTIES = {
    TOOL_SMINA: ("minimizedAffinity", "affinity", "score"),
    TOOL_GNINA: ("minimizedAffinity", "affinity", "score"),
    TOOL_MATCHA: ("minimizedAffinity", "affinity", "score"),
}

_STRUCTURAL_ELEMENT_COLUMNS = {
    "N": "n_N_atoms",
    "O": "n_O_atoms",
    "S": "n_S_atoms",
}

_STRUCTURAL_DIRECT_COLUMNS = {
    "max_n_or_o_atoms": "n_NO_atoms",
    "max_small_rings_3_4": "n_small_rings_3_4",
    "max_acyclic_chain_length": "max_acyclic_chain_length",
}


def validate_target_coverage_percent(percentile: float) -> float:
    """Validate and normalize requested target coverage."""
    if not isinstance(percentile, (int, float)) or isinstance(percentile, bool):
        raise ValueError(
            "Target coverage percent must be a number greater than 0 and at most 100."
        )
    value = float(percentile)
    if not math.isfinite(value) or value <= 0.0 or value > 100.0:
        raise ValueError(
            "Target coverage percent must be greater than 0 and at most 100."
        )
    return value


# Backward-compatible Python API alias; config files should use target_coverage_percent.
validate_alignment_percentile = validate_target_coverage_percent


def _dump_yaml(data: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)


def _copy_master_configs(master: dict[str, Any], destination: Path) -> dict[str, Any]:
    """Copy all referenced config files and point a master copy at the copies."""
    copied = copy.deepcopy(master)
    destination.mkdir(parents=True, exist_ok=True)
    for key, raw_path in list(master.items()):
        if (
            not key.startswith("config_")
            or not isinstance(raw_path, str)
            or not raw_path
        ):
            continue
        source = Path(raw_path)
        if not source.exists() or not source.is_file():
            continue
        suffix = source.suffix or ".yml"
        target = destination / f"{key}{suffix}"
        shutil.copyfile(source, target)
        copied[key] = str(target.resolve())
    return copied


def _update_yaml(path: str | None, updater) -> None:
    if not path:
        return
    config_path = Path(path)
    if not config_path.exists():
        return
    config = load_config(str(config_path))
    updater(config)
    _dump_yaml(config, config_path)


def _target_atom_symbols(
    molecules: pd.DataFrame, smiles_column: str = "smiles"
) -> list[str]:
    """Return the periodic-table ordered union of atoms in target molecules."""
    column = smiles_column if smiles_column in molecules.columns else None
    if column is None:
        by_lower = {str(name).lower(): str(name) for name in molecules.columns}
        column = by_lower.get("smiles")
    if column is None:
        raise ValueError("Target molecules do not contain a SMILES column.")

    symbols: set[str] = set()
    invalid_count = 0
    for raw_smiles in molecules[column].dropna():
        smiles = str(raw_smiles).strip()
        if not smiles:
            invalid_count += 1
            continue
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            invalid_count += 1
            continue
        symbols.update(atom.GetSymbol() for atom in mol.GetAtoms())

    if not symbols:
        raise ValueError("No atom symbols could be read from the target molecules.")
    if invalid_count:
        logger.warning(
            "Ignored %d invalid target molecule(s) while deriving MolPrep allowed atoms.",
            invalid_count,
        )

    periodic_table = Chem.GetPeriodicTable()
    return sorted(symbols, key=lambda symbol: periodic_table.GetAtomicNumber(symbol))


def set_probe_molprep_allowed_atoms(
    probe: dict[str, Any], molecules: pd.DataFrame
) -> list[str] | None:
    """Set the probe MolPrep allowed atoms from all supplied target molecules."""
    raw_path = probe.get(_CONFIG_MOL_PREP)
    if not isinstance(raw_path, str) or not Path(raw_path).is_file():
        return None
    config = load_config(raw_path)
    if config.get("set_allowed_atoms_from_targets", True) is False:
        return None

    columns = config.get("columns")
    smiles_column = (
        str(columns.get("smiles", "smiles")) if isinstance(columns, dict) else "smiles"
    )
    allowed_atoms = _target_atom_symbols(molecules, smiles_column)
    filters = config.setdefault("filters", {})
    if not isinstance(filters, dict):
        raise ValueError("MolPrep filters must be a mapping.")
    filters["allowed_atoms"] = allowed_atoms
    _dump_yaml(config, Path(raw_path))
    logger.info("Target alignment MolPrep allowed atoms: %s", ", ".join(allowed_atoms))
    return allowed_atoms


def create_probe_config(
    master: dict[str, Any], target_mols_path: str, alignment_root: Path
) -> dict[str, Any]:
    """Create a copied, non-filtering config for observing target metrics."""
    probe_dir = alignment_root / PROBE_CONFIGS_DIR_NAME
    probe = _copy_master_configs(master, probe_dir)
    probe["generated_mols_path"] = str(Path(target_mols_path).resolve())
    probe["target_mols_path"] = str(Path(target_mols_path).resolve())
    probe["folder_to_save"] = str(
        (alignment_root / TARGET_CALIBRATION_RUN_DIR_NAME).resolve()
    )
    probe["sample_size"] = None
    probe["save_sampled_mols"] = True
    probe["large_dataset_mode"] = False
    # The initial descriptor stage already records every metric needed for
    # alignment. Do not calculate the same descriptors again at stage 7.
    probe[KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS] = True

    def relax_descriptors(config: dict[str, Any]) -> None:
        config["filter_data"] = True
        config["borders"] = {}
        constraints = config.get("structural_constraints")
        if isinstance(constraints, dict):
            constraints["enabled"] = False

    def enable_structural_filters(config: dict[str, Any]) -> None:
        """Evaluate every categorical rule without filtering target molecules."""
        config["run"] = True
        config["filter_data"] = False
        config["write_per_filter_outputs"] = True
        config["combine_in_memory"] = True
        config["generate_plots"] = False
        config["generate_failure_analysis"] = False
        for key in list(config):
            if key.startswith("calculate_"):
                config[key] = True

        alerts_path = config.get("alerts_data_path")
        if isinstance(alerts_path, str) and Path(alerts_path).exists():
            alerts = _read_csv(Path(alerts_path))
            if alerts is not None and "rule_set_name" in alerts.columns:
                config["include_rulesets"] = sorted(
                    alerts["rule_set_name"].dropna().astype(str).unique().tolist()
                )
        config["exclude_descriptions"] = {}

    def relax_synthesis(config: dict[str, Any]) -> None:
        # Calibration needs every configured score for every target molecule,
        # but it must not reject molecules using the source thresholds. Keep
        # those numeric thresholds visible and disable filtering explicitly.
        config["alignment_measurement_mode"] = True
        config["apply_score_filters"] = False
        config["filter_solved_only"] = False
        config["run_retrosynthesis"] = False

    def relax_docking(config: dict[str, Any]) -> None:
        calculate_thresholds = (
            config.get("calculate_score_thresholds_from_targets") is True
        )
        config["run"] = bool(config.get("run", False) and calculate_thresholds)
        config["score_thresholds"] = {}

    def disable_docking_filters(config: dict[str, Any]) -> None:
        # Target alignment uses docking only to derive raw score thresholds.
        config["run"] = False

    _update_yaml(probe.get(_CONFIG_DESCRIPTORS), relax_descriptors)
    _update_yaml(probe.get(_CONFIG_STRUCT_FILTERS), enable_structural_filters)
    _update_yaml(probe.get(_CONFIG_SYNTHESIS), relax_synthesis)
    _update_yaml(probe.get(_CONFIG_DOCKING), relax_docking)
    _update_yaml(probe.get(_CONFIG_DOCKING_FILTERS), disable_docking_filters)

    _dump_yaml(probe, probe_dir / "probe_config.yml")
    return probe


def _read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        logger.warning("Could not read alignment metrics %s: %s", path, exc)
        return None


def _read_docking_score_metrics(
    path: Path, selected_tools: list[str]
) -> pd.DataFrame | None:
    """Read one best docking score per target molecule and docking tool."""
    if not path.exists():
        return None
    rows: list[dict[str, Any]] = []
    try:
        for mol in Chem.SDMolSupplier(str(path), removeHs=False):
            if mol is None:
                continue
            tool = (
                mol.GetProp("docking_tool").strip().lower()
                if mol.HasProp("docking_tool")
                else ""
            )
            if tool not in selected_tools:
                continue
            mol_idx = ""
            for property_name in ("source_mol_idx", "mol_idx", "_Name"):
                if mol.HasProp(property_name):
                    mol_idx = mol.GetProp(property_name).strip()
                    if mol_idx:
                        break
            if not mol_idx:
                continue
            score = None
            for property_name in _DOCKING_SCORE_PROPERTIES[tool]:
                if not mol.HasProp(property_name):
                    continue
                try:
                    score = float(mol.GetProp(property_name))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(score):
                    break
                score = None
            if score is not None:
                rows.append({"mol_idx": mol_idx, "tool": tool, "score": score})
    except (OSError, RuntimeError, ValueError) as exc:
        logger.warning("Could not read alignment docking scores %s: %s", path, exc)
        return None
    if not rows:
        return None

    long_metrics = pd.DataFrame(rows)
    # Every docking score used here is an affinity where lower is better. If a
    # tool emitted multiple poses, retain its best pose for the molecule.
    long_metrics = long_metrics.groupby(
        ["mol_idx", "tool"], as_index=False, sort=False
    )["score"].min()
    metrics = long_metrics.pivot(index="mol_idx", columns="tool", values="score")

    # Retention is defined over all input target molecules, including molecules
    # for which one of the tools failed to produce a score.
    for identities_path in (
        path.parent / "input_molecules.csv",
        path.parent / "docking_results.csv",
        path.parent / "ligands.csv",
    ):
        identities = _read_csv(identities_path)
        if identities is None or "mol_idx" not in identities.columns:
            continue
        ordered_ids = list(dict.fromkeys(identities["mol_idx"].astype(str)))
        metrics = metrics.reindex(ordered_ids)
        break

    metrics.index.name = "mol_idx"
    return metrics.reset_index()


_STRUCTURAL_IDENTITY_COLUMNS = ("model_name", "mol_idx", "smiles")


def _boolean_pass_values(values: pd.Series) -> pd.Series:
    """Normalize pass columns loaded from structural-filter CSV artifacts."""
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return (
        values.fillna(False)
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes"})
    )


def _structural_identity_keys(data: pd.DataFrame) -> pd.Series:
    """Build stable keys for joining independently evaluated structural rules."""
    model = data.get("model_name", pd.Series("", index=data.index)).fillna("")
    mol_idx = data.get("mol_idx", pd.Series("", index=data.index)).fillna("")
    smiles = data.get("smiles", pd.Series("", index=data.index)).fillna("")
    return (
        model.astype(str) + "\x1f" + mol_idx.astype(str) + "\x1f" + smiles.astype(str)
    )


_STRUCTURAL_METRIC_PREFIX = "_metric_"
_RING_ALLENE = Chem.MolFromSmarts("[R]=[R]=[R]")
_DOUBLE_BOND_SMALL_RING = Chem.MolFromSmarts("[r3,r4]=[r3,r4]")


def _add_structural_numeric_metrics(
    data: pd.DataFrame, rule_names: set[str]
) -> pd.DataFrame:
    """Calculate target measurements used by parameterized structural rules."""
    score_symmetry = None
    if "symmetry" in rule_names:
        try:
            from medchem.utils.graph import score_symmetry as medchem_score_symmetry

            score_symmetry = medchem_score_symmetry
        except Exception as exc:  # pragma: no cover - optional dependency failure
            logger.warning("Could not load structural symmetry scorer: %s", exc)

    rows: list[dict[str, Any]] = []
    for raw_smiles in data.get("smiles", pd.Series("", index=data.index)):
        row = {
            "stereo_centers": math.nan,
            "stereo_undefined": math.nan,
            "halogen_F": math.nan,
            "halogen_Br": math.nan,
            "halogen_Cl": math.nan,
            "symmetry": math.nan,
            "ring_hard_failure": True,
            "ring_problem_size": math.inf,
        }
        try:
            mol = Chem.MolFromSmiles(str(raw_smiles))
        except (TypeError, ValueError):
            mol = None
        if mol is None:
            rows.append(row)
            continue

        if "stereo_center" in rule_names:
            prepared = Chem.Mol(mol)
            Chem.AssignStereochemistry(prepared, cleanIt=True, force=True)
            stereo_centers = Chem.FindMolChiralCenters(
                prepared,
                includeUnassigned=True,
                useLegacyImplementation=False,
            )
            row["stereo_centers"] = len(stereo_centers)
            row["stereo_undefined"] = sum(
                1 for _atom_index, label in stereo_centers if label == "?"
            )

        if "halogenicity" in rule_names:
            symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
            row["halogen_F"] = symbols.count("F")
            row["halogen_Br"] = symbols.count("Br")
            row["halogen_Cl"] = symbols.count("Cl")

        if "ring_infraction" in rule_names:
            row["ring_hard_failure"] = bool(
                (_RING_ALLENE is not None and mol.HasSubstructMatch(_RING_ALLENE))
                or (
                    _DOUBLE_BOND_SMALL_RING is not None
                    and mol.HasSubstructMatch(_DOUBLE_BOND_SMALL_RING)
                )
            )
            problem_sizes: list[int] = []
            for ring in mol.GetRingInfo().BondRings():
                bonds = [mol.GetBondWithIdx(index) for index in ring]
                atom_symbols = {
                    symbol
                    for bond in bonds
                    for symbol in (
                        bond.GetBeginAtom().GetSymbol(),
                        bond.GetEndAtom().GetSymbol(),
                    )
                }
                heteroatom_types = sum(
                    symbol not in {"C", "H"} for symbol in atom_symbols
                )
                if heteroatom_types > 1 or any(
                    bond.GetBondType() != Chem.BondType.SINGLE for bond in bonds
                ):
                    problem_sizes.append(len(ring))
            if problem_sizes:
                row["ring_problem_size"] = min(problem_sizes)

        if score_symmetry is not None:
            try:
                row["symmetry"] = float(score_symmetry(mol))
            except Exception as exc:
                logger.debug("Could not calculate target symmetry score: %s", exc)
        rows.append(row)

    enriched = data.copy()
    numeric = pd.DataFrame(rows, index=enriched.index)
    for column in numeric:
        enriched[f"{_STRUCTURAL_METRIC_PREFIX}{column}"] = numeric[column]
    return enriched


def _structural_rule_columns(masks: pd.DataFrame) -> list[str]:
    identity_columns = {
        column for column in _STRUCTURAL_IDENTITY_COLUMNS if column in masks.columns
    }
    return [
        column
        for column in masks.columns
        if column not in identity_columns
        and not column.startswith(_STRUCTURAL_METRIC_PREFIX)
    ]


def _read_structural_rule_masks(stage_dir: Path) -> pd.DataFrame | None:
    """Collect per-molecule pass masks for all structural rules that ran."""
    rule_masks: dict[str, pd.Series] = {}
    identities: dict[str, dict[str, Any]] = {}

    for extended_path in sorted(stage_dir.glob("*/extended.csv")):
        data = _read_csv(extended_path)
        if data is None or data.empty:
            continue
        keys = _structural_identity_keys(data)
        for row_index, identity_key in keys.items():
            identities.setdefault(
                identity_key,
                {
                    column: data.at[row_index, column] if column in data else ""
                    for column in _STRUCTURAL_IDENTITY_COLUMNS
                },
            )

        filter_name = extended_path.parent.name
        if filter_name == "common_alerts":
            pass_columns = [
                column
                for column in data.columns
                if column.startswith("pass_") and column != "pass_any"
            ]
            for column in pass_columns:
                rule_name = f"common_alerts:{column.removeprefix('pass_')}"
                mask = pd.Series(
                    _boolean_pass_values(data[column]).to_numpy(),
                    index=keys.to_numpy(),
                )
                rule_masks[rule_name] = mask.groupby(level=0).last()
            continue

        pass_column = "pass" if "pass" in data.columns else "pass_filter"
        if pass_column not in data.columns:
            continue
        mask = pd.Series(
            _boolean_pass_values(data[pass_column]).to_numpy(),
            index=keys.to_numpy(),
        )
        rule_masks[filter_name] = mask.groupby(level=0).last()

    if not rule_masks or not identities:
        return None

    identity_frame = pd.DataFrame.from_dict(identities, orient="index")
    for rule_name, mask in rule_masks.items():
        identity_frame[rule_name] = mask.reindex(identity_frame.index).fillna(False)
    return _add_structural_numeric_metrics(identity_frame, set(rule_masks)).reset_index(
        drop=True
    )


def _align_structural_numeric_parameters(
    config: dict[str, Any],
    masks: pd.DataFrame,
    percentile: float,
) -> dict[str, int | float]:
    """Calibrate parameterized rule limits and refresh their pass masks."""
    parameters: dict[str, int | float] = {}
    rule_columns = set(_structural_rule_columns(masks))

    stereo_metrics = {
        f"{_STRUCTURAL_METRIC_PREFIX}stereo_centers": "stereo_max_centers",
        f"{_STRUCTURAL_METRIC_PREFIX}stereo_undefined": "stereo_max_undefined",
    }
    if "stereo_center" in rule_columns:
        retained = _select_stage_subset(
            masks,
            {column: {"max"} for column in stereo_metrics},
            percentile,
        )
        for column, key in stereo_metrics.items():
            values = _numeric_values(retained, column)
            if key not in config or values is None:
                continue
            threshold = math.ceil(float(values.max())) + 1
            config[key] = threshold
            parameters[key] = threshold
        if all(key in parameters for key in stereo_metrics.values()):
            centers = pd.to_numeric(
                masks[f"{_STRUCTURAL_METRIC_PREFIX}stereo_centers"],
                errors="coerce",
            )
            undefined = pd.to_numeric(
                masks[f"{_STRUCTURAL_METRIC_PREFIX}stereo_undefined"],
                errors="coerce",
            )
            masks["stereo_center"] = (centers < parameters["stereo_max_centers"]) & (
                undefined < parameters["stereo_max_undefined"]
            )

    halogen_metrics = {
        f"{_STRUCTURAL_METRIC_PREFIX}halogen_F": "halogenicity_thresh_F",
        f"{_STRUCTURAL_METRIC_PREFIX}halogen_Br": "halogenicity_thresh_Br",
        f"{_STRUCTURAL_METRIC_PREFIX}halogen_Cl": "halogenicity_thresh_Cl",
    }
    if "halogenicity" in rule_columns:
        retained = _select_stage_subset(
            masks,
            {column: {"max"} for column in halogen_metrics},
            percentile,
        )
        for column, key in halogen_metrics.items():
            if key not in config:
                continue
            threshold = _observed_bound(retained, column, "max")
            if threshold is None:
                continue
            config[key] = threshold
            parameters[key] = threshold
        if all(key in parameters for key in halogen_metrics.values()):
            masks["halogenicity"] = pd.Series(True, index=masks.index)
            for column, key in halogen_metrics.items():
                values = pd.to_numeric(masks[column], errors="coerce")
                masks["halogenicity"] &= values <= parameters[key]

    symmetry_column = f"{_STRUCTURAL_METRIC_PREFIX}symmetry"
    if "symmetry" in rule_columns and "symmetry_threshold" in config:
        retained = _select_stage_subset(
            masks,
            {symmetry_column: {"max"}},
            percentile,
        )
        threshold = _observed_bound(retained, symmetry_column, "max")
        if threshold is not None:
            config["symmetry_threshold"] = threshold
            parameters["symmetry_threshold"] = threshold
            symmetry = pd.to_numeric(masks[symmetry_column], errors="coerce")
            masks["symmetry"] = symmetry <= threshold

    ring_size_column = f"{_STRUCTURAL_METRIC_PREFIX}ring_problem_size"
    ring_hard_column = f"{_STRUCTURAL_METRIC_PREFIX}ring_hard_failure"
    ring_key = "ring_infraction_hetcycle_min_size"
    if "ring_infraction" in rule_columns and ring_key in config:
        required_count = math.ceil(len(masks) * percentile / 100.0)
        problem_sizes = pd.to_numeric(masks[ring_size_column], errors="coerce").fillna(
            math.inf
        )
        hard_failures = _boolean_pass_values(masks[ring_hard_column])
        finite_sizes = problem_sizes[problem_sizes.map(math.isfinite)]
        largest_size = int(finite_sizes.max()) if not finite_sizes.empty else 0
        threshold = 0
        ring_pass = (~hard_failures) & (problem_sizes > threshold)
        for candidate in range(largest_size + 1):
            proposed = (~hard_failures) & (problem_sizes > candidate)
            if int(proposed.sum()) >= required_count:
                threshold = candidate
                ring_pass = proposed
        config[ring_key] = threshold
        parameters[ring_key] = threshold
        masks["ring_infraction"] = ring_pass

    return parameters


def _align_structural_filter_config(
    config: dict[str, Any],
    masks: pd.DataFrame | None,
    percentile: float,
) -> dict[str, Any]:
    """Enable the largest low-impact greedy rule set meeting target retention."""
    if masks is None or masks.empty:
        return {}
    parameters = _align_structural_numeric_parameters(config, masks, percentile)
    rule_columns = _structural_rule_columns(masks)
    required_count = math.ceil(len(masks) * percentile / 100.0)
    retained = pd.Series(True, index=masks.index)
    selected: list[str] = []
    rule_audit: dict[str, Any] = {}

    ordered_rules = sorted(
        rule_columns,
        key=lambda rule: (int((~masks[rule].astype(bool)).sum()), rule.lower()),
    )
    for rule in ordered_rules:
        pass_mask = masks[rule].astype(bool)
        proposed = retained & pass_mask
        enabled = int(proposed.sum()) >= required_count
        if enabled:
            retained = proposed
            selected.append(rule)
        rule_audit[rule] = {
            "enabled": enabled,
            "failed_molecules": int((~pass_mask).sum()),
            "failed_percent": float((~pass_mask).mean() * 100.0),
            "combined_retained_if_enabled": int(proposed.sum()),
        }

    calculation_keys = {
        key.removeprefix("calculate_").lower(): key
        for key in config
        if key.startswith("calculate_")
    }
    for key in calculation_keys.values():
        config[key] = False

    common_rulesets: list[str] = []
    for rule in selected:
        if rule.startswith("common_alerts:"):
            common_rulesets.append(rule.split(":", 1)[1])
            continue
        config_key = calculation_keys.get(rule.lower())
        if config_key is not None:
            config[config_key] = True

    common_key = calculation_keys.get("common_alerts")
    if common_key is not None:
        config[common_key] = bool(common_rulesets)
    config["include_rulesets"] = sorted(common_rulesets)
    config["exclude_descriptions"] = {}
    config["run"] = True

    return {
        "target_molecules": len(masks),
        "required_retained_molecules": required_count,
        "retained_molecules": int(retained.sum()),
        "retained_percent": float(retained.mean() * 100.0),
        "parameters": parameters,
        "enabled_rules": selected,
        "disabled_rules": [rule for rule in ordered_rules if rule not in selected],
        "rules": rule_audit,
    }


def _write_structural_failure_audit(masks: pd.DataFrame, path: Path) -> None:
    """Write one row for every target molecule/rule failure."""
    identity_columns = [
        column for column in _STRUCTURAL_IDENTITY_COLUMNS if column in masks.columns
    ]
    rule_columns = _structural_rule_columns(masks)
    rows: list[dict[str, Any]] = []
    for rule in rule_columns:
        failed = masks.loc[~masks[rule].astype(bool), identity_columns]
        for record in failed.to_dict(orient="records"):
            rows.append({**record, "rule": rule})
    pd.DataFrame(rows, columns=[*identity_columns, "rule"]).to_csv(path, index=False)


def _numeric_values(df: pd.DataFrame | None, column: str) -> pd.Series | None:
    if df is None or column not in df.columns:
        return None
    values = pd.to_numeric(df[column], errors="coerce").dropna()
    if values.empty:
        return None
    return values


def _ring_size_extrema(df: pd.DataFrame | None, side: str) -> pd.Series | None:
    """Return one per-molecule ring-size extreme without changing list metrics."""
    if df is None or "ring_size" not in df.columns:
        return None

    def row_extreme(raw: Any) -> float:
        values = raw
        if isinstance(raw, str):
            try:
                values = ast.literal_eval(raw)
            except (SyntaxError, ValueError):
                return math.nan
        if not isinstance(values, (list, tuple)):
            return math.nan
        numeric = pd.to_numeric(pd.Series(values, dtype="object"), errors="coerce")
        numeric = numeric.dropna()
        if numeric.empty:
            return math.nan
        return float(numeric.min() if side == "min" else numeric.max())

    result = df["ring_size"].map(row_extreme)
    return result if result.notna().any() else None


def _yaml_number(
    value: Any,
    *,
    outward: str | None = None,
    integer_values: bool = False,
) -> int | float:
    """Normalize generated thresholds without narrowing target retention."""
    number = float(value)
    if integer_values:
        return math.ceil(number)

    decimal_value = Decimal(str(number))
    quantum = Decimal("0.01")
    if outward == "min":
        rounded = float(decimal_value.quantize(quantum, rounding=ROUND_FLOOR))
    elif outward == "max":
        rounded = float(decimal_value.quantize(quantum, rounding=ROUND_CEILING))
    else:
        rounded = round(number, 2)
    return 0.0 if rounded == 0 else rounded


def _add_threshold_side(specs: dict[str, set[str]], column: str, side: str) -> None:
    specs.setdefault(column, set()).add(side)


def _select_stage_subset(
    metrics: pd.DataFrame | None,
    specs: dict[str, set[str]],
    percentile: float,
) -> pd.DataFrame | None:
    """Select one shared subset used to derive every threshold in a stage.

    Rows are ranked by their worst normalized extremeness across the configured
    metrics. Two-sided metrics favor central values, while one-sided metrics
    favor the configured better direction. Deriving all bounds from the same
    subset guarantees that at least the requested discrete target count passes.
    """
    if metrics is None or metrics.empty:
        return None
    usable = {
        column: sides
        for column, sides in specs.items()
        if (
            _numeric_values(metrics, column) is not None
            or (
                column == "ring_size"
                and _ring_size_extrema(metrics, "min") is not None
                and _ring_size_extrema(metrics, "max") is not None
            )
        )
    }
    if not usable:
        return None

    worst = pd.Series(0.0, index=metrics.index)
    for column, sides in usable.items():
        if column == "ring_size" and sides == {"min", "max"}:
            minimums = _ring_size_extrema(metrics, "min")
            maximums = _ring_size_extrema(metrics, "max")
            if minimums is None or maximums is None:
                continue
            valid_count = int((minimums.notna() & maximums.notna()).sum())
            if valid_count <= 1:
                penalty = minimums * 0.0 + 0.5
            else:
                min_ranks = (minimums.rank(method="average") - 1.0) / (
                    valid_count - 1.0
                )
                max_ranks = (maximums.rank(method="average") - 1.0) / (
                    valid_count - 1.0
                )
                penalty = pd.concat([1.0 - min_ranks, max_ranks], axis=1).max(axis=1)
            # Acyclic molecules have an empty ring-size list and pass the runtime
            # all-rings check, so missing ring extrema carry no penalty.
            worst = pd.concat([worst, penalty.fillna(0.0)], axis=1).max(axis=1)
            continue

        values = pd.to_numeric(metrics[column], errors="coerce")
        raw_ranks = values.rank(method="average")
        valid_count = int(values.notna().sum())
        if valid_count <= 1:
            ranks = raw_ranks * 0.0 + 0.5
        else:
            ranks = (raw_ranks - 1.0) / (valid_count - 1.0)
        if sides == {"min", "max"}:
            penalty = (ranks - 0.5).abs() * 2.0
        elif "min" in sides:
            penalty = 1.0 - ranks
        else:
            penalty = ranks
        worst = pd.concat([worst, penalty.fillna(0.0)], axis=1).max(axis=1)

    keep_count = max(1, math.ceil(len(metrics) * percentile / 100.0))
    selected_index = worst.sort_values(kind="stable").index[:keep_count]
    return metrics.loc[selected_index].copy()


def _observed_bound(
    retained: pd.DataFrame | None, column: str, side: str
) -> int | float | None:
    values = (
        _ring_size_extrema(retained, side)
        if column == "ring_size"
        else _numeric_values(retained, column)
    )
    if values is None:
        return None
    values = values.dropna()
    if values.empty:
        return None
    integer_values = bool(((values % 1) == 0).all())
    raw = values.min() if side == "min" else values.max()
    return _yaml_number(
        raw,
        outward=side,
        integer_values=integer_values,
    )


def _align_molprep_config(
    config: dict[str, Any], metrics: pd.DataFrame | None, _percentile: float
) -> dict[str, Any]:
    """Write the complete target atom union into MolPrep allowed_atoms."""
    if config.get("set_allowed_atoms_from_targets", True) is False:
        return {}
    if metrics is None or "atomic_symbol" not in metrics.columns:
        return {}
    symbols = {str(value) for value in metrics["atomic_symbol"].dropna() if value}
    if not symbols:
        return {}
    periodic_table = Chem.GetPeriodicTable()
    allowed_atoms = sorted(
        symbols, key=lambda symbol: periodic_table.GetAtomicNumber(symbol)
    )
    filters = config.setdefault("filters", {})
    if not isinstance(filters, dict):
        raise ValueError("MolPrep filters must be a mapping.")
    filters["allowed_atoms"] = allowed_atoms
    return {"filters.allowed_atoms": allowed_atoms}


def _align_descriptor_config(
    config: dict[str, Any], metrics: pd.DataFrame | None, percentile: float
) -> dict[str, Any]:
    aligned: dict[str, Any] = {}
    borders = config.get("borders")
    if not isinstance(borders, dict):
        return aligned

    constraints = config.get("structural_constraints")
    structural_range_columns: set[str] = set()
    if isinstance(constraints, dict):
        type_limits = constraints.get("type_limits")
        if isinstance(type_limits, dict):
            structural_range_columns.update(str(key) for key in type_limits)

        element_limits = constraints.get("element_limits")
        if isinstance(element_limits, dict):
            structural_range_columns.update(
                column
                for element, column in _STRUCTURAL_ELEMENT_COLUMNS.items()
                if element in element_limits
            )

        structural_range_columns.update(
            column
            for key, column in _STRUCTURAL_DIRECT_COLUMNS.items()
            if key in constraints
        )

    specs: dict[str, set[str]] = {}
    for key in borders:
        if key.endswith("_min"):
            _add_threshold_side(specs, key[: -len("_min")], "min")
        elif key.endswith("_max"):
            _add_threshold_side(specs, key[: -len("_max")], "max")
    for column in structural_range_columns:
        _add_threshold_side(specs, column, "min")
        _add_threshold_side(specs, column, "max")

    # Select one common reference subset for the complete descriptor stage,
    # then derive every range from that same subset.
    retained = _select_stage_subset(metrics, specs, percentile)
    for key in list(borders):
        if key.endswith("_min"):
            column, side = key[: -len("_min")], "min"
        elif key.endswith("_max"):
            column, side = key[: -len("_max")], "max"
        else:
            continue
        value = _observed_bound(retained, column, side)
        if value is None:
            continue
        borders[key] = value
        aligned[key] = value

    for column in sorted(structural_range_columns):
        for side in ("min", "max"):
            value = _observed_bound(retained, column, side)
            if value is None:
                continue
            key = f"{column}_{side}"
            borders[key] = value
            aligned[key] = value

    if isinstance(constraints, dict):
        constraints["enabled"] = False
        aligned["structural_constraints.enabled"] = False
    return aligned


def _align_synthesis_config(
    config: dict[str, Any], metrics: pd.DataFrame | None, percentile: float
) -> dict[str, Any]:
    aligned: dict[str, Any] = {}
    enabled_columns = _enabled_synthesis_score_columns(config)
    specs: dict[str, set[str]] = {}
    for column, (min_key, max_key) in _SYNTHESIS_LEGACY_FILTERS.items():
        if column not in enabled_columns:
            continue
        if min_key in config:
            _add_threshold_side(specs, column, "min")
        if max_key in config:
            _add_threshold_side(specs, column, "max")

    nested = config.get("score_filters")
    if isinstance(nested, dict):
        for column, thresholds in nested.items():
            if column not in enabled_columns or not isinstance(thresholds, dict):
                continue
            if "min" in thresholds:
                _add_threshold_side(specs, column, "min")
            if "max" in thresholds:
                _add_threshold_side(specs, column, "max")

    retained = _select_stage_subset(metrics, specs, percentile)
    for column, (min_key, max_key) in _SYNTHESIS_LEGACY_FILTERS.items():
        if column not in enabled_columns:
            continue
        if min_key in config:
            value = _observed_bound(retained, column, "min")
            config[min_key] = value
            if value is not None:
                aligned[min_key] = value
        if max_key in config:
            value = _observed_bound(retained, column, "max")
            config[max_key] = value
            if value is not None:
                aligned[max_key] = value

    if isinstance(nested, dict):
        for column, thresholds in nested.items():
            if column not in enabled_columns or not isinstance(thresholds, dict):
                continue
            applied: dict[str, int | float] = {}
            for side in ("min", "max"):
                if side not in thresholds:
                    continue
                value = _observed_bound(retained, column, side)
                thresholds[side] = value
                if value is not None:
                    applied[side] = value
            if applied:
                aligned[column] = applied
    return aligned


def _align_docking_config(
    config: dict[str, Any], metrics: pd.DataFrame | None, percentile: float
) -> dict[str, Any]:
    """Derive one affinity cutoff per configured tool from one shared subset."""
    if config.get("calculate_score_thresholds_from_targets") is not True:
        return {}
    if metrics is None or metrics.empty:
        return {}

    selected_tools = [
        tool
        for tool in _parse_tools_config(config)
        if tool in _DOCKING_SCORE_PROPERTIES
    ]
    if not selected_tools:
        return {}

    target_count = len(metrics)
    required_count = math.ceil(target_count * percentile / 100.0)
    numeric = pd.DataFrame(
        {
            tool: pd.to_numeric(metrics.get(tool), errors="coerce")
            for tool in selected_tools
        },
        index=metrics.index,
    )
    available_by_tool = {
        tool: int(numeric[tool].notna().sum()) for tool in selected_tools
    }
    complete_mask = numeric.notna().all(axis=1)
    complete_count = int(complete_mask.sum())
    if complete_count < required_count:
        raise ValueError(
            "Docking target alignment cannot retain "
            f"{required_count}/{target_count} molecules: only {complete_count} have "
            "scores from every configured docking tool "
            f"({available_by_tool})."
        )

    complete = numeric.loc[complete_mask]
    penalties: list[pd.Series] = []
    for tool in selected_tools:
        values = complete[tool]
        valid_count = len(values)
        if valid_count <= 1:
            penalty = pd.Series(0.0, index=values.index)
        else:
            penalty = (values.rank(method="average") - 1.0) / (valid_count - 1.0)
        penalties.append(penalty.rename(tool))
    worst = pd.concat(penalties, axis=1).max(axis=1)
    selected_index = worst.sort_values(kind="stable").index[:required_count]
    retained_subset = complete.loc[selected_index]

    thresholds: dict[str, dict[str, str | int | float]] = {}
    pass_mask = pd.Series(True, index=numeric.index)
    for tool in selected_tools:
        maximum = _observed_bound(retained_subset, tool, "max")
        if maximum is None:
            continue
        thresholds[tool] = {
            "score_property": "minimizedAffinity",
            "max": maximum,
        }
        pass_mask &= numeric[tool].notna() & (numeric[tool] <= float(maximum))

    config["score_thresholds"] = thresholds
    retained_count = int(pass_mask.sum())
    return {
        "target_molecules": target_count,
        "required_retained_molecules": required_count,
        "retained_molecules": retained_count,
        "retained_percent": float(retained_count / target_count * 100.0),
        "combination": "all_configured_tools_must_pass",
        "available_scores": available_by_tool,
        "score_thresholds": thresholds,
    }


def _new_threshold_summary(
    target_run: Path,
    target_mols_path: str,
    percentile: float,
) -> dict[str, Any]:
    """Build the audit document used by incremental stage alignment."""
    return {
        "target_mols_path": str(Path(target_mols_path).resolve()),
        "target_run": str(target_run.resolve()),
        "target_coverage_percent": percentile,
        "selection_method": "stage_specific_target_coverage",
        "stages": {
            "mol_prep": {
                "thresholds": {},
                "status": "pending_atoms",
                "note": "Allowed atoms are derived from all target molecules.",
            },
            "descriptors": {
                "thresholds": {},
                "status": "pending_metrics",
                "note": "Min/max ranges use one shared stage-level reference subset.",
            },
            "struct_filters": {
                "thresholds": {},
                "status": "pending_metrics",
            },
            "synthesis": {"thresholds": {}, "status": "pending_metrics"},
            "docking": {"thresholds": {}, "status": "pending_metrics"},
            "docking_filters": {
                "thresholds": {},
                "status": "not_aligned",
                "note": "Docking alignment calibrates raw score thresholds only.",
            },
            "final_descriptors": {
                "thresholds": {},
                "status": "skipped_redundant",
                "note": "Initial descriptor measurements are reused.",
            },
        },
    }


def create_aligned_stage_config(
    master: dict[str, Any],
    target_run: Path,
    alignment_root: Path,
    target_mols_path: str,
    percentile: float,
    stage: str,
) -> tuple[dict[str, Any], Path, Path] | None:
    """Create only the aligned config belonging to one completed target stage."""
    percentile = validate_target_coverage_percent(percentile)
    aligned_dir = alignment_root / "aligned_configs"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    master_path = aligned_dir / ALIGNED_CONFIG_NAME
    thresholds_path = aligned_dir / THRESHOLDS_NAME

    aligned = (
        load_config(str(master_path)) if master_path.exists() else copy.deepcopy(master)
    )
    summary = (
        load_config(str(thresholds_path))
        if thresholds_path.exists()
        else _new_threshold_summary(target_run, target_mols_path, percentile)
    )
    summary.pop("retention_percentile", None)
    summary["target_coverage_percent"] = percentile

    if stage == "mol_prep":
        molprep_config_path = master.get(_CONFIG_MOL_PREP)
        molprep_config = (
            load_config(str(molprep_config_path))
            if isinstance(molprep_config_path, str)
            and Path(molprep_config_path).is_file()
            else {}
        )
        if molprep_config.get("set_allowed_atoms_from_targets", True) is False:
            summary["stages"][stage]["status"] = "disabled_by_config"
            _dump_yaml(summary, thresholds_path)
            return None
        target_molecules = _read_csv(target_run / "input" / "sampled_molecules.csv")
        if target_molecules is None:
            metrics = None
        else:
            allowed_atoms = _target_atom_symbols(target_molecules)
            metrics = pd.DataFrame({"atomic_symbol": allowed_atoms})
        specs = [(_CONFIG_MOL_PREP, _align_molprep_config)]
    elif stage == "descriptors":
        metrics = _read_csv(
            target_run
            / "stages"
            / "02_descriptors_initial"
            / "metrics"
            / "descriptors_all.csv"
        )
        specs = [(_CONFIG_DESCRIPTORS, _align_descriptor_config)]
    elif stage == "struct_filters":
        metrics = _read_structural_rule_masks(
            target_run / "stages" / "03_structural_filters_post"
        )
        specs = [(_CONFIG_STRUCT_FILTERS, _align_structural_filter_config)]
    elif stage == "synthesis":
        metrics = _read_csv(
            target_run / "stages" / "04_synthesis" / "synthesis_scores.csv"
        )
        specs = [(_CONFIG_SYNTHESIS, _align_synthesis_config)]
    elif stage == "docking":
        docking_config_path = master.get(_CONFIG_DOCKING)
        docking_config = (
            load_config(str(docking_config_path))
            if isinstance(docking_config_path, str)
            and Path(docking_config_path).is_file()
            else {}
        )
        if docking_config.get("calculate_score_thresholds_from_targets") is not True:
            summary["stages"][stage]["status"] = "disabled_by_config"
            summary["stages"][stage]["note"] = (
                "Set calculate_score_thresholds_from_targets: true in the "
                "docking config to calibrate docking scores."
            )
            _dump_yaml(summary, thresholds_path)
            return None
        metrics = _read_docking_score_metrics(
            target_run / "stages" / "05_docking" / "docking_out.sdf",
            _parse_tools_config(docking_config),
        )
        specs = [(_CONFIG_DOCKING, _align_docking_config)]
    else:
        return None

    if metrics is None:
        summary["stages"][stage]["status"] = "missing_metrics"
        _dump_yaml(summary, thresholds_path)
        logger.warning("No alignment metrics found after target stage '%s'.", stage)
        return None

    for config_key, updater in specs:
        raw_path = master.get(config_key)
        if not isinstance(raw_path, str):
            continue
        source = Path(raw_path)
        if not source.exists() or not source.is_file():
            continue
        target = aligned_dir / f"{config_key}{source.suffix or '.yml'}"
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
        config = load_config(str(target))
        thresholds = updater(config, metrics, percentile)
        if stage == "struct_filters":
            failure_audit_path = aligned_dir / "structural_filter_failures.csv"
            _write_structural_failure_audit(metrics, failure_audit_path)
            thresholds["failure_audit_path"] = str(failure_audit_path.resolve())
        _dump_yaml(config, target)
        aligned[config_key] = str(target.resolve())
        summary["stages"][stage]["thresholds"] = thresholds
        summary["stages"][stage]["status"] = "ready"

    aligned["target_mols_path"] = str(Path(target_mols_path).resolve())
    aligned["alignment"] = {
        "enabled": False,
        "target_coverage_percent": percentile,
        "thresholds_path": str(thresholds_path.resolve()),
        "target_run": str(target_run.resolve()),
    }
    _dump_yaml(summary, thresholds_path)
    _dump_yaml(aligned, master_path)
    return aligned, master_path, thresholds_path
