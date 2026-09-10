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
from yaml.nodes import MappingNode, ScalarNode, SequenceNode

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
SOURCE_CONFIGS_DIR_NAME = "source_configs"
TARGET_CALIBRATION_RUN_DIR_NAME = "calibration_target_run"
SOURCE_MASTER_PATH_KEY = "_source_master_path"
_RUNTIME_STAGE_OVERRIDE_KEYS = (
    "_run_single_stage_override",
    "_run_stage_selection_override",
)

DESCRIPTOR_BOUNDS_MODE_TARGET = "target"
DESCRIPTOR_BOUNDS_MODE_EXPAND = "expand"
DEFAULT_DESCRIPTOR_BOUNDS_MODE = DESCRIPTOR_BOUNDS_MODE_EXPAND
DESCRIPTOR_BOUNDS_MODES = {
    DESCRIPTOR_BOUNDS_MODE_TARGET,
    DESCRIPTOR_BOUNDS_MODE_EXPAND,
}

_CONFIG_MOL_PREP = "config_mol_prep"
_CONFIG_DESCRIPTORS = "config_descriptors"
_CONFIG_STRUCT_FILTERS = "config_structFilters"
_CONFIG_SYNTHESIS = "config_synthesis"
_CONFIG_DOCKING = "config_docking"
_CONFIG_DOCKING_FILTERS = "config_docking_filters"

_DOCKING_SCORE_PROPERTIES = {
    TOOL_SMINA: ("minimizedAffinity", "affinity", "score"),
    TOOL_GNINA: ("minimizedAffinity", "affinity", "score"),
    TOOL_MATCHA: ("minimizedAffinity", "affinity", "score"),
}
_DEFAULT_DOCKING_SCORE_PROPERTY = "minimizedAffinity"


def _configured_docking_score_properties(selected_tools: list[str]) -> dict[str, str]:
    return {tool: _DEFAULT_DOCKING_SCORE_PROPERTY for tool in selected_tools}


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


class GeneratedConfigShapeError(ValueError):
    """Raised when generated data would change a source config's schema."""


def _source_shaped_value(source: Any, updated: Any, path: tuple[str, ...] = ()) -> Any:
    """Project updated values onto the exact mapping shape of the source YAML."""
    if isinstance(source, dict):
        if not isinstance(updated, dict):
            raise GeneratedConfigShapeError(
                f"Generated config changes mapping {'.'.join(path) or '<root>'} "
                f"into {type(updated).__name__}."
            )
        # An explicitly empty mapping is a supported placeholder for generated
        # values such as score_thresholds.
        if not source:
            return copy.deepcopy(updated)
        return {
            key: _source_shaped_value(
                source_value,
                updated.get(key, source_value),
                (*path, str(key)),
            )
            for key, source_value in source.items()
        }
    if isinstance(source, list):
        if not isinstance(updated, list):
            raise GeneratedConfigShapeError(
                f"Generated config changes list {'.'.join(path)} into "
                f"{type(updated).__name__}."
            )
        return copy.deepcopy(updated)
    if isinstance(updated, (dict, list)):
        raise GeneratedConfigShapeError(
            f"Generated config changes scalar {'.'.join(path)} into a collection."
        )
    return copy.deepcopy(updated)


def _render_yaml_node(value: Any, node: ScalarNode | SequenceNode | MappingNode) -> str:
    """Render one replacement value without rewriting the surrounding YAML."""
    flow_style = (
        bool(node.flow_style) if isinstance(node, (SequenceNode, MappingNode)) else None
    )
    rendered = yaml.safe_dump(
        value,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=flow_style,
        width=4096,
    )
    lines = rendered.rstrip().splitlines()
    if lines and lines[-1] == "...":
        lines.pop()
    if not lines:
        return "null"
    padding = " " * node.start_mark.column
    replacement = ("\n" + padding).join(lines)
    if isinstance(node, (SequenceNode, MappingNode)) and not node.flow_style:
        replacement += "\n"
    return replacement


def _yaml_replacements(
    source: Any,
    updated: Any,
    node: ScalarNode | SequenceNode | MappingNode,
    path: tuple[str, ...] = (),
) -> list[tuple[int, int, str]]:
    """Return source-text slices needed to materialize updated values."""
    if isinstance(node, MappingNode):
        if not isinstance(source, dict) or not isinstance(updated, dict):
            return [
                (
                    node.start_mark.index,
                    node.end_mark.index,
                    _render_yaml_node(updated, node),
                )
            ]
        source_keys = list(source)
        updated_keys = list(updated)
        if source_keys != updated_keys:
            if source:
                raise GeneratedConfigShapeError(
                    f"Generated config changes keys under {'.'.join(path) or '<root>'}."
                )
            return [
                (
                    node.start_mark.index,
                    node.end_mark.index,
                    _render_yaml_node(updated, node),
                )
            ]
        replacements: list[tuple[int, int, str]] = []
        for (_key_node, value_node), key in zip(node.value, source_keys, strict=True):
            replacements.extend(
                _yaml_replacements(
                    source[key],
                    updated[key],
                    value_node,
                    (*path, str(key)),
                )
            )
        return replacements

    if isinstance(node, SequenceNode):
        if not isinstance(source, list) or not isinstance(updated, list):
            return [
                (
                    node.start_mark.index,
                    node.end_mark.index,
                    _render_yaml_node(updated, node),
                )
            ]
        if len(source) != len(updated):
            return [
                (
                    node.start_mark.index,
                    node.end_mark.index,
                    _render_yaml_node(updated, node),
                )
            ]
        replacements = []
        for index, (source_item, updated_item, item_node) in enumerate(
            zip(source, updated, node.value, strict=True)
        ):
            replacements.extend(
                _yaml_replacements(
                    source_item,
                    updated_item,
                    item_node,
                    (*path, str(index)),
                )
            )
        return replacements

    if not isinstance(node, ScalarNode):
        raise GeneratedConfigShapeError(
            f"Unsupported YAML node under {'.'.join(path) or '<root>'}."
        )
    if source == updated and type(source) is type(updated):
        return []
    return [
        (
            node.start_mark.index,
            node.end_mark.index,
            _render_yaml_node(updated, node),
        )
    ]


def _dump_generated_yaml(
    data: dict[str, Any],
    target: Path,
    *,
    source: Path | None,
) -> dict[str, Any]:
    """Write a generated config by patching values into source YAML text."""
    if source is None or not source.is_file():
        _dump_yaml(data, target)
        return copy.deepcopy(data)

    source_text = source.read_text(encoding="utf-8")
    source_data = yaml.safe_load(source_text)
    source_node = yaml.compose(source_text)
    if not isinstance(source_data, dict) or not isinstance(source_node, MappingNode):
        raise GeneratedConfigShapeError(f"Source config must be a mapping: {source}")

    shaped = _source_shaped_value(source_data, data)
    replacements = _yaml_replacements(source_data, shaped, source_node)
    generated = source_text
    for start, end, replacement in sorted(replacements, reverse=True):
        # PyYAML includes the indentation before the next mapping key in the
        # end mark of an indentless block sequence. Keep that indentation in
        # the source text instead of consuming it with the replaced list.
        while end > start and source_text[end - 1] in " \t":
            end -= 1
        generated = generated[:start] + replacement + generated[end:]

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(generated, encoding="utf-8")
    return shaped


def _changed_value_paths(
    source: Any,
    updated: Any,
    path: tuple[str, ...] = (),
) -> set[tuple[str, ...]]:
    """Return the smallest source paths whose rendered values would change."""
    if isinstance(source, dict) and isinstance(updated, dict):
        if not source:
            return {path} if source != updated else set()
        changed: set[tuple[str, ...]] = set()
        for key, source_value in source.items():
            changed.update(
                _changed_value_paths(
                    source_value,
                    updated.get(key, source_value),
                    (*path, str(key)),
                )
            )
        return changed
    if isinstance(source, list) and isinstance(updated, list):
        return {path} if source != updated else set()
    return {path} if source != updated or type(source) is not type(updated) else set()


def _aligned_threshold_paths(
    config_key: str,
    source: dict[str, Any],
) -> set[tuple[str, ...]]:
    """Return source fields that target alignment may replace."""
    if config_key == _CONFIG_MOL_PREP:
        return {("filters", "allowed_atoms")}
    if config_key == _CONFIG_DESCRIPTORS:
        borders = source.get("borders")
        if not isinstance(borders, dict):
            return set()
        return {
            ("borders", str(key))
            for key in borders
            if str(key).endswith(("_min", "_max"))
        }
    if config_key == _CONFIG_STRUCT_FILTERS:
        return set()
    if config_key == _CONFIG_DOCKING and isinstance(
        source.get("score_thresholds"), dict
    ):
        return {("score_thresholds",)}
    return set()


def _dump_aligned_stage_yaml(
    data: dict[str, Any],
    target: Path,
    *,
    source: Path,
    config_key: str,
) -> dict[str, Any]:
    """Copy a source stage config while replacing threshold values only."""
    source_data = load_config(str(source))
    shaped = _source_shaped_value(source_data, data)
    changed_paths = _changed_value_paths(source_data, shaped)
    allowed_paths = _aligned_threshold_paths(config_key, source_data)
    unexpected = sorted(
        changed_path
        for changed_path in changed_paths
        if not any(
            changed_path[: len(allowed_path)] == allowed_path
            for allowed_path in allowed_paths
        )
    )
    if unexpected:
        formatted = ", ".join(".".join(path) or "<root>" for path in unexpected)
        raise GeneratedConfigShapeError(
            f"Aligned {config_key} changes non-threshold source fields: {formatted}."
        )
    return _dump_generated_yaml(shaped, target, source=source)


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
    source_dir = alignment_root / SOURCE_CONFIGS_DIR_NAME
    source_master_raw = master.get(SOURCE_MASTER_PATH_KEY)
    clean_master = copy.deepcopy(master)
    clean_master.pop(SOURCE_MASTER_PATH_KEY, None)
    runtime_stage_overrides = {
        key: copy.deepcopy(clean_master[key])
        for key in _RUNTIME_STAGE_OVERRIDE_KEYS
        if key in clean_master
    }
    source = _copy_master_configs(clean_master, source_dir)
    source_master_path = source_dir / "source_config.yml"
    raw_template = (
        Path(str(source_master_raw))
        if isinstance(source_master_raw, str) and Path(str(source_master_raw)).is_file()
        else None
    )
    persisted_source = copy.deepcopy(source)
    for key in _RUNTIME_STAGE_OVERRIDE_KEYS:
        persisted_source.pop(key, None)
    source = _write_generated_master(
        persisted_source,
        source_master_path,
        source=raw_template,
    )
    # Runtime CLI stage selection is intentionally absent from reusable source
    # YAML, but the current target probe must obey the user's selected stages.
    source.update(runtime_stage_overrides)

    probe_dir = alignment_root / PROBE_CONFIGS_DIR_NAME
    probe = _copy_master_configs(source, probe_dir)
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
        config["generate_plots"] = False
        config["generate_failure_analysis"] = False
        for key in list(config):
            if key.startswith("calculate_"):
                config[key] = True

        alerts_path = config.get("alerts_data_path")
        if isinstance(alerts_path, str) and Path(alerts_path).exists():
            alerts = _read_csv(Path(alerts_path))
            if alerts is not None and "rule_set_name" in alerts.columns:
                config["include_rulesets"] = "all"
        config["exclude_smarts"] = []

    def skip_synthesis(config: dict[str, Any]) -> None:
        """Bypass synthesis during target calibration without changing its policy."""
        config["run"] = False

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
    _update_yaml(probe.get(_CONFIG_SYNTHESIS), skip_synthesis)
    _update_yaml(probe.get(_CONFIG_DOCKING), relax_docking)
    _update_yaml(probe.get(_CONFIG_DOCKING_FILTERS), disable_docking_filters)

    _drop_synthesis_bounds_mode(probe)
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
    path: Path,
    selected_tools: list[str],
    score_properties: dict[str, str] | None = None,
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
            configured_property = (score_properties or {}).get(tool)
            property_names = (
                (configured_property,)
                if configured_property
                else _DOCKING_SCORE_PROPERTIES[tool]
            )
            for property_name in property_names:
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
    # All configured docking scores use lower-is-better semantics. GNINA can
    # emit multiple poses, so retain the lowest minimizedAffinity per molecule
    # before deriving its target threshold.
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
        if filter_name == "stereo_center" and "undefined_stereo_pass" in data:
            undefined_mask = pd.Series(
                _boolean_pass_values(data["undefined_stereo_pass"]).to_numpy(),
                index=keys.to_numpy(),
            )
            rule_masks["undefined_stereo_center"] = (
                undefined_mask.groupby(level=0).last()
            )

    if not rule_masks or not identities:
        return None

    identity_frame = pd.DataFrame.from_dict(identities, orient="index")
    for rule_name, mask in rule_masks.items():
        identity_frame[rule_name] = mask.reindex(identity_frame.index).fillna(False)
    return _add_structural_numeric_metrics(identity_frame, set(rule_masks)).reset_index(
        drop=True
    )


def _configured_structural_rule_columns(
    config: dict[str, Any], masks: pd.DataFrame
) -> list[str]:
    """Return the hard structural rules enabled by the source config."""

    def calculation_enabled(policy_name: str) -> bool:
        calculation_name = (
            "stereo_center"
            if policy_name == "undefined_stereo_center"
            else policy_name
        )
        return bool(config.get(f"calculate_{calculation_name}", False))

    def hard_enabled(policy_name: str) -> bool:
        key = f"filter_{policy_name}"
        if key in config:
            return bool(config[key])
        return calculation_enabled(policy_name)

    raw_include = config.get("include_rulesets")
    if isinstance(raw_include, str) and raw_include.strip().lower() == "all":
        # Scalar "all" means every catalog ruleset is calculated.
        calculated_rulesets: set[str] = set()
        calculated_unrestricted = True
    else:
        calculated_rulesets = {
            str(value)
            for value in (raw_include or [])
            if value is not None
        }
        calculated_unrestricted = not calculated_rulesets
    filter_rulesets = {
        str(value)
        for value in config.get("common_alerts_filter_include_rulesets", []) or []
        if value is not None
    }
    excluded_rulesets = {
        str(value)
        for value in config.get("common_alerts_filter_exclude_rulesets", []) or []
        if value is not None
    }

    selected: list[str] = []
    for rule in _structural_rule_columns(masks):
        if rule.startswith("common_alerts:"):
            ruleset = rule.split(":", 1)[1]
            enabled = (
                calculation_enabled("common_alerts")
                and hard_enabled("common_alerts")
                and (calculated_unrestricted or ruleset in calculated_rulesets)
                and (not filter_rulesets or ruleset in filter_rulesets)
                and ruleset not in excluded_rulesets
            )
        else:
            enabled = calculation_enabled(rule) and hard_enabled(rule)
        if enabled:
            selected.append(rule)
    return selected


def _align_structural_filter_config(
    config: dict[str, Any],
    masks: pd.DataFrame | None,
    percentile: float,
) -> dict[str, Any]:
    """Audit the fixed source structural policy without changing its config."""
    if masks is None or masks.empty:
        return {}

    rule_columns = _configured_structural_rule_columns(config, masks)
    required_count = math.ceil(len(masks) * percentile / 100.0)
    retained = pd.Series(True, index=masks.index)
    rule_audit: dict[str, Any] = {}
    for rule in rule_columns:
        pass_mask = masks[rule].astype(bool)
        retained &= pass_mask
        rule_audit[rule] = {
            "enabled": True,
            "failed_molecules": int((~pass_mask).sum()),
            "failed_percent": float((~pass_mask).mean() * 100.0),
            "combined_retained_if_enabled": int(retained.sum()),
        }

    retained_count = int(retained.sum())
    return {
        "target_molecules": len(masks),
        "required_retained_molecules": required_count,
        "retained_molecules": retained_count,
        "retained_percent": float(retained.mean() * 100.0),
        "coverage_met": retained_count >= required_count,
        "policy": "source_config_preserved",
        "parameters": {},
        "enabled_rules": rule_columns,
        "disabled_rules": [],
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


def _validate_bounds_mode(value: object, config_key: str) -> str:
    """Validate how target-derived bounds modify source bounds."""
    if not isinstance(value, str) or value not in DESCRIPTOR_BOUNDS_MODES:
        choices = ", ".join(sorted(DESCRIPTOR_BOUNDS_MODES))
        raise ValueError(f"alignment.{config_key} must be one of: {choices}")
    return value


def validate_descriptor_bounds_mode(value: object) -> str:
    """Validate how target-derived descriptor borders modify source borders."""
    return _validate_bounds_mode(value, "descriptor_bounds_mode")


def descriptor_bounds_mode_from_master(master: dict[str, Any]) -> str:
    """Return the descriptor mode, defaulting to expand when unset."""
    alignment = master.get("alignment")
    if not isinstance(alignment, dict):
        return DEFAULT_DESCRIPTOR_BOUNDS_MODE
    return validate_descriptor_bounds_mode(
        alignment.get("descriptor_bounds_mode", DEFAULT_DESCRIPTOR_BOUNDS_MODE)
    )


def _drop_synthesis_bounds_mode(config: dict[str, Any]) -> None:
    """Remove leftover synthesis_bounds_mode so it is never written back out."""
    alignment = config.get("alignment")
    if isinstance(alignment, dict):
        alignment.pop("synthesis_bounds_mode", None)


def _drop_synthesis_bounds_mode_from_yaml_text(text: str) -> str:
    """Drop leftover synthesis_bounds_mode lines copied from a source YAML template."""
    return "".join(
        line
        for line in text.splitlines(keepends=True)
        if not line.lstrip().startswith("synthesis_bounds_mode:")
    )


def _write_generated_master(
    aligned: dict[str, Any],
    master_path: Path,
    source: Path | None,
) -> dict[str, Any]:
    """Write an aligned/probe master without leftover synthesis_bounds_mode."""
    _drop_synthesis_bounds_mode(aligned)
    written = _dump_generated_yaml(aligned, master_path, source=source)
    _drop_synthesis_bounds_mode(written)
    text = master_path.read_text(encoding="utf-8")
    stripped = _drop_synthesis_bounds_mode_from_yaml_text(text)
    if stripped != text:
        master_path.write_text(stripped, encoding="utf-8")
    return written


def _merge_numeric_bound(
    source_value: object,
    target_value: int | float,
    side: str,
    label: str,
) -> object:
    """Keep the source bound unless the target requires an outward expansion."""
    if isinstance(source_value, bool):
        raise ValueError(f"{label} bounds must be numeric, not boolean")
    try:
        source_numeric = float(source_value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} bound must be numeric: {source_value!r}") from exc
    if side == "min":
        return target_value if float(target_value) < source_numeric else source_value
    return target_value if float(target_value) > source_numeric else source_value


def _merge_descriptor_bound(
    source_value: object,
    target_value: int | float,
    side: str,
) -> object:
    """Keep the descriptor source bound unless target requires expansion."""
    return _merge_numeric_bound(source_value, target_value, side, "Descriptor")


def _align_molprep_config(
    config: dict[str, Any], metrics: pd.DataFrame | None
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
    config: dict[str, Any],
    metrics: pd.DataFrame | None,
    percentile: float,
    bounds_mode: str = DEFAULT_DESCRIPTOR_BOUNDS_MODE,
) -> dict[str, Any]:
    """Align configured borders to target values, optionally only widening them."""
    bounds_mode = validate_descriptor_bounds_mode(bounds_mode)
    aligned: dict[str, Any] = {}
    borders = config.get("borders")
    if not isinstance(borders, dict):
        return aligned

    specs: dict[str, set[str]] = {}
    for key in borders:
        if key.endswith("_min"):
            _add_threshold_side(specs, key.removesuffix("_min"), "min")
        elif key.endswith("_max"):
            _add_threshold_side(specs, key.removesuffix("_max"), "max")

    retained = _select_stage_subset(metrics, specs, percentile)
    for key in borders:
        if key.endswith("_min"):
            column, side = key.removesuffix("_min"), "min"
        elif key.endswith("_max"):
            column, side = key.removesuffix("_max"), "max"
        else:
            continue
        value = _observed_bound(retained, column, side)
        if value is None:
            continue
        if bounds_mode == DESCRIPTOR_BOUNDS_MODE_EXPAND:
            value = _merge_descriptor_bound(borders[key], value, side)
        borders[key] = value
        aligned[key] = value
    return aligned


def _align_docking_config(
    config: dict[str, Any], metrics: pd.DataFrame | None, percentile: float
) -> dict[str, Any]:
    """Derive one affinity cutoff per configured tool from one shared subset.

    Target calibration may relax a configured upper bound, but never tighten it.
    For lower-is-better affinity scores this is ``max(calibrated, configured)``.
    """
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

    configured_thresholds = config.get("score_thresholds")
    if not isinstance(configured_thresholds, dict):
        configured_thresholds = {}

    thresholds: dict[str, dict[str, str | int | float]] = {}
    calibrated_thresholds: dict[str, dict[str, str | int | float]] = {}
    pass_mask = pd.Series(True, index=numeric.index)
    for tool in selected_tools:
        calibrated_maximum = _observed_bound(retained_subset, tool, "max")
        if calibrated_maximum is None:
            continue
        maximum = calibrated_maximum
        configured_tool_threshold = configured_thresholds.get(tool)
        if isinstance(configured_tool_threshold, dict):
            configured_maximum = configured_tool_threshold.get("max")
            if isinstance(configured_maximum, (int, float)) and not isinstance(
                configured_maximum, bool
            ):
                maximum = max(float(calibrated_maximum), float(configured_maximum))
                maximum = _yaml_number(maximum)

        score_property = _DEFAULT_DOCKING_SCORE_PROPERTY
        calibrated_thresholds[tool] = {
            "score_property": score_property,
            "max": calibrated_maximum,
        }
        thresholds[tool] = {
            "score_property": score_property,
            "max": maximum,
        }
        pass_mask &= numeric[tool].notna() & (numeric[tool] <= float(maximum))

    if isinstance(config.get("score_thresholds"), dict):
        config["score_thresholds"] = thresholds
    retained_count = int(pass_mask.sum())
    return {
        "target_molecules": target_count,
        "required_retained_molecules": required_count,
        "retained_molecules": retained_count,
        "retained_percent": float(retained_count / target_count * 100.0),
        "combination": "all_configured_tools_must_pass",
        "available_scores": available_by_tool,
        "calibrated_score_thresholds": calibrated_thresholds,
        "score_thresholds": thresholds,
    }


def _config_from_master(
    master: dict[str, Any], config_key: str
) -> tuple[dict[str, Any], Path] | None:
    raw_path = master.get(config_key)
    if not isinstance(raw_path, str):
        return None
    path = Path(raw_path)
    if not path.is_file():
        return None
    return load_config(str(path)), path


def _descriptor_threshold_specs(config: dict[str, Any]) -> dict[str, set[str]]:
    specs: dict[str, set[str]] = {}
    borders = config.get("borders")
    if isinstance(borders, dict):
        for key in borders:
            if key.endswith("_min"):
                _add_threshold_side(specs, key.removesuffix("_min"), "min")
            elif key.endswith("_max"):
                _add_threshold_side(specs, key.removesuffix("_max"), "max")

    constraints = config.get("structural_constraints")
    if not isinstance(constraints, dict):
        return specs
    type_limits = constraints.get("type_limits")
    if isinstance(type_limits, dict):
        for column in type_limits:
            _add_threshold_side(specs, str(column), "min")
            _add_threshold_side(specs, str(column), "max")
    element_limits = constraints.get("element_limits")
    if isinstance(element_limits, dict):
        for element, column in _STRUCTURAL_ELEMENT_COLUMNS.items():
            if element in element_limits:
                _add_threshold_side(specs, column, "min")
                _add_threshold_side(specs, column, "max")
    for key, column in _STRUCTURAL_DIRECT_COLUMNS.items():
        if key in constraints:
            _add_threshold_side(specs, column, "min")
            _add_threshold_side(specs, column, "max")
    return specs


def _metric_extremeness(metrics: pd.DataFrame, specs: dict[str, set[str]]) -> pd.Series:
    """Return each molecule's worst normalized threshold extremeness."""
    penalties: list[pd.Series] = []
    for column, sides in specs.items():
        if column == "ring_size" and sides == {"min", "max"}:
            minimums = _ring_size_extrema(metrics, "min")
            maximums = _ring_size_extrema(metrics, "max")
            if minimums is None or maximums is None:
                continue
            valid_count = int((minimums.notna() & maximums.notna()).sum())
            if valid_count <= 1:
                penalty = pd.Series(0.5, index=metrics.index)
            else:
                min_ranks = (minimums.rank(method="average") - 1.0) / (
                    valid_count - 1.0
                )
                max_ranks = (maximums.rank(method="average") - 1.0) / (
                    valid_count - 1.0
                )
                penalty = pd.concat([1.0 - min_ranks, max_ranks], axis=1).max(axis=1)
            penalties.append(penalty.fillna(0.0).rename(column))
            continue

        if column not in metrics:
            continue
        values = pd.to_numeric(metrics[column], errors="coerce")
        valid_count = int(values.notna().sum())
        if valid_count <= 1:
            ranks = pd.Series(0.5, index=metrics.index)
        else:
            ranks = (values.rank(method="average") - 1.0) / (valid_count - 1.0)
        if sides == {"min", "max"}:
            penalty = (ranks - 0.5).abs() * 2.0
        elif "min" in sides:
            penalty = 1.0 - ranks
        else:
            penalty = ranks
        penalties.append(penalty.fillna(1.0).rename(column))
    if not penalties:
        return pd.Series(0.0, index=metrics.index)
    return pd.concat(penalties, axis=1).max(axis=1)


def _metrics_by_target_id(
    metrics: pd.DataFrame,
    target_ids: pd.Index,
    stage: str,
) -> pd.DataFrame:
    """Index one-row-per-target metrics by stable mol_idx."""
    indexed = metrics.copy()
    if "mol_idx" in indexed.columns:
        indexed["mol_idx"] = indexed["mol_idx"].astype(str)
        if indexed["mol_idx"].duplicated().any():
            raise ValueError(
                f"{stage} alignment metrics contain duplicate mol_idx values."
            )
        indexed = indexed.set_index("mol_idx", drop=False)
    elif len(indexed) == len(target_ids):
        logger.warning(
            "%s alignment metrics have no mol_idx; using preserved target row order.",
            stage,
        )
        indexed.index = target_ids
        indexed["mol_idx"] = target_ids
    else:
        raise ValueError(
            f"{stage} alignment metrics have no mol_idx and contain "
            f"{len(indexed)} rows for {len(target_ids)} targets."
        )
    return indexed.reindex(target_ids)


def _select_global_protected_cohort(
    master: dict[str, Any],
    target_run: Path,
    target_molecules: pd.DataFrame,
    percentile: float,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Select one fixed target cohort from configurable numeric score stages."""
    raw_selected_stages = master.get("_run_stage_selection_override")
    selected_stages = (
        {str(stage) for stage in raw_selected_stages}
        if isinstance(raw_selected_stages, (list, tuple, set))
        else None
    )

    def stage_selected(stage: str) -> bool:
        return selected_stages is None or stage in selected_stages

    target_ids = pd.Index(target_molecules["mol_idx"].astype(str), name="mol_idx")
    required_count = math.ceil(len(target_ids) * percentile / 100.0)
    penalties: dict[str, pd.Series] = {}
    stage_metrics: dict[str, pd.DataFrame] = {}
    eligible = pd.Series(True, index=target_ids)

    descriptor_source = _config_from_master(master, _CONFIG_DESCRIPTORS)
    if (
        stage_selected("descriptors")
        and descriptor_source is not None
        and descriptor_source[0].get("run", True)
    ):
        raw = _read_csv(
            target_run
            / "stages"
            / "02_descriptors_initial"
            / "metrics"
            / "descriptors_all.csv"
        )
        if raw is None:
            raise ValueError("Descriptor target metrics are missing.")
        metrics = _metrics_by_target_id(raw, target_ids, "Descriptor")
        specs = _descriptor_threshold_specs(descriptor_source[0])
        penalties["descriptors"] = _metric_extremeness(metrics, specs)
        stage_metrics["descriptors"] = metrics

    structural_source = _config_from_master(master, _CONFIG_STRUCT_FILTERS)
    if (
        stage_selected("struct_filters")
        and structural_source is not None
        and structural_source[0].get("run", False)
    ):
        raw = _read_structural_rule_masks(
            target_run / "stages" / "03_structural_filters_post"
        )
        if raw is None:
            raise ValueError("Structural-filter target metrics are missing.")
        # Structural filters are a fixed policy. Keep their measurements for
        # diagnostics, but never use them to choose the percentile cohort.
        stage_metrics["struct_filters"] = _metrics_by_target_id(
            raw, target_ids, "Structural-filter"
        )

    docking_source = _config_from_master(master, _CONFIG_DOCKING)
    if (
        stage_selected("docking")
        and docking_source is not None
        and docking_source[0].get("run", False)
        and docking_source[0].get("calculate_score_thresholds_from_targets") is True
    ):
        selected_tools = [
            tool
            for tool in _parse_tools_config(docking_source[0])
            if tool in _DOCKING_SCORE_PROPERTIES
        ]
        raw = _read_docking_score_metrics(
            target_run / "stages" / "05_docking" / "docking_out.sdf",
            selected_tools,
            _configured_docking_score_properties(selected_tools),
        )
        if raw is None:
            raise ValueError("Docking target metrics are missing.")
        metrics = _metrics_by_target_id(raw, target_ids, "Docking")
        specs = {tool: {"max"} for tool in selected_tools}
        penalties["docking"] = _metric_extremeness(metrics, specs)
        eligible &= metrics[selected_tools].notna().all(axis=1)
        stage_metrics["docking"] = metrics

    if not penalties:
        raise ValueError("No target metrics are available for global alignment.")
    if int(eligible.sum()) < required_count:
        raise ValueError(
            "Global target alignment cannot protect "
            f"{required_count}/{len(target_ids)} molecules: only "
            f"{int(eligible.sum())} have all required stage measurements."
        )

    penalty_frame = pd.DataFrame(penalties, index=target_ids).fillna(1.0)
    ranking = pd.DataFrame(
        {
            "eligible": eligible,
            "worst": penalty_frame.max(axis=1),
            "mean": penalty_frame.mean(axis=1),
            "source_order": range(len(target_ids)),
        },
        index=target_ids,
    )
    selected_ids = ranking.sort_values(
        ["eligible", "worst", "mean", "source_order"],
        ascending=[False, True, True, True],
        kind="stable",
    ).index[:required_count]
    if not bool(eligible.loc[selected_ids].all()):
        raise ValueError(
            "Global target cohort contains molecules with missing metrics."
        )

    protected = target_molecules.copy()
    protected["mol_idx"] = protected["mol_idx"].astype(str)
    protected = protected.set_index("mol_idx", drop=False).loc[selected_ids].copy()
    protected["alignment_worst_extremeness"] = ranking.loc[selected_ids, "worst"]
    protected["alignment_mean_extremeness"] = ranking.loc[selected_ids, "mean"]
    return protected.reset_index(drop=True), stage_metrics


def _protected_numeric_pass_count(
    metrics: pd.DataFrame,
    specs: dict[str, set[str]],
    config: dict[str, Any],
) -> int:
    passed = pd.Series(True, index=metrics.index)
    borders = config.get("borders", {})
    for column, sides in specs.items():
        if column == "ring_size":
            values_by_side = {
                "min": _ring_size_extrema(metrics, "min"),
                "max": _ring_size_extrema(metrics, "max"),
            }
        else:
            values = pd.to_numeric(metrics.get(column), errors="coerce")
            values_by_side = {"min": values, "max": values}
        for side in sides:
            threshold = borders.get(f"{column}_{side}")
            if threshold is None:
                continue
            values = values_by_side[side]
            if values is None:
                passed &= False
                continue
            comparison = (
                values >= float(threshold)
                if side == "min"
                else values <= float(threshold)
            )
            passed &= values.isna() | comparison
    return int(passed.sum())


def finalize_global_alignment(
    master: dict[str, Any],
    target_run: Path,
    alignment_root: Path,
    target_mols_path: str,
    percentile: float,
) -> tuple[dict[str, Any], Path, Path]:
    """Build production configs around one globally protected target cohort."""
    percentile = validate_target_coverage_percent(percentile)
    source_master_path = alignment_root / SOURCE_CONFIGS_DIR_NAME / "source_config.yml"
    source_master = (
        load_config(str(source_master_path))
        if source_master_path.is_file()
        else copy.deepcopy(master)
    )
    for runtime_key in _RUNTIME_STAGE_OVERRIDE_KEYS:
        if runtime_key in master:
            source_master[runtime_key] = copy.deepcopy(master[runtime_key])
    descriptor_bounds_mode = descriptor_bounds_mode_from_master(source_master)
    _drop_synthesis_bounds_mode(source_master)
    target_molecules = _read_csv(target_run / "input" / "sampled_molecules.csv")
    if target_molecules is None or target_molecules.empty:
        raise ValueError("Saved target molecules are missing for global alignment.")
    if "mol_idx" not in target_molecules:
        raise ValueError("Saved target molecules do not contain stable mol_idx values.")
    target_molecules["mol_idx"] = target_molecules["mol_idx"].astype(str)
    if target_molecules["mol_idx"].duplicated().any():
        raise ValueError("Saved target molecules contain duplicate mol_idx values.")

    protected, stage_metrics = _select_global_protected_cohort(
        source_master, target_run, target_molecules, percentile
    )
    protected_ids = pd.Index(protected["mol_idx"].astype(str), name="mol_idx")
    required_count = len(protected)
    aligned_dir = alignment_root / "aligned_configs"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    protected_path = aligned_dir / "protected_target_molecules.csv"
    protected.to_csv(protected_path, index=False)
    aligned = copy.deepcopy(source_master)
    for runtime_key in (
        "_continue_mode",
        "_continue_completed_stages",
        "_run_single_stage_override",
        "_run_stage_selection_override",
    ):
        aligned.pop(runtime_key, None)
    summary = _new_threshold_summary(target_run, target_mols_path, percentile)
    summary["selection_method"] = "global_protected_target_cohort"
    summary["descriptor_bounds_mode"] = descriptor_bounds_mode
    verified_stages: list[str] = []

    def write_stage(
        stage: str,
        config_key: str,
        updater,
        metrics: pd.DataFrame,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        source = _config_from_master(source_master, config_key)
        if source is None:
            raise ValueError(f"Missing source config for aligned stage {stage}.")
        config, source_path = source
        target_path = aligned_dir / f"{config_key}{source_path.suffix or '.yml'}"
        protected_metrics = metrics.reindex(protected_ids).copy()
        thresholds = updater(config, protected_metrics, 100.0)
        config = _dump_aligned_stage_yaml(
            config,
            target_path,
            source=source_path,
            config_key=config_key,
        )
        aligned[config_key] = str(target_path.resolve())
        summary["stages"][stage] = {
            "thresholds": thresholds,
            "status": "globally_aligned",
            "protected_molecules": required_count,
        }
        verified_stages.append(stage)
        return config, thresholds

    molprep_source = _config_from_master(source_master, _CONFIG_MOL_PREP)
    if molprep_source is not None:
        molprep_config, source_path = molprep_source
        allowed_atoms = _target_atom_symbols(target_molecules)
        filters = molprep_config.setdefault("filters", {})
        filters["allowed_atoms"] = allowed_atoms
        target_path = aligned_dir / f"{_CONFIG_MOL_PREP}{source_path.suffix or '.yml'}"
        molprep_config = _dump_aligned_stage_yaml(
            molprep_config,
            target_path,
            source=source_path,
            config_key=_CONFIG_MOL_PREP,
        )
        aligned[_CONFIG_MOL_PREP] = str(target_path.resolve())
        summary["stages"]["mol_prep"] = {
            "thresholds": {"filters.allowed_atoms": allowed_atoms},
            "status": "globally_aligned",
            "protected_molecules": required_count,
        }
        verified_stages.append("mol_prep")

    if "descriptors" in stage_metrics:
        config, thresholds = write_stage(
            "descriptors",
            _CONFIG_DESCRIPTORS,
            lambda config, metrics, value: _align_descriptor_config(
                config,
                metrics,
                value,
                bounds_mode=descriptor_bounds_mode,
            ),
            stage_metrics["descriptors"],
        )
        summary["stages"]["descriptors"]["bounds_mode"] = descriptor_bounds_mode
        protected_metrics = stage_metrics["descriptors"].reindex(protected_ids)
        retained = _protected_numeric_pass_count(
            protected_metrics, _descriptor_threshold_specs(config), config
        )
        summary["stages"]["descriptors"]["protected_retained_molecules"] = retained
        if retained < required_count:
            raise ValueError(
                f"Descriptor alignment protects only {retained}/{required_count} molecules."
            )

    if "struct_filters" in stage_metrics:
        source = _config_from_master(source_master, _CONFIG_STRUCT_FILTERS)
        if source is None:
            raise ValueError("Missing source config for structural filters.")
        config, source_path = source
        target_path = aligned_dir / (
            f"{_CONFIG_STRUCT_FILTERS}{source_path.suffix or '.yml'}"
        )
        shutil.copyfile(source_path, target_path)
        aligned[_CONFIG_STRUCT_FILTERS] = str(target_path.resolve())
        thresholds = _align_structural_filter_config(
            config,
            stage_metrics["struct_filters"],
            percentile,
        )
        failure_audit_path = aligned_dir / "structural_filter_failures.csv"
        _write_structural_failure_audit(
            stage_metrics["struct_filters"], failure_audit_path
        )
        thresholds["failure_audit_path"] = str(failure_audit_path.resolve())
        summary["stages"]["struct_filters"] = {
            "thresholds": thresholds,
            "status": "source_config_preserved",
            "note": (
                "Structural hard filters and their parameters are copied unchanged "
                "and are not calibrated by target_coverage_percent."
            ),
        }

    if "docking" in stage_metrics:
        _config, thresholds = write_stage(
            "docking",
            _CONFIG_DOCKING,
            _align_docking_config,
            stage_metrics["docking"],
        )
        retained = int(thresholds.get("retained_molecules", 0))
        summary["stages"]["docking"]["protected_retained_molecules"] = retained
        if retained < required_count:
            raise ValueError(
                f"Docking alignment protects only {retained}/{required_count} molecules."
            )

    docking_filters_source = _config_from_master(source_master, _CONFIG_DOCKING_FILTERS)
    if docking_filters_source is not None:
        config, source_path = docking_filters_source
        target_path = (
            aligned_dir / f"{_CONFIG_DOCKING_FILTERS}{source_path.suffix or '.yml'}"
        )
        shutil.copyfile(source_path, target_path)
        aligned[_CONFIG_DOCKING_FILTERS] = str(target_path.resolve())
        summary["stages"]["docking_filters"] = {
            "thresholds": {},
            "status": "source_config_preserved",
            "note": (
                "Copied unchanged because pose-filter thresholds were not measured "
                "during target calibration."
            ),
        }

    summary["global_guarantee"] = {
        "status": "verified",
        "target_molecules": len(target_molecules),
        "required_retained_molecules": required_count,
        "guaranteed_retained_molecules": required_count,
        "guaranteed_retained_percent": float(
            required_count / len(target_molecules) * 100.0
        ),
        "protected_cohort_path": str(protected_path.resolve()),
        "verified_stages": verified_stages,
    }
    thresholds_path = aligned_dir / THRESHOLDS_NAME
    master_path = aligned_dir / ALIGNED_CONFIG_NAME
    if "target_mols_path" in aligned:
        aligned["target_mols_path"] = str(Path(target_mols_path).resolve())
    alignment = aligned.get("alignment")
    if isinstance(alignment, dict):
        if "enabled" in alignment:
            alignment["enabled"] = False
        if "target_coverage_percent" in alignment:
            alignment["target_coverage_percent"] = percentile
        alignment.pop("synthesis_bounds_mode", None)
    _dump_yaml(summary, thresholds_path)
    source_master_template = (
        source_master_path if source_master_path.is_file() else None
    )
    aligned = _write_generated_master(
        aligned,
        master_path,
        source=source_master_template,
    )
    logger.info(
        "Verified global target coverage: %d/%d molecules (%.2f%%) are protected.",
        required_count,
        len(target_molecules),
        required_count / len(target_molecules) * 100.0,
    )
    return aligned, master_path, thresholds_path


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
    descriptor_bounds_mode = descriptor_bounds_mode_from_master(master)
    _drop_synthesis_bounds_mode(master)
    aligned_dir = alignment_root / "aligned_configs"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    master_path = aligned_dir / ALIGNED_CONFIG_NAME
    thresholds_path = aligned_dir / THRESHOLDS_NAME

    aligned = (
        load_config(str(master_path)) if master_path.exists() else copy.deepcopy(master)
    )
    _drop_synthesis_bounds_mode(aligned)
    summary = (
        load_config(str(thresholds_path))
        if thresholds_path.exists()
        else _new_threshold_summary(target_run, target_mols_path, percentile)
    )
    summary.pop("retention_percentile", None)
    summary["target_coverage_percent"] = percentile
    summary["descriptor_bounds_mode"] = descriptor_bounds_mode
    summary.pop("synthesis_bounds_mode", None)

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
        specs = [
            (
                _CONFIG_MOL_PREP,
                lambda config, metrics, _coverage: _align_molprep_config(
                    config, metrics
                ),
            )
        ]
    elif stage == "descriptors":
        metrics = _read_csv(
            target_run
            / "stages"
            / "02_descriptors_initial"
            / "metrics"
            / "descriptors_all.csv"
        )
        specs = [
            (
                _CONFIG_DESCRIPTORS,
                lambda config, values, coverage: _align_descriptor_config(
                    config,
                    values,
                    coverage,
                    bounds_mode=descriptor_bounds_mode,
                ),
            )
        ]
    elif stage == "struct_filters":
        metrics = _read_structural_rule_masks(
            target_run / "stages" / "03_structural_filters_post"
        )
        specs = [(_CONFIG_STRUCT_FILTERS, _align_structural_filter_config)]
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
        selected_tools = _parse_tools_config(docking_config)
        metrics = _read_docking_score_metrics(
            target_run / "stages" / "05_docking" / "docking_out.sdf",
            selected_tools,
            _configured_docking_score_properties(selected_tools),
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
        config = load_config(str(source))
        thresholds = updater(config, metrics, percentile)
        if stage == "struct_filters":
            failure_audit_path = aligned_dir / "structural_filter_failures.csv"
            _write_structural_failure_audit(metrics, failure_audit_path)
            thresholds["failure_audit_path"] = str(failure_audit_path.resolve())
        config = _dump_aligned_stage_yaml(
            config,
            target,
            source=source,
            config_key=config_key,
        )
        aligned[config_key] = str(target.resolve())
        summary["stages"][stage]["thresholds"] = thresholds
        if stage == "struct_filters":
            summary["stages"][stage]["status"] = "source_config_preserved"
            summary["stages"][stage]["note"] = (
                "Structural hard filters and their parameters are copied unchanged "
                "and are not calibrated by target_coverage_percent."
            )
        else:
            summary["stages"][stage]["status"] = "ready"
        if stage == "descriptors":
            summary["stages"][stage]["bounds_mode"] = descriptor_bounds_mode

    if "target_mols_path" in aligned:
        aligned["target_mols_path"] = str(Path(target_mols_path).resolve())
    alignment = aligned.get("alignment")
    if isinstance(alignment, dict):
        if "enabled" in alignment:
            alignment["enabled"] = False
        if "target_coverage_percent" in alignment:
            alignment["target_coverage_percent"] = percentile
        alignment.pop("synthesis_bounds_mode", None)
    _dump_yaml(summary, thresholds_path)
    master_template = master_path if master_path.is_file() else None
    if master_template is None:
        candidate = alignment_root / SOURCE_CONFIGS_DIR_NAME / "source_config.yml"
        master_template = candidate if candidate.is_file() else None
    aligned = _write_generated_master(
        aligned,
        master_path,
        source=master_template,
    )
    return aligned, master_path, thresholds_path
