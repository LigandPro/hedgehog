from __future__ import annotations

import json
import os
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog.molprep.filters import (
    _allowed_atoms_ok,
    _has_isotopes,
    _has_radicals,
    _is_single_fragment,
)
from hedgehog.molprep.orchestrator import _get_cfg
from hedgehog.molprep.types import MolPrepFailure
from hedgehog.utils.datamol_import import import_datamol_quietly

_DATAMOL_WORKER_WARMED = False
_MOLPREP_CFG: dict[str, Any] | None = None
dm = import_datamol_quietly()


def _is_missing_smiles_value(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and pd.isna(value):
        return True
    if not isinstance(value, str):
        return False
    token = value.strip().lower()
    return token in {"", "nan", "none"}


def _safe_to_mol(
    smiles: str,
    to_mol_cfg: dict[str, Any],
) -> Chem.Mol | None:
    try:
        return dm.to_mol(
            smiles,
            ordered=bool(to_mol_cfg.get("ordered", True)),
            sanitize=bool(to_mol_cfg.get("sanitize", False)),
            allow_cxsmiles=bool(to_mol_cfg.get("allow_cxsmiles", True)),
            strict_cxsmiles=bool(to_mol_cfg.get("strict_cxsmiles", True)),
            remove_hs=bool(to_mol_cfg.get("remove_hs", True)),
        )
    except Exception:
        return None


def _molprep_one(
    smiles_raw: str,
    cfg: dict[str, Any],
) -> tuple[str | None, str | None, str | None, str | None]:
    """Return (smiles, reason, step, reason_detail_json) on success/failure."""
    if _is_missing_smiles_value(smiles_raw):
        return None, "parse_failed", "to_mol", None

    smiles_raw = str(smiles_raw).strip()
    to_mol_cfg = _get_cfg(cfg, ["steps", "to_mol"], {}) or {}
    mol = _safe_to_mol(smiles_raw, to_mol_cfg)
    if mol is None:
        return None, "parse_failed", "to_mol", None

    if _get_cfg(cfg, ["steps", "fix_mol", "enabled"], True):
        try:
            fix_cfg = _get_cfg(cfg, ["steps", "fix_mol"], {}) or {}
            mol = dm.fix_mol(
                mol,
                n_iter=int(fix_cfg.get("n_iter", 1)),
                remove_singleton=bool(fix_cfg.get("remove_singleton", True)),
                largest_only=bool(fix_cfg.get("largest_only", False)),
                inplace=False,
            )
        except Exception:
            return None, "fix_failed", "fix_mol", None

    if _get_cfg(cfg, ["steps", "sanitize_mol", "enabled"], True):
        try:
            mol = dm.sanitize_mol(mol)
        except Exception:
            mol = None
        if mol is None:
            return None, "sanitize_failed", "sanitize_mol", None

    if _get_cfg(cfg, ["steps", "remove_salts_solvents", "enabled"], True):
        try:
            rss_cfg = _get_cfg(cfg, ["steps", "remove_salts_solvents"], {}) or {}
            mol = dm.remove_salts_solvents(
                mol,
                defn_data=rss_cfg.get("defn_data"),
                defn_format=str(rss_cfg.get("defn_format", "smarts")),
                dont_remove_everything=bool(
                    rss_cfg.get("dont_remove_everything", True)
                ),
                sanitize=bool(rss_cfg.get("sanitize", True)),
            )
        except Exception:
            return None, "remove_salts_failed", "remove_salts_solvents", None
        if mol is None:
            return None, "remove_salts_failed", "remove_salts_solvents", None

    if bool(_get_cfg(cfg, ["steps", "keep_largest_fragment"], True)):
        try:
            mol = dm.keep_largest_fragment(mol)
        except Exception:
            return None, "largest_fragment_failed", "keep_largest_fragment", None

        if mol is None:
            return None, "largest_fragment_failed", "keep_largest_fragment", None

        require_single_fragment = bool(
            _get_cfg(cfg, ["filters", "require_single_fragment"], True)
        )
        if require_single_fragment and not _is_single_fragment(mol):
            return None, "multifragment_after_largest", "keep_largest_fragment", None

    if _get_cfg(cfg, ["steps", "standardize_mol", "enabled"], True):
        try:
            std_cfg = _get_cfg(cfg, ["steps", "standardize_mol"], {}) or {}
            mol = dm.standardize_mol(
                mol,
                disconnect_metals=bool(std_cfg.get("disconnect_metals", True)),
                normalize=bool(std_cfg.get("normalize", True)),
                reionize=bool(std_cfg.get("reionize", True)),
                uncharge=bool(std_cfg.get("uncharge", True)),
                stereo=bool(std_cfg.get("stereo", True)),
            )
        except Exception:
            return None, "standardize_mol_failed", "standardize_mol", None
        if mol is None:
            return None, "standardize_mol_failed", "standardize_mol", None

    if bool(_get_cfg(cfg, ["steps", "remove_stereochemistry"], True)):
        try:
            Chem.RemoveStereochemistry(mol)
        except Exception:
            return None, "remove_stereo_failed", "remove_stereochemistry", None

    try:
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None, "to_smiles_failed", "to_smiles", None

    if _get_cfg(cfg, ["steps", "standardize_smiles", "enabled"], True):
        try:
            smiles = dm.standardize_smiles(smiles)
        except Exception:
            return None, "standardize_smiles_failed", "standardize_smiles", None

        if _is_missing_smiles_value(smiles):
            return None, "standardize_smiles_failed", "standardize_smiles", None
        smiles = str(smiles).strip()

        mol2 = dm.to_mol(smiles, sanitize=True)
        if mol2 is None:
            return None, "post_standardize_parse_failed", "post_standardize_parse", None
        mol = mol2

    # Final strict filters: evaluate independently so downstream pass rates are criterion-wise.
    allowed_atoms = set(_get_cfg(cfg, ["filters", "allowed_atoms"], []) or [])
    reject_radicals = bool(_get_cfg(cfg, ["filters", "reject_radicals"], True))
    reject_isotopes = bool(_get_cfg(cfg, ["filters", "reject_isotopes"], True))
    require_single_fragment = bool(
        _get_cfg(cfg, ["filters", "require_single_fragment"], True)
    )
    filter_flags = {
        "filter_allowed_atoms_pass": _allowed_atoms_ok(mol, allowed_atoms),
        "filter_radicals_pass": (not _has_radicals(mol)) if reject_radicals else True,
        "filter_isotopes_pass": (not _has_isotopes(mol)) if reject_isotopes else True,
        "filter_single_fragment_pass": _is_single_fragment(mol)
        if require_single_fragment
        else True,
    }
    if not filter_flags["filter_allowed_atoms_pass"]:
        return (
            None,
            "disallowed_atoms",
            "filters",
            json.dumps(filter_flags, ensure_ascii=False),
        )
    if not filter_flags["filter_radicals_pass"]:
        return None, "radicals", "filters", json.dumps(filter_flags, ensure_ascii=False)
    if not filter_flags["filter_isotopes_pass"]:
        return None, "isotopes", "filters", json.dumps(filter_flags, ensure_ascii=False)
    if not filter_flags["filter_single_fragment_pass"]:
        return (
            None,
            "multifragment",
            "filters",
            json.dumps(filter_flags, ensure_ascii=False),
        )

    # Normalize final SMILES once more after filters (canonical, no stereo)
    try:
        smiles_final = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None, "to_smiles_failed", "to_smiles_final", None

    if _is_missing_smiles_value(smiles_final):
        return None, "to_smiles_failed", "to_smiles_final", None
    smiles_final = str(smiles_final).strip()

    if "." in smiles_final and bool(
        _get_cfg(cfg, ["filters", "require_single_fragment"], True)
    ):
        return (
            None,
            "multifragment",
            "filters",
            json.dumps(filter_flags, ensure_ascii=False),
        )

    return smiles_final, None, None, json.dumps(filter_flags, ensure_ascii=False)


def _process_molprep_item(
    item: tuple[str, str | None, str | None],
) -> tuple[dict[str, Any] | None, MolPrepFailure | None]:
    global _DATAMOL_WORKER_WARMED
    if not _DATAMOL_WORKER_WARMED:
        stdout_fd = os.dup(1)
        stderr_fd = os.dup(2)
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        try:
            os.dup2(devnull_fd, 1)
            os.dup2(devnull_fd, 2)
            try:
                dm.to_mol("C")
            except Exception:
                pass
        finally:
            os.dup2(stdout_fd, 1)
            os.dup2(stderr_fd, 2)
            os.close(devnull_fd)
            os.close(stdout_fd)
            os.close(stderr_fd)
        _DATAMOL_WORKER_WARMED = True

    cfg = _MOLPREP_CFG
    if cfg is None:
        raise RuntimeError("MolPrep worker configuration was not initialized")

    smiles_raw, model_name, mol_idx = item
    smiles_std, reason, step, reason_detail = _molprep_one(smiles_raw, cfg)
    if smiles_std is None:
        return None, MolPrepFailure(
            smiles_raw=smiles_raw,
            model_name=model_name,
            mol_idx=mol_idx,
            reason=str(reason or "unknown"),
            step=str(step or "unknown"),
            reason_detail=reason_detail,
        )

    return {
        "smiles": smiles_std,
        "smiles_raw": smiles_raw,
        "model_name": model_name,
        "mol_idx": mol_idx,
    }, None


def _init_molprep_worker(cfg: dict[str, Any]) -> None:
    global _MOLPREP_CFG
    _MOLPREP_CFG = cfg
