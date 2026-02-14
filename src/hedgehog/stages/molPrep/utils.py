from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import datamol as dm
import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import logger


@dataclass(frozen=True)
class MolPrepFailure:
    smiles_raw: str
    model_name: str | None
    mol_idx: str | None
    reason: str
    step: str
    reason_detail: str | None = None


def _get_cfg(cfg: dict[str, Any], path: list[str], default: Any) -> Any:
    cur: Any = cfg
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _is_single_fragment(mol: Chem.Mol) -> bool:
    try:
        frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
        return len(frags) <= 1
    except Exception:
        return False


def _has_radicals(mol: Chem.Mol) -> bool:
    return any(a.GetNumRadicalElectrons() > 0 for a in mol.GetAtoms())


def _is_neutral(mol: Chem.Mol) -> bool:
    return not any(a.GetFormalCharge() != 0 for a in mol.GetAtoms())


def _has_isotopes(mol: Chem.Mol) -> bool:
    return any(a.GetIsotope() != 0 for a in mol.GetAtoms())


def _allowed_atoms_ok(mol: Chem.Mol, allowed: set[str]) -> bool:
    if not allowed:
        return True
    return all(a.GetSymbol() in allowed for a in mol.GetAtoms())


def _safe_to_mol(
    smiles: str,
    to_mol_cfg: dict[str, Any],
) -> Chem.Mol | None:
    try:
        return dm.mol.to_mol(
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
) -> tuple[str | None, str | None, str | None]:
    """Return (smiles, reason, step) where smiles is standardized or None on failure."""
    to_mol_cfg = _get_cfg(cfg, ["steps", "to_mol"], {}) or {}
    mol = _safe_to_mol(smiles_raw, to_mol_cfg)
    if mol is None:
        return None, "parse_failed", "to_mol"

    if _get_cfg(cfg, ["steps", "fix_mol", "enabled"], True):
        try:
            fix_cfg = _get_cfg(cfg, ["steps", "fix_mol"], {}) or {}
            mol = dm.mol.fix_mol(
                mol,
                n_iter=int(fix_cfg.get("n_iter", 1)),
                remove_singleton=bool(fix_cfg.get("remove_singleton", True)),
                largest_only=bool(fix_cfg.get("largest_only", False)),
                inplace=False,
            )
        except Exception:
            return None, "fix_failed", "fix_mol"

    if _get_cfg(cfg, ["steps", "sanitize_mol", "enabled"], True):
        try:
            mol = dm.mol.sanitize_mol(mol)
        except Exception:
            mol = None
        if mol is None:
            return None, "sanitize_failed", "sanitize_mol"

    if _get_cfg(cfg, ["steps", "remove_salts_solvents", "enabled"], True):
        try:
            rss_cfg = _get_cfg(cfg, ["steps", "remove_salts_solvents"], {}) or {}
            mol = dm.mol.remove_salts_solvents(
                mol,
                defn_data=rss_cfg.get("defn_data"),
                defn_format=str(rss_cfg.get("defn_format", "smarts")),
                dont_remove_everything=bool(
                    rss_cfg.get("dont_remove_everything", True)
                ),
                sanitize=bool(rss_cfg.get("sanitize", True)),
            )
        except Exception:
            return None, "remove_salts_failed", "remove_salts_solvents"
        if mol is None:
            return None, "remove_salts_failed", "remove_salts_solvents"

    if bool(_get_cfg(cfg, ["steps", "keep_largest_fragment"], True)):
        try:
            mol = dm.mol.keep_largest_fragment(mol)
        except Exception:
            return None, "largest_fragment_failed", "keep_largest_fragment"

        if mol is None:
            return None, "largest_fragment_failed", "keep_largest_fragment"

        require_single_fragment = bool(
            _get_cfg(cfg, ["filters", "require_single_fragment"], True)
        )
        if require_single_fragment and not _is_single_fragment(mol):
            return None, "multifragment_after_largest", "keep_largest_fragment"

    if _get_cfg(cfg, ["steps", "standardize_mol", "enabled"], True):
        try:
            std_cfg = _get_cfg(cfg, ["steps", "standardize_mol"], {}) or {}
            mol = dm.mol.standardize_mol(
                mol,
                disconnect_metals=bool(std_cfg.get("disconnect_metals", True)),
                normalize=bool(std_cfg.get("normalize", True)),
                reionize=bool(std_cfg.get("reionize", True)),
                uncharge=bool(std_cfg.get("uncharge", True)),
                stereo=bool(std_cfg.get("stereo", True)),
            )
        except Exception:
            return None, "standardize_mol_failed", "standardize_mol"
        if mol is None:
            return None, "standardize_mol_failed", "standardize_mol"

    if bool(_get_cfg(cfg, ["steps", "remove_stereochemistry"], True)):
        try:
            Chem.RemoveStereochemistry(mol)
        except Exception:
            return None, "remove_stereo_failed", "remove_stereochemistry"

    try:
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None, "to_smiles_failed", "to_smiles"

    if _get_cfg(cfg, ["steps", "standardize_smiles", "enabled"], True):
        try:
            smiles = dm.mol.standardize_smiles(smiles)
        except Exception:
            return None, "standardize_smiles_failed", "standardize_smiles"

        mol2 = dm.to_mol(smiles, sanitize=True)
        if mol2 is None:
            return None, "post_standardize_parse_failed", "post_standardize_parse"
        mol = mol2

    # Final strict filters
    allowed_atoms = set(_get_cfg(cfg, ["filters", "allowed_atoms"], []) or [])
    if not _allowed_atoms_ok(mol, allowed_atoms):
        return None, "disallowed_atoms", "filters"

    if bool(_get_cfg(cfg, ["filters", "reject_radicals"], True)) and _has_radicals(mol):
        return None, "radicals", "filters"

    if bool(_get_cfg(cfg, ["filters", "reject_isotopes"], True)) and _has_isotopes(mol):
        return None, "isotopes", "filters"

    if bool(_get_cfg(cfg, ["filters", "require_neutral"], True)) and not _is_neutral(
        mol
    ):
        return None, "charged", "filters"

    if bool(_get_cfg(cfg, ["filters", "require_single_fragment"], True)) and (
        not _is_single_fragment(mol)
    ):
        return None, "multifragment", "filters"

    # Normalize final SMILES once more after filters (canonical, no stereo)
    try:
        smiles_final = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None, "to_smiles_failed", "to_smiles_final"

    if "." in smiles_final and bool(
        _get_cfg(cfg, ["filters", "require_single_fragment"], True)
    ):
        return None, "multifragment", "filters"

    return smiles_final, None, None


def run_mol_prep(
    df: pd.DataFrame,
    cfg: dict[str, Any],
    out_dir: Path,
) -> pd.DataFrame:
    """Run MolPrep on a dataframe and write stage outputs."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cols_cfg = _get_cfg(cfg, ["columns"], {}) or {}
    smiles_col = str(cols_cfg.get("smiles", "smiles"))
    model_col = str(cols_cfg.get("model_name", "model_name"))
    idx_col = str(cols_cfg.get("mol_idx", "mol_idx"))
    smiles_raw_col = str(cols_cfg.get("smiles_raw", "smiles_raw"))

    passed_rows: list[dict[str, Any]] = []
    failures: list[MolPrepFailure] = []

    for _, row in df.iterrows():
        smiles_raw = str(row.get(smiles_col, "")).strip()
        model_name = row.get(model_col)
        mol_idx = row.get(idx_col)

        smiles_std, reason, step = _molprep_one(smiles_raw, cfg)
        if smiles_std is None:
            failures.append(
                MolPrepFailure(
                    smiles_raw=smiles_raw,
                    model_name=str(model_name) if model_name is not None else None,
                    mol_idx=str(mol_idx) if mol_idx is not None else None,
                    reason=str(reason or "unknown"),
                    step=str(step or "unknown"),
                    reason_detail=None,
                )
            )
            continue

        passed_rows.append(
            {
                "smiles": smiles_std,
                smiles_raw_col: smiles_raw,
                "model_name": model_name,
                "mol_idx": mol_idx,
            }
        )

    passed_df = pd.DataFrame(passed_rows)
    failed_df = pd.DataFrame(
        [
            {
                "smiles": f.smiles_raw,  # Keep compatibility: 'smiles' column present
                smiles_raw_col: f.smiles_raw,
                "model_name": f.model_name,
                "mol_idx": f.mol_idx,
                "reason": f.reason,
                "step": f.step,
                "reason_detail": f.reason_detail,
            }
            for f in failures
        ]
    )

    # Deduplicate after standardization (within model_name)
    dedup_removed = 0
    duplicates_df = pd.DataFrame()
    if not passed_df.empty:
        before = len(passed_df)
        dup_mask = passed_df.duplicated(subset=["smiles", "model_name"], keep="first")
        duplicates_df = passed_df[dup_mask].copy() if dup_mask.any() else pd.DataFrame()
        passed_df = passed_df[~dup_mask].reset_index(drop=True)
        dedup_removed = before - len(passed_df)

    # Write outputs
    filtered_path = out_dir / "filtered_molecules.csv"
    failed_path = out_dir / "failed_molecules.csv"
    metrics_path = out_dir / "metrics.csv"

    stable_cols = ["smiles", smiles_raw_col, "model_name", "mol_idx"]
    for c in stable_cols:
        if c not in passed_df.columns:
            passed_df[c] = pd.NA
        if c not in failed_df.columns:
            failed_df[c] = pd.NA

    passed_df[stable_cols].to_csv(filtered_path, index=False)
    if failed_df.empty:
        failed_df = pd.DataFrame(
            columns=[*stable_cols, "reason", "step", "reason_detail"]
        )
    failed_df.to_csv(failed_path, index=False)

    reasons = Counter([f.reason for f in failures])
    metrics = {
        "total_in": int(len(df)),
        "passed": int(len(passed_df)),
        "failed": int(len(failures)),
        "dedup_removed": int(dedup_removed),
        "failed_by_reason_json": json.dumps(
            dict(sorted(reasons.items())), ensure_ascii=False
        ),
    }
    pd.DataFrame([metrics]).to_csv(metrics_path, index=False)

    if (
        bool(_get_cfg(cfg, ["output", "write_duplicates_removed"], True))
        and not duplicates_df.empty
    ):
        duplicates_df.to_csv(out_dir / "duplicates_removed.csv", index=False)

    logger.info(
        "MolPrep: %d in, %d passed, %d failed, %d dedup removed",
        len(df),
        len(passed_df),
        len(failures),
        dedup_removed,
    )

    # Also log top failure reasons (human-friendly)
    if reasons:
        top = reasons.most_common(10)
        logger.info("MolPrep failure reasons (top): %s", json.dumps(top))

    return passed_df[stable_cols]
