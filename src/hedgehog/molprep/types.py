from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MolPrepFailure:
    smiles_raw: str
    model_name: str | None
    mol_idx: str | None
    reason: str
    step: str
    reason_detail: str | None = None
