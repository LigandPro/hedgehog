"""Isolated Shepherd-Score worker process."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from rdkit import Chem


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute Shepherd-Score shape similarity in an isolated process"
    )
    parser.add_argument("--input-sdf", required=True)
    parser.add_argument("--reference-sdf", required=True)
    parser.add_argument("--params-json", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--meta-json", required=True)
    return parser.parse_args()


def _load_first_mol(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise RuntimeError(f"No valid molecule found in {path}")
    return mol


def _shape_score(
    mol: Chem.Mol,
    ref_positions_list: list[list[float]],
    alpha: float,
) -> float:
    import torch
    from shepherd_score.score.gaussian_overlap import shape_tanimoto

    if mol.GetNumConformers() == 0:
        raise RuntimeError("Input molecule has no conformer")

    ref_positions = torch.tensor(ref_positions_list, dtype=torch.float32)
    mol_h = Chem.AddHs(mol, addCoords=True)
    mol_conf = mol_h.GetConformer()
    mol_positions = torch.tensor(mol_conf.GetPositions(), dtype=torch.float32)
    return float(shape_tanimoto(ref_positions, mol_positions, alpha=alpha).item())


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def main() -> int:
    args = _parse_args()
    started = time.time()

    input_sdf = Path(args.input_sdf)
    reference_sdf = Path(args.reference_sdf)
    params_json = Path(args.params_json)
    output_json = Path(args.output_json)
    meta_json = Path(args.meta_json)

    try:
        with open(params_json, "r", encoding="utf-8") as f:
            params = json.load(f)
        alpha = float(params.get("alpha", 0.81))
        min_shape_score = float(params.get("min_shape_score", 0.5))

        reference_mol = _load_first_mol(reference_sdf)
        if reference_mol.GetNumConformers() == 0:
            raise RuntimeError("Reference molecule has no conformer")

        ref_h = Chem.AddHs(reference_mol, addCoords=True)
        ref_conf = ref_h.GetConformer()
        ref_positions_list = ref_conf.GetPositions().tolist()

        supplier = Chem.SDMolSupplier(str(input_sdf), sanitize=False, removeHs=False)
        results: list[dict[str, object]] = []

        for idx, mol in enumerate(supplier):
            if mol is None:
                results.append(
                    {
                        "mol_idx": idx,
                        "shape_score": 0.0,
                        "pass_shepherd_score": False,
                        "error": "invalid molecule",
                    }
                )
                continue

            try:
                score = _shape_score(mol, ref_positions_list, alpha)
                results.append(
                    {
                        "mol_idx": idx,
                        "shape_score": score,
                        "pass_shepherd_score": score >= min_shape_score,
                        "error": None,
                    }
                )
            except Exception as exc:  # noqa: BLE001
                results.append(
                    {
                        "mol_idx": idx,
                        "shape_score": 0.0,
                        "pass_shepherd_score": False,
                        "error": str(exc),
                    }
                )

        _write_json(output_json, results)
        _write_json(
            meta_json,
            {
                "ok": True,
                "count": len(results),
                "elapsed_sec": round(time.time() - started, 3),
            },
        )
        return 0
    except Exception as exc:  # noqa: BLE001
        _write_json(
            meta_json,
            {
                "ok": False,
                "error": str(exc),
                "elapsed_sec": round(time.time() - started, 3),
            },
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
