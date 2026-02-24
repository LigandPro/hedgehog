"""Isolated legacy RAScore worker process.

This worker exists to evaluate the MolScore-provided ``model.pkl`` with an
older xgboost runtime in a subprocess, while the main project environment can
keep a newer xgboost version.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem import AllChem


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Legacy RAScore batch predictor")
    parser.add_argument("--model-pkl", required=True)
    parser.add_argument("--input-json", required=True)
    parser.add_argument("--output-json", required=True)
    return parser.parse_args()


def _ecfp6_counts(smiles: str, n_bits: int = 2048) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprint(mol, 3, useFeatures=False, useCounts=True)
    arr = np.zeros((n_bits,), dtype=np.int32)
    for idx, value in fp.GetNonzeroElements().items():
        arr[idx % n_bits] += value
    return arr


def _write_output(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f)


def main() -> int:
    args = _parse_args()
    model_path = Path(args.model_pkl)
    input_path = Path(args.input_json)
    output_path = Path(args.output_json)

    try:
        payload = json.loads(input_path.read_text(encoding="utf-8"))
        smiles_list = payload.get("smiles", [])
        with open(model_path, "rb") as f:
            model = pickle.load(f)

        fps = []
        valid_indices = []
        for idx, smiles in enumerate(smiles_list):
            fp = _ecfp6_counts(smiles)
            if fp is not None:
                valid_indices.append(idx)
                fps.append(fp)

        scores: list[float | None] = [None] * len(smiles_list)
        if fps:
            booster = model.get_booster() if hasattr(model, "get_booster") else model
            probs = booster.predict(xgb.DMatrix(np.array(fps)))
            for idx, prob in zip(valid_indices, probs):
                scores[idx] = float(prob)

        _write_output(output_path, {"scores": scores})
        return 0
    except Exception as exc:  # noqa: BLE001
        _write_output(output_path, {"scores": [], "error": str(exc)})
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
