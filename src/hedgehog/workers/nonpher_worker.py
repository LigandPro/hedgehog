"""Run Nonpher complexity scoring in an external Python environment."""

from __future__ import annotations

import argparse
import csv
from collections.abc import Sequence
from typing import Any


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nonpher isolated worker launcher")
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--smiles-column", default="smiles")
    parser.add_argument("--score-column", default="nonpher_complexity_score")
    return parser.parse_args(argv)


def _create_complexity_filter() -> Any:
    from nonpher import nonpher as nonpher_mod

    return nonpher_mod.ComplexityFilter()


def _mol_from_smiles(smiles: str) -> Any:
    from rdkit import Chem

    return Chem.MolFromSmiles(smiles)


def _score_smiles(smiles: str, complexity_filter: Any) -> str:
    if not smiles:
        return ""

    mol = _mol_from_smiles(smiles)
    if mol is None:
        return ""

    try:
        return "1.0" if complexity_filter.isTooComplex(mol) else "0.0"
    except Exception:
        return ""


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)

    try:
        with open(args.input_csv, encoding="utf-8", newline="") as source:
            reader = csv.DictReader(source)
            fieldnames = list(reader.fieldnames or [])
            if args.smiles_column not in fieldnames:
                raise ValueError(
                    f"Missing required smiles column: {args.smiles_column}"
                )
            if args.score_column not in fieldnames:
                fieldnames.append(args.score_column)

            complexity_filter = _create_complexity_filter()
            rows: list[dict[str, str]] = []
            for row in reader:
                row[args.score_column] = _score_smiles(
                    row.get(args.smiles_column, ""),
                    complexity_filter,
                )
                rows.append(row)

        with open(args.output_csv, "w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
