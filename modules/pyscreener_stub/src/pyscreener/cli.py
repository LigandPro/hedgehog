from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path
from typing import Iterable


def _parse_ini(lines: Iterable[str]) -> dict[str, str]:
    data: dict[str, str] = {}
    for raw in lines:
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip().lower()] = value.strip()
    return data


def main() -> None:
    """Minimal stub for the `pyscreener` CLI."""
    parser = argparse.ArgumentParser(description="Stub pyscreener CLI")
    parser.add_argument("--config", required=True, help="Path to SMINA INI file")
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Seed for deterministic dummy scores",
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = Path.cwd() / config_path

    try:
        lines = config_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:  # pragma: no cover
        raise SystemExit(f"Failed to read config file: {exc}") from exc

    config_data = _parse_ini(lines)
    output_dir_name = config_data.get("output-dir", "smina_results")
    output_dir = Path(output_dir_name)
    if not output_dir.is_absolute():
        output_dir = Path.cwd() / output_dir

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "scores.csv"

    random.seed(args.seed)
    rows = [[f"ligand_{idx+1}", f"{random.uniform(-10, -5):.2f}"] for idx in range(5)]

    try:
        with output_file.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["ligand", "score"])
            writer.writerows(rows)
    except OSError as exc:  # pragma: no cover
        raise SystemExit(f"Failed to write stub docking results: {exc}") from exc

    print(f"Stub docking results written to {output_file}")


if __name__ == "__main__":  # pragma: no cover
    main()
