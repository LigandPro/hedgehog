from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> None:
    """Minimal stub implementation of the `aizynthcli` entry point."""
    parser = argparse.ArgumentParser(description="Stub AiZynthFinder CLI")
    parser.add_argument("--config", required=True, help="Path to config YAML")
    parser.add_argument("--smiles", required=True, help="Input SMILES file")
    parser.add_argument("--output", required=True, help="Destination JSON file")
    parser.add_argument(
        "--max-routes",
        type=int,
        default=3,
        help="Retained for compatibility; ignored in the stub",
    )
    args = parser.parse_args()

    smiles_path = Path(args.smiles)
    output_path = Path(args.output)

    entries: list[dict[str, object]] = []
    try:
        with smiles_path.open() as handle:
            for idx, line in enumerate(handle):
                smi = line.strip()
                if not smi:
                    continue
                entries.append(
                    {
                        "index": idx,
                        "target": smi,
                        "is_solved": True,
                        "search_time": 0.1,
                    }
                )
    except OSError as exc:  # pragma: no cover
        raise SystemExit(f"Failed to read SMILES file: {exc}") from exc

    payload = {
        "config": str(Path(args.config).resolve()),
        "generated_at": time.time(),
        "data": entries,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
    except OSError as exc:  # pragma: no cover
        raise SystemExit(f"Failed to write output JSON: {exc}") from exc


if __name__ == "__main__":  # pragma: no cover
    main()
