"""Convert MolScore RAScore pickle model into portable XGBoost JSON model."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="RAScore pickle-to-json converter")
    parser.add_argument("--model-pkl", required=True)
    parser.add_argument("--output-json", required=True)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    model_path = Path(args.model_pkl)
    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(model_path, "rb") as f:
        model = pickle.load(f)

    booster = model.get_booster() if hasattr(model, "get_booster") else model
    booster.save_model(str(output_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
