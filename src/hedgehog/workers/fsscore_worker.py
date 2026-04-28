"""Run FSScore scoring in an isolated Python environment."""

from __future__ import annotations

import argparse
import os
import subprocess


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="FSScore isolated worker launcher")
    parser.add_argument("--worker-python", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--smiles-column", default="smiles")
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=None)
    parser.add_argument("--graph-datapath", default=None)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    command = [
        args.worker_python,
        "-m",
        "fsscore.score",
        "--model_path",
        args.model_path,
        "--data_path",
        args.input_csv,
        "--compound_cols",
        args.smiles_column,
        "--save_filepath",
        args.output_csv,
        "--batch_size",
        str(args.batch_size),
    ]
    if args.num_workers is not None:
        command.extend(["--num_workers", str(args.num_workers)])
    if args.graph_datapath:
        command.extend(["--graph_datapath", args.graph_datapath])

    try:
        env = os.environ.copy()
        env.setdefault("CUDA_VISIBLE_DEVICES", "")
        subprocess.run(
            command,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
            env=env,
        )
        return 0
    except subprocess.CalledProcessError as exc:
        output = (exc.stderr or exc.stdout or "").strip()
        if output:
            print(output)
        return exc.returncode if exc.returncode != 0 else 1
    except Exception as exc:  # noqa: BLE001
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
