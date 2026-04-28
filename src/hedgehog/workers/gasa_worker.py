"""Run GASA scoring in an isolated Python environment."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

_GASA_SCORE_SCRIPT = (
    "import math\n"
    "import pandas as pd\n"
    "import sys\n"
    "from gasa import GASA\n"
    "input_csv, output_csv, smiles_col = sys.argv[1], sys.argv[2], sys.argv[3]\n"
    "df = pd.read_csv(input_csv)\n"
    "if smiles_col not in df.columns:\n"
    "    raise ValueError(f'Missing SMILES column: {smiles_col}')\n"
    "raw_smiles = df[smiles_col].astype(str).tolist()\n"
    "valid_pairs = [(idx, smi) for idx, smi in enumerate(raw_smiles) if smi and smi.lower() not in {'nan', 'none'}]\n"
    "scores = [float('nan')] * len(raw_smiles)\n"
    "if valid_pairs:\n"
    "    valid_smiles = [smi for _, smi in valid_pairs]\n"
    "    sys.argv = [sys.argv[0]]\n"
    "    _, _, hs_probs = GASA(valid_smiles)\n"
    "    for (idx, _), hs_prob in zip(valid_pairs, hs_probs):\n"
    "        try:\n"
    "            scores[idx] = float(hs_prob)\n"
    "        except Exception:\n"
    "            scores[idx] = float('nan')\n"
    "# Upstream returns class-1 (HS=hard-to-synthesize) probability as `neg`.\n"
    "df['gasa_score'] = scores\n"
    "df.to_csv(output_csv, index=False)\n"
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GASA isolated worker launcher")
    parser.add_argument("--worker-python", required=True)
    parser.add_argument("--repo-path", required=True)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--smiles-column", default="smiles")
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=None)
    return parser.parse_args()


def _worker_env(repo_path: str) -> dict[str, str]:
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "").strip()
    if existing:
        env["PYTHONPATH"] = f"{repo_path}{os.pathsep}{existing}"
    else:
        env["PYTHONPATH"] = repo_path
    return env


def main() -> int:
    args = _parse_args()
    del args.batch_size, args.num_workers

    repo_path = Path(args.repo_path).expanduser()
    if not repo_path.exists():
        print(f"GASA repo path does not exist: {repo_path}")
        return 1

    command = [
        args.worker_python,
        "-c",
        _GASA_SCORE_SCRIPT,
        args.input_csv,
        args.output_csv,
        args.smiles_column,
    ]

    try:
        subprocess.run(
            command,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
            cwd=str(repo_path),
            env=_worker_env(str(repo_path)),
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
