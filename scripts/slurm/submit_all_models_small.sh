#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MASTER_INPUT="${MASTER_INPUT:-/mnt/ligandpro/shared_storage/hedgehog_prepared_all_models_v1/all_models.csv}"
MODEL_LIST="${MODEL_LIST:-/mnt/ligandpro/shared_storage/hedgehog_prepared_all_models_v1/model_names.txt}"
MAX_CONCURRENT="${MAX_CONCURRENT:-6}"

cd "${REPO_ROOT}"

mkdir -p logs/slurm

if [[ ! -f "${MASTER_INPUT}" ]]; then
  echo "ERROR: MASTER_INPUT not found: ${MASTER_INPUT}" 1>&2
  exit 2
fi

if [[ ! -f "${MODEL_LIST}" ]]; then
  echo "MODEL_LIST not found, generating: ${MODEL_LIST}"
  python - <<PY
import csv
from pathlib import Path

master = Path("${MASTER_INPUT}")
out = Path("${MODEL_LIST}")

names = set()
with master.open("r", encoding="utf-8", newline="") as f:
    r = csv.DictReader(f)
    if r.fieldnames is None or "model_name" not in r.fieldnames:
        raise SystemExit(f"Bad schema: {master} fields={r.fieldnames}")
    for row in r:
        m = (row.get("model_name") or "").strip()
        if m:
            names.add(m)

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", encoding="utf-8") as f:
    for n in sorted(names):
        f.write(n + "\\n")
print(f"Wrote {out} models={len(names)}")
PY
fi

N="$(wc -l < "${MODEL_LIST}" | tr -d ' ')"
if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: MODEL_LIST is empty: ${MODEL_LIST}" 1>&2
  exit 2
fi

MAX_INDEX="$(( N - 1 ))"

echo "Submitting job array for ${N} inputs (0..${MAX_INDEX})"
echo "MAX_CONCURRENT=${MAX_CONCURRENT}"
sbatch \
  --export=ALL,MASTER_INPUT="${MASTER_INPUT}",MODEL_LIST="${MODEL_LIST}" \
  --array="0-${MAX_INDEX}%${MAX_CONCURRENT}" \
  scripts/slurm/hedgehog_array_small.sbatch
