#!/bin/bash
#SBATCH -p main
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-0%1
#SBATCH -J hedgehog_10k
#SBATCH -o logs/slurm/%x_%A_%a.out
#SBATCH -e logs/slurm/%x_%A_%a.err

set -euo pipefail

REPO_ROOT="/home/nikolenko/work/Projects/hedgehog"
MASTER_INPUT="${MASTER_INPUT:-/mnt/ligandpro/shared_storage/hedgehog_prepared_all_models_v1/all_models.csv}"
MODEL_LIST="${MODEL_LIST:-/mnt/ligandpro/shared_storage/hedgehog_prepared_all_models_v1/model_names.txt}"
RUNS_ROOT="${HEDGEHOG_RUNS_ROOT:-/mnt/ligandpro/shared_storage/hedgehog_runs}"

mkdir -p "${RUNS_ROOT}"

if [[ ! -f "${MASTER_INPUT}" ]]; then
  echo "ERROR: MASTER_INPUT not found: ${MASTER_INPUT}" 1>&2
  exit 2
fi

if [[ ! -f "${MODEL_LIST}" ]]; then
  echo "ERROR: MODEL_LIST not found: ${MODEL_LIST}" 1>&2
  exit 2
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is not set" 1>&2
  exit 2
fi

N_MODELS="$(wc -l < "${MODEL_LIST}" | tr -d ' ')"
if [[ "${N_MODELS}" -eq 0 ]]; then
  echo "ERROR: MODEL_LIST is empty: ${MODEL_LIST}" 1>&2
  exit 2
fi

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= N_MODELS )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (0..$(( ${#INPUTS[@]} - 1 )))" 1>&2
  exit 2
fi

MODEL_NAME="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MODEL_LIST}" | sed 's/\r$//')"
if [[ -z "${MODEL_NAME}" ]]; then
  echo "ERROR: Empty model name for task ${SLURM_ARRAY_TASK_ID}" 1>&2
  exit 2
fi

MODEL_DIR="${RUNS_ROOT}/${MODEL_NAME}"
OUT="${MODEL_DIR}"
MODEL_INPUT="${MODEL_DIR}/model_input.csv"

rm -rf "${MODEL_DIR}"
mkdir -p "${MODEL_DIR}"

python - <<PY
import csv
from pathlib import Path

master = Path("${MASTER_INPUT}")
out = Path("${MODEL_INPUT}")
model = "${MODEL_NAME}"

with master.open("r", encoding="utf-8", newline="") as f_in, out.open("w", encoding="utf-8", newline="") as f_out:
    r = csv.DictReader(f_in)
    if r.fieldnames is None or "smiles" not in r.fieldnames or "model_name" not in r.fieldnames:
        raise SystemExit(f"Bad schema in master input: {master} fields={r.fieldnames}")
    w = csv.writer(f_out)
    w.writerow(["smiles", "model_name"])
    n = 0
    for row in r:
        if (row.get("model_name") or "").strip() != model:
            continue
        smi = (row.get("smiles") or "").strip()
        if not smi:
            continue
        w.writerow([smi, model])
        n += 1

if n == 0:
    raise SystemExit(f"No rows found for model={model} in {master}")
print(f"Wrote {out} rows={n}")
PY

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MOLSCORE_NJOBS="${SLURM_CPUS_PER_TASK:-16}"

cd "${REPO_ROOT}"

echo "MASTER_INPUT=${MASTER_INPUT}"
echo "MODEL_LIST=${MODEL_LIST}"
echo "MODEL_NAME=${MODEL_NAME}"
echo "MODEL_INPUT=${MODEL_INPUT}"
echo "OUT=${OUT}"
echo "CPUS=${SLURM_CPUS_PER_TASK:-}"
echo "HOST=$(hostname)"
echo "START=$(date -Is)"

uv run hedge run \
  --config src/hedgehog/configs/config_slurm_10k.yml \
  --mols "${MODEL_INPUT}" \
  --out "${OUT}"

echo "END=$(date -Is)"
