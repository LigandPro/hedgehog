#!/bin/bash
#SBATCH -p main
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -J hedgehog_all_models
#SBATCH -o logs/slurm/%x_%A.out
#SBATCH -e logs/slurm/%x_%A.err

set -euo pipefail

REPO_ROOT="/home/nikolenko/work/Projects/hedgehog"
MASTER_INPUT="/mnt/ligandpro/shared_storage/hedgehog_prepared_all_models_v1/all_models.csv"
RUNS_ROOT="/mnt/ligandpro/shared_storage/hedgehog_runs"

mkdir -p "${RUNS_ROOT}"
mkdir -p "${REPO_ROOT}/logs/slurm"

OUT="${RUNS_ROOT}/all_models"

rm -rf "${OUT}"
mkdir -p "${OUT}"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MOLSCORE_NJOBS="${SLURM_CPUS_PER_TASK:-32}"
export AIZYNTH_NPROC="${AIZYNTH_NPROC:-8}"
if [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "${SLURM_CPUS_PER_TASK}" -lt "${AIZYNTH_NPROC}" ]; then
  export AIZYNTH_NPROC="${SLURM_CPUS_PER_TASK}"
fi

cd "${REPO_ROOT}"

echo "MASTER_INPUT=${MASTER_INPUT}"
echo "OUT=${OUT}"
echo "CPUS=${SLURM_CPUS_PER_TASK:-}"
echo "AIZYNTH_NPROC=${AIZYNTH_NPROC}"
echo "HOST=$(hostname)"
echo "START=$(date -Is)"

uv run hedge run \
  --config src/hedgehog/configs/config_slurm_10k_32c.yml \
  --mols "${MASTER_INPUT}" \
  --out "${OUT}"

echo "END=$(date -Is)"
