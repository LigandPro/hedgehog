#!/bin/bash
#SBATCH -p main
#SBATCH -c 32
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -J hedgehog_resume_from_synthesis
#SBATCH -o logs/slurm/%x_%A.out
#SBATCH -e logs/slurm/%x_%A.err

set -euo pipefail

REPO_ROOT="/home/nikolenko/work/Projects/hedgehog"
MASTER_INPUT="/mnt/ligandpro/shared_storage/hedgehog_prepared_all_models_v1/all_models.csv"
OUT="/mnt/ligandpro/shared_storage/hedgehog_runs/all_models"

mkdir -p "${REPO_ROOT}/logs/slurm"
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
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
echo "AIZYNTH_NPROC=${AIZYNTH_NPROC}"
echo "HOST=$(hostname)"
echo "START=$(date -Is)"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi -L || true
fi

# Resume from existing synthesis results:
# - Re-run docking without a hard timeout
# - Re-run docking filters
# - Re-run final descriptors
# - Regenerate the full HTML report
rm -rf "${OUT}/stages/05_docking" "${OUT}/stages/06_docking_filters" "${OUT}/stages/07_descriptors_final" "${OUT}/output"

uv run hedge run \
  --config src/hedgehog/configs/config_slurm_10k_32c.yml \
  --mols "${MASTER_INPUT}" \
  --out "${OUT}" \
  --stage docking

uv run hedge run \
  --config src/hedgehog/configs/config_slurm_10k_32c.yml \
  --mols "${MASTER_INPUT}" \
  --out "${OUT}" \
  --stage docking_filters

uv run hedge run \
  --config src/hedgehog/configs/config_slurm_10k_32c.yml \
  --mols "${MASTER_INPUT}" \
  --out "${OUT}" \
  --stage final_descriptors

uv run hedgehog report "${OUT}"

echo "END=$(date -Is)"
