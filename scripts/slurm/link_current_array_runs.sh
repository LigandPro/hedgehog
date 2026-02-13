#!/bin/bash
set -euo pipefail

# Create symlinks under RUNS_ROOT/<array_job_id>/ pointing to per-task job-id folders.
#
# Useful for already-submitted arrays where the output folder used SLURM_JOB_ID
# (each task got its own numeric directory), but you want a stable structure:
#   RUNS_ROOT/<array_job_id>/<model_name> -> RUNS_ROOT/<task_job_id>/<taskid>_<model_name>
#
# Example:
#   scripts/slurm/link_current_array_runs.sh 32 /mnt/ligandpro/shared_storage/hedgehog_runs

ARRAY_ID="${1:-}"
RUNS_ROOT="${2:-/mnt/ligandpro/shared_storage/hedgehog_runs}"

if [[ -z "${ARRAY_ID}" ]]; then
  echo "Usage: $0 <array_job_id> [runs_root]" 1>&2
  exit 2
fi

DEST="${RUNS_ROOT}/${ARRAY_ID}"
mkdir -p "${DEST}"

shopt -s nullglob
for job_dir in "${RUNS_ROOT}"/[0-9]*; do
  [[ -d "${job_dir}" ]] || continue
  job_base="$(basename "${job_dir}")"
  [[ "${job_base}" =~ ^[0-9]+$ ]] || continue
  # Skip destination directory itself if it happens to be numeric.
  if [[ "${job_base}" == "${ARRAY_ID}" ]]; then
    continue
  fi

  for run_dir in "${job_dir}"/*_*; do
    [[ -d "${run_dir}" ]] || continue
    run_base="$(basename "${run_dir}")"
    model="${run_base#*_}"
    link="${DEST}/${model}"
    target="../${job_base}/${run_base}"

    if [[ -e "${link}" || -L "${link}" ]]; then
      continue
    fi

    ln -s "${target}" "${link}"
    echo "Linked ${link} -> ${target}"
  done
done
