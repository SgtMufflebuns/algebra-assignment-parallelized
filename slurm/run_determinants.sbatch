#!/usr/bin/env bash
#SBATCH --job-name=determinants
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=results/%x_%j.out
#
# One run of build-<SIZE>/src/algebra_determinants.
#
# On the cluster, submit it with slurm/submit.sh (which names the job and the
# output file). On a laptop without Slurm, run it directly:
#
#   NP=3 bash slurm/run_determinants.sbatch 16
#
# Optional environment variables (they carry into the job when you submit):
#   LAUNCHER     mpirun (default) or srun
#   MPIRUN_ARGS  extra mpirun options, e.g. "--oversubscribe"
#   SRUN_MPI     MPI type for srun (default pmix; see 'srun --mpi=list')
#   MPI_MODULE   module that provides MPI (default openmpi)
#   NP           number of ranks when running without Slurm (default 1)
set -eo pipefail

size="${1:?usage: run_determinants.sbatch SIZE}"

# Under Slurm this file runs from a spool copy, so locate the repo from where
# the job was submitted. Without Slurm, use this file's own location.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    repo="$SLURM_SUBMIT_DIR"
else
    repo="$(cd "$(dirname "$0")/.." && pwd)"
fi
# shellcheck source=common.sh
source "$repo/slurm/common.sh"

exe="$REPO_ROOT/build-$size/src/algebra_determinants"
[ -x "$exe" ] || die "$exe not found. Build it first: slurm/build.sh $size"
if [ "$REPO_ROOT/src/main.cpp" -nt "$exe" ]; then
    echo "warning: src/main.cpp changed after build-$size was built. Rebuild with: slurm/build.sh $size" >&2
fi
load_mpi

np="${SLURM_NTASKS:-${NP:-1}}"
launcher="${LAUNCHER:-mpirun}"
read -r -a mpirun_args <<< "${MPIRUN_ARGS:-}"
case "$launcher" in
    mpirun) launch=(mpirun -np "$np" "${mpirun_args[@]}") ;;
    srun)   launch=(srun --mpi="${SRUN_MPI:-pmix}") ;;
    *)      die "LAUNCHER must be 'mpirun' or 'srun', got '$launcher'" ;;
esac

# Header lines are "key value" so slurm/collect.sh can turn them into CSV columns.
echo "== Job"
echo "job_id ${SLURM_JOB_ID:-local}"
echo "branch $(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
echo "commit $(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "size $size"
echo "nodes ${SLURM_JOB_NUM_NODES:-1}"
echo "tasks_per_node ${SLURM_NTASKS_PER_NODE:-$np}"
echo "ranks $np"
echo "node_list ${SLURM_JOB_NODELIST:-$(hostname)}"
echo "launcher ${launch[*]}"
echo "mpi $(mpirun --version 2>&1 | head -1)"
echo "started $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Record which host each rank ran on, using the same launch settings as the real run.
echo "== Placement (rank -> host)"
if [ "$launcher" = mpirun ]; then
    "${launch[@]}" --tag-output hostname | sort -t, -k2,2n || echo "(placement check failed)"
else
    "${launch[@]}" -l hostname | sort -n || echo "(placement check failed)"
fi

echo "== Program output"
set +e
"${launch[@]}" "$exe"
status=$?
set -e

echo "== Done"
echo "finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "exit_status $status"
exit "$status"
