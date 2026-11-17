#!/usr/bin/env bash
# Submit one run to Slurm. Output goes to results/det_n<size>_N<nodes>x<tasks>_<jobid>.out
#
#   slurm/submit.sh --size 16 --nodes 2                       quick check across 2 machines
#   slurm/submit.sh --size 1000 --nodes 4 --time 02:00:00
#   slurm/submit.sh --size 512 --nodes 2 -- --mail-type=END    extra sbatch options after --
set -eo pipefail
# shellcheck source=common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

usage() {
    cat <<'EOF'
Usage: slurm/submit.sh --size N [options] [-- extra sbatch options]

  --size N              matrix size; build it first with slurm/build.sh N   (required)
  --nodes N             machines to use                                      (default 2)
  --tasks-per-node N    MPI ranks on each machine                            (default 1)
  --time HH:MM:SS       time limit; Slurm stops the job when it is reached   (default 01:00:00)
  --partition NAME      Slurm partition, see slurm/check_cluster.sh          (default: cluster default)
EOF
}

size=""
nodes=2
tasks_per_node=1
time_limit="01:00:00"
partition="${PARTITION:-}"
while [ $# -gt 0 ]; do
    case "$1" in
        --size|--nodes|--tasks-per-node|--time|--partition)
            [ $# -ge 2 ] || die "$1 needs a value"
            case "$1" in
                --size) size="$2" ;;
                --nodes) nodes="$2" ;;
                --tasks-per-node) tasks_per_node="$2" ;;
                --time) time_limit="$2" ;;
                --partition) partition="$2" ;;
            esac
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        *) usage >&2; die "unknown option '$1'" ;;
    esac
done

[ -n "$size" ] || { usage >&2; die "--size is required"; }
for value in "$size" "$nodes" "$tasks_per_node"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || die "sizes and counts must be positive integers, got '$value'"
done
command -v sbatch >/dev/null 2>&1 || die "sbatch not found. Run this on the cluster's head node."
[ -x "$REPO_ROOT/build-$size/src/algebra_determinants" ] || die "build-$size is not built. Run: slurm/build.sh $size"

cd "$REPO_ROOT"   # the job runs from the folder it was submitted from
mkdir -p results
name="det_n${size}_N${nodes}x${tasks_per_node}"
args=(--job-name="$name" --nodes="$nodes" --ntasks-per-node="$tasks_per_node"
      --time="$time_limit" --output="results/${name}_%j.out" --parsable)
[ -n "$partition" ] && args+=(--partition="$partition")

job_id="$(sbatch "${args[@]}" "$@" slurm/run_determinants.sbatch "$size")"
job_id="${job_id%%;*}"
echo "Submitted job $job_id: size $size on $nodes node(s) x $tasks_per_node rank(s)"
echo "  output: results/${name}_${job_id}.out"
echo "  status: squeue -j $job_id      cancel: scancel $job_id"
