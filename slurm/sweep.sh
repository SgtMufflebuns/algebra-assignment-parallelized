#!/usr/bin/env bash
# Submit every combination of sizes and node counts, e.g. for the break-even study:
#
#   slurm/sweep.sh --sizes "64 128 256 512 1000" --nodes "1 2 4" --repeat 3
#
# 1 node x 1 task is the single-process baseline, using the same program.
# Options other than --sizes, --nodes and --repeat are passed on to slurm/submit.sh.
set -eo pipefail
# shellcheck source=common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

sizes=""
node_counts=""
repeat=1
pass=()
while [ $# -gt 0 ]; do
    case "$1" in
        --sizes)  [ $# -ge 2 ] || die "--sizes needs a value";  sizes="$2";       shift 2 ;;
        --nodes)  [ $# -ge 2 ] || die "--nodes needs a value";  node_counts="$2"; shift 2 ;;
        --repeat) [ $# -ge 2 ] || die "--repeat needs a value"; repeat="$2";      shift 2 ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        *) pass+=("$1"); shift ;;
    esac
done
[ -n "$sizes" ] && [ -n "$node_counts" ] || die "usage: slurm/sweep.sh --sizes \"64 128\" --nodes \"1 2 4\" [--repeat N] [submit.sh options]"
[[ "$repeat" =~ ^[1-9][0-9]*$ ]] || die "--repeat must be a positive integer"

# Check every build exists before submitting anything.
missing=""
for size in $sizes; do
    [ -x "$REPO_ROOT/build-$size/src/algebra_determinants" ] || missing="$missing $size"
done
[ -z "$missing" ] || die "not built:$missing. Run: slurm/build.sh$missing"

count=0
for size in $sizes; do
    for nodes in $node_counts; do
        for _ in $(seq "$repeat"); do
            "$REPO_ROOT/slurm/submit.sh" --size "$size" --nodes "$nodes" "${pass[@]}"
            count=$((count + 1))
        done
    done
done
echo "Submitted $count job(s). Watch them with: squeue -u $USER"
