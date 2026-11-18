#!/usr/bin/env bash
# Turn every results/*.out file into one CSV row:
#
#   slurm/collect.sh > results/summary.csv
#
# Runs that failed or were cut off still get a row; their missing columns are empty.
set -eo pipefail
# shellcheck source=common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

shopt -s nullglob
files=("$REPO_ROOT"/results/*.out)
[ ${#files[@]} -gt 0 ] || die "no results/*.out files yet"

awk '
function reset() {
    job = branch = commit = size = nodes = tpn = ranks = ""
    gen = prod = phase = total = verified = status = ""
    ns = 0
    split("", secs)
}
function emit() {
    print job "," branch "," commit "," size "," nodes "," tpn "," ranks "," \
          gen "," prod "," secs[1] "," secs[2] "," secs[3] "," phase "," \
          total "," verified "," status "," file
}
BEGIN {
    print "job_id,branch,commit,size,nodes,tasks_per_node,ranks,generation_s,product_s," \
          "det_a_s,det_b_s,det_c_s,det_phase_s,total_s,verified,exit_status,file"
    reset()
}
FNR == 1 { if (NR > 1) emit(); reset(); file = FILENAME; sub(/.*\//, "", file) }
$1 == "job_id"         { job = $2 }
$1 == "branch"         { branch = $2 }
$1 == "commit"         { commit = $2 }
$1 == "size"           { size = $2 }
$1 == "nodes"          { nodes = $2 }
$1 == "tasks_per_node" { tpn = $2 }
$1 == "ranks"          { ranks = $2 }
$1 == "exit_status"    { status = $2 }
/^Generation seconds / { gen = $3 }
/^Product seconds /    { prod = $3 }
/^Seconds /            { secs[++ns] = $2 }
/^Determinant seconds:/ { phase = $3 }
/^Total seconds /      { total = $3 }
/^Identity verified /  { verified = $3 }
END { emit() }
' "${files[@]}"
