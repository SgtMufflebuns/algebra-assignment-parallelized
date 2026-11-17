#!/usr/bin/env bash
# Build one executable per matrix size. The size is fixed at compile time,
# so each size gets its own folder:
#
#   slurm/build.sh 16 512 1000     ->  build-16/  build-512/  build-1000/
#
# Rebuild after pulling changes or switching branches.
set -eo pipefail
# shellcheck source=common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

[ $# -ge 1 ] || die "usage: slurm/build.sh SIZE [SIZE ...]"
for size in "$@"; do
    [[ "$size" =~ ^[1-9][0-9]*$ ]] || die "size must be a positive integer, got '$size'"
done

load_mpi

# The algebra libraries are listed with git@github.com: (SSH) addresses, which
# fail on machines without a GitHub SSH key. They are public, so fetch them over
# HTTPS instead. This applies to this script only; your git config is untouched.
git_version="$(git --version | awk '{print $3}')"
if [ "$(printf '%s\n' 2.31 "$git_version" | sort -V | head -1)" != 2.31 ]; then
    die "git $git_version is too old for this script (need 2.31+). Instead run once:
       git config --global url.\"https://github.com/\".insteadOf \"git@github.com:\""
fi
export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0="url.https://github.com/.insteadOf"
export GIT_CONFIG_VALUE_0="git@github.com:"

for size in "$@"; do
    dir="$REPO_ROOT/build-$size"
    echo "== Building size $size in $dir"
    # shellcheck disable=SC2086  # EXTRA_CMAKE_ARGS is intentionally split into words
    cmake -S "$REPO_ROOT" -B "$dir" -DCMAKE_BUILD_TYPE=Release \
          -DALGEBRA_ASSIGNMENT_SIZE="$size" ${EXTRA_CMAKE_ARGS:-}
    cmake --build "$dir" --parallel
    echo "== Built $dir/src/algebra_determinants"
done
