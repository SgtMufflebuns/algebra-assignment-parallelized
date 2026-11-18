#!/usr/bin/env bash
# Build one program per matrix size. The size is fixed at compile time,
# so each size gets its own folder:
#
#   slurm/build.sh 16 512 1000     ->  build-16/  build-512/  build-1000/
#
# Rebuild after pulling changes or switching branches.
#
# This builds with mpicxx directly instead of CMake, because our cluster's CMake
# (3.22) is older than the 3.24 the project requires. The algebra libraries are
# header-only, so CMake's whole build comes down to one compiler command, and this
# script produces the same program:
#   1. reads cmake/Dependencies.cmake and fetches the same library versions CMake
#      would (including the versions those libraries pin in turn) into build-deps/
#   2. compiles src/main.cpp with mpicxx, using the same flags as CMake's Release build
# If CMakeLists.txt ever gains more source files or settings, or a library with .cpp
# files, update this script to match (or build with CMake 3.24+ instead).
set -eo pipefail
# shellcheck source=common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

[ $# -ge 1 ] || die "usage: slurm/build.sh SIZE [SIZE ...]"
for size in "$@"; do
    [[ "$size" =~ ^[1-9][0-9]*$ ]] || die "size must be a positive integer, got '$size'"
done
load_mpi
command -v mpicxx >/dev/null 2>&1 || die "mpicxx not found. It comes with MPI; check 'module load openmpi'."

deps="$REPO_ROOT/build-deps"    # build-*/ is already in .gitignore
mkdir -p "$deps"
seen=" "
include_flags=()

# Print "name repository commit" for each FetchContent_Declare in a CMake file.
declarations() {
    awk '
        /FetchContent_Declare\(/        { inblock = 1; name = repo = tag = ""; next }
        inblock && /GIT_REPOSITORY/     { repo = $2; next }
        inblock && /GIT_TAG/            { tag = $2; next }
        inblock && /\)/                 { if (name != "") print name, repo, tag; inblock = 0; next }
        inblock && name == "" && NF > 0 { name = $1 }
    ' "$1"
}

# Fetch each declared library at its pinned commit, then its own dependencies.
# Like CMake's FetchContent, the first version declared for a name wins.
fetch() {
    local name repo tag dir url
    while read -r name repo tag; do
        [ -n "$name" ] || continue
        [ "$name" = googletest ] && continue    # only used by the libraries' own unit tests
        case "$seen" in *" $name "*) continue ;; esac
        seen="$seen$name "
        dir="$deps/$name"
        url="${repo/git@github.com:/https://github.com/}"   # public repos, so no GitHub key needed
        if [ "$(git -C "$dir" rev-parse HEAD 2>/dev/null)" != "$tag" ]; then
            echo "== Fetching $name at ${tag:0:7}"
            rm -rf "$dir"
            git clone --quiet "$url" "$dir"
            git -C "$dir" checkout --quiet "$tag"
        fi
        include_flags+=(-I "$dir/include")
        if [ -f "$dir/cmake/Dependencies.cmake" ]; then
            fetch < <(declarations "$dir/cmake/Dependencies.cmake")
        fi
    done
}

fetch < <(declarations "$REPO_ROOT/cmake/Dependencies.cmake")
[ ${#include_flags[@]} -gt 0 ] || die "no libraries found in cmake/Dependencies.cmake"

for size in "$@"; do
    out="$REPO_ROOT/build-$size/src"
    mkdir -p "$out"
    echo "== Compiling size $size"
    mpicxx -O3 -DNDEBUG -std=gnu++20 -DALGEBRA_ASSIGNMENT_SIZE="$size" \
           "${include_flags[@]}" "$REPO_ROOT/src/main.cpp" -o "$out/algebra_determinants"
    echo "== Built $out/algebra_determinants"
done
