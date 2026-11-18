#!/usr/bin/env bash
# Print what the other slurm/ scripts need to know about this machine and cluster.
# Read-only: changes nothing and submits no jobs.
#
#   slurm/check_cluster.sh

# shellcheck source=common.sh
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

section() { printf '\n== %s ==\n' "$1"; }

section "This machine"
echo "host $(hostname)   arch $(uname -m)"
if [ -r /etc/os-release ]; then
    grep '^PRETTY_NAME=' /etc/os-release | cut -d= -f2 | tr -d '"'
else
    uname -sr
fi

section "Slurm"
if command -v sinfo >/dev/null 2>&1; then
    echo "partition | nodes | CPUs per node | memory per node (MB) | state | node names"
    sinfo -o "%P | %D | %c | %m | %t | %N"
    echo
    echo "MPI types srun supports (only matters if you use LAUNCHER=srun):"
    srun --mpi=list 2>&1 | head -10 | sed 's/^/  /'
else
    echo "sinfo not found: Slurm is not available on this machine."
fi

section "MPI"
if load_mpi 2>/dev/null; then
    echo "mpirun: $(command -v mpirun)"
    mpirun --version 2>&1 | head -1
else
    echo "mpirun not found, and no '${MPI_MODULE:-openmpi}' module could be loaded."
fi
if type module >/dev/null 2>&1; then
    echo "MPI modules available:"
    module avail 2>&1 | grep -io '[^ ]*mpi[^ ]*' | sort -u | sed 's/^/  /'
fi

section "Build tools (slurm/build.sh needs mpicxx, a C++20 compiler such as GCC 10+, and git; CMake is not used)"
for tool in mpicxx c++ git cmake; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "$tool: $("$tool" --version 2>&1 | head -1)"
    else
        echo "$tool: NOT FOUND"
    fi
done

section "Home folder (compute nodes only see your builds if this is shared, e.g. nfs or lustre)"
if [ "$(uname -s)" = Linux ]; then
    df -hT "$HOME" | tail -1
else
    df -h "$HOME" | tail -1
fi
