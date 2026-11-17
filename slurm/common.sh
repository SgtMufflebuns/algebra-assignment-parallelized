# Shared helpers for the slurm/ scripts. Sourced by them, not run directly.

# Repository root: the folder that contains slurm/.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

die() {
    echo "error: $*" >&2
    exit 1
}

# Make mpirun available. Clusters usually provide MPI as an environment
# module ("module load openmpi"). Set MPI_MODULE if yours has another name.
load_mpi() {
    if command -v mpirun >/dev/null 2>&1; then
        return 0
    fi
    if ! type module >/dev/null 2>&1; then
        for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
            if [ -f "$init" ]; then
                # shellcheck disable=SC1090
                . "$init"
                break
            fi
        done
    fi
    if type module >/dev/null 2>&1; then
        module load "${MPI_MODULE:-openmpi}" >/dev/null 2>&1 || true
    fi
    if ! command -v mpirun >/dev/null 2>&1; then
        echo "error: mpirun not found. Load your MPI module first (e.g. 'module load openmpi')," >&2
        echo "       or set MPI_MODULE to the name shown by 'module avail'." >&2
        return 1
    fi
}
