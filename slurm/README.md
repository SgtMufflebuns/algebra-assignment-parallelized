# Running on the Slurm cluster

These scripts build the determinant program and run it across several machines
(nodes) on our AWS Slurm cluster. They work with any branch of this repo, because
every branch builds the same `build-<size>/src/algebra_determinants` program.

| Script | What it does |
|---|---|
| `slurm/check_cluster.sh` | Shows partitions, node sizes, MPI version and build tools. Changes nothing. |
| `slurm/build.sh SIZE...` | Builds one program per matrix size, into `build-<size>/` (no CMake needed) |
| `slurm/submit.sh` | Submits one run. Output goes to `results/` |
| `slurm/sweep.sh` | Submits every combination of sizes and node counts |
| `slurm/collect.sh` | Turns all `results/*.out` files into one CSV |
| `slurm/run_determinants.sbatch` | The job itself. `submit.sh` uses it; you rarely call it directly |

## 1. One-time setup

On the cluster's head node (where you log in):

```sh
git clone https://github.com/SgtMufflebuns/algebra-assignment-parallelized.git
cd algebra-assignment-parallelized
git checkout 1D-Col-Blocked          # or whichever branch you want to run
slurm/check_cluster.sh
```

Cloning with `https://` means you don't need a GitHub key on the cluster.

In the `check_cluster.sh` output, look at:
- **Slurm:** the partition names and how many nodes each has.
- **MPI:** if it says `mpirun not found`, run `module avail`, find the MPI module,
  and add `export MPI_MODULE=<its name>` to your `~/.bashrc`.
- **Home folder:** it should be a shared filesystem (`nfs`, `lustre`, ...).
  Compute nodes run the program from your home folder.

## 2. Build

```sh
slurm/build.sh 16 512 1000
```

This makes `build-16/`, `build-512/` and `build-1000/`. The matrix size is fixed
when the program is compiled, so each size needs its own build. **Rebuild after
pulling changes or switching branches**; the job prints a warning if the source
is newer than the build.

To compare two branches, clone the repo twice into two folders, so their builds
and results stay separate.

**Why this doesn't use CMake:** the project needs CMake 3.24 or newer, and our
cluster has 3.22. The algebra libraries are header-only, so `build.sh` builds
the same program directly:

1. it reads `cmake/Dependencies.cmake` and downloads the same library versions
   CMake would use into `build-deps/` (once; later builds reuse them)
2. it compiles `src/main.cpp` with `mpicxx`, using CMake's Release settings
   (`-O3 -DNDEBUG -std=gnu++20`)

The result is the same program, with the same library versions and settings, as
the README's CMake build. Those CMake commands still work anywhere CMake 3.24+
is installed. In the report, say that cluster builds were compiled directly with
`mpicxx` using the same settings and library versions as the CMake build.

**Keep it in sync:** `build.sh` mirrors today's CMake build. If `CMakeLists.txt`
gains more source files or settings, or a library starts shipping `.cpp` files
(for example, a cached field that isn't header-only), update `build.sh` to match,
or build those runs with CMake 3.24+.

## 3. Quick test across two machines

```sh
slurm/submit.sh --size 16 --nodes 2
```

When it finishes, open the output file it names under `results/` and check:
- `Identity verified true` near the end
- the `Placement` section lists **two different hostnames**, which proves the run
  used two machines

## 4. Real runs

```sh
slurm/submit.sh --size 1000 --nodes 4 --time 02:00:00
```

| Option | Meaning | Default |
|---|---|---|
| `--size N` | matrix size (build it first) | required |
| `--nodes N` | number of machines | 2 |
| `--tasks-per-node N` | MPI ranks on each machine | 1 |
| `--time HH:MM:SS` | time limit; **Slurm kills the job when it runs out** | 01:00:00 |
| `--partition NAME` | partition from `check_cluster.sh` | cluster default |
| `-- <sbatch options>` | anything after `--` goes straight to `sbatch` | |

Total ranks = nodes × tasks-per-node. One rank per node is the cleanest setup.
It spreads the work across machines, and each machine only holds one copy of the
matrices.

**Setting `--time`:** run time grows with the cube of the size, so doubling the
size takes about 8× longer. A single-process run at size 2000 took about 65
minutes on one of our machines, and more than half of that was `C = A * B`, which
every rank still computes. Give large sizes plenty of time.

**Memory:** every rank currently holds A, B, C and a working copy, roughly
180 bytes × size² per rank: about 0.7 GB at size 2000 and 3 GB at size 4000.
More ranks per node means more memory per node.

## 5. Break-even study (many runs at once)

```sh
slurm/sweep.sh --sizes "64 128 256 512 1000" --nodes "1 2 4" --repeat 3 --time 02:00:00
```

This submits every size × node-count combination, 3 times each. `--nodes 1` with
the default 1 rank per node is the **single-process baseline**, using the same
program. Report the median of the repeats, because timings vary between runs.

## 6. Watching jobs

```sh
squeue -u $USER                  # your jobs
scancel <job id>                 # cancel one
tail -f results/<file>.out       # follow a running job's output
```

The `ST` column in `squeue`:
- `PD`: waiting for nodes. `squeue -u $USER -o "%i %T %r"` shows why.
- `CF`: nodes are starting up. On AWS this can take several minutes the first
  time, while the machines boot.
- `R`: running.

## 7. Collecting results

```sh
slurm/collect.sh > results/summary.csv
```

This gives one row per output file: branch, commit, size, nodes, ranks, the
generation, product and determinant times, total time, whether the identity
verified, and the exit status. Failed or timed-out runs still get a row, with
blanks where values are missing.

How to read the determinant columns, which differ by branch:
- `1D-Col-Blocked`: `det_a_s`, `det_b_s` and `det_c_s` are each the time for **all
  ranks together** to compute that one determinant.
- `feature/mpi-orchestration`: each `det_*_s` is the time on the single rank that
  computed it. `det_phase_s` is the time for all three together.

## 8. Running on your own computer (no Slurm)

The same job script runs locally with plain `mpirun`:

```sh
slurm/build.sh 16
mkdir -p results
NP=3 bash slurm/run_determinants.sbatch 16 > results/local_n16_np3.out
```

On a laptop, asking for more ranks than you have cores needs
`MPIRUN_ARGS=--oversubscribe`.

## Fair timing checklist

- Always submit with `submit.sh` or `sweep.sh`, never run on the head node. Jobs
  use `--exclusive`, so no one else's job shares your nodes during a run.
- Only compare runs from the same partition (same machine type).
- Record the branch and commit. Every output file starts with them.
- Repeat each setup and use the median.

## Troubleshooting

| Problem | Fix |
|---|---|
| `build-N is not built` | `slurm/build.sh N` |
| `CMake 3.24 or higher is required` | You ran the README's CMake commands. On the cluster, use `slurm/build.sh N` instead (see section 2) |
| `mpicxx not found` | Same fix as `mpirun not found` below; both come with MPI |
| `mpirun not found` | `export MPI_MODULE=<name>` (see `module avail`), then submit again |
| Job stays `PD` | Ask for fewer nodes than the partition has (`check_cluster.sh` shows the count), or wait |
| Job ends with `DUE TO TIME LIMIT` | Submit again with a larger `--time` |
| Placement shows only one hostname | The cluster's MPI isn't using the Slurm allocation. Try `LAUNCHER=srun slurm/submit.sh ...` |
| Run hangs at the start across nodes | Try `MPIRUN_ARGS="--mca btl_tcp_if_include ens5" slurm/submit.sh ...` (find the interface name with `ip addr` on a compute node) |

Environment variables such as `LAUNCHER`, `MPIRUN_ARGS` and `MPI_MODULE` carry
into the job when you submit.
