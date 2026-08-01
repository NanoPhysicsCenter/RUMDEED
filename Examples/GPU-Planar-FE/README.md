# Planar field emission on a GPU (OpenACC)

The same simulation as `Planar-FE/2.0eV`, set up to run on an NVIDIA GPU
through the OpenACC offload path. The input file is a normal mode-10 deck
with one addition: `MH_BATCH = .True.` (see below).

## Building

The GPU path needs the NVIDIA HPC SDK compiler (`nvfortran`):

```sh
cd build
make FC=nvfortran ACC=gpu
```

On a machine without a GPU visible at compile time (e.g. a cluster login
node) `nvfortran` cannot auto-detect the target architecture, so give the
compute capability explicitly, for example for an A100:

```sh
make FC=nvfortran ACC=gpu GPU_CC=cc80
```

Other useful variants:

- `make FC=nvfortran ACC=gpu OPENMP=yes` — OpenACC kernels on the GPU,
  remaining loops threaded with OpenMP on the CPU (recommended).
- `make FC=nvfortran ACC=multicore OPENMP=no` — runs the same OpenACC
  kernels on the CPU cores; useful to test the OpenACC code path on a
  machine without a GPU.

## Running

Run exactly like the CPU examples:

```sh
cd Examples/GPU-Planar-FE
../../build/RUMDEED.out
```

Useful environment variables from the NVIDIA OpenACC runtime:

- `NV_ACC_TIME=1` — per-kernel timing report on stderr at exit.
- `NV_ACC_NOTIFY=1` — print every kernel launch (debugging).
- `CUDA_VISIBLE_DEVICES=0` — select a GPU on multi-GPU nodes.

### Example Slurm job script (A100 cluster)

```sh
#!/bin/bash -l
#SBATCH --job-name=rumdeed-gpu
#SBATCH --partition=gpu-2xA100
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

# Record the GPU state; degraded GPUs run silently ~3x slower.
nvidia-smi --query-gpu=name,clocks.sm,clocks.max.sm --format=csv

./RUMDEED.out
```

Build on the login node (with `GPU_CC=...`), not on the compute nodes.

## What runs on the GPU

The O(N^2) particle-particle acceleration kernel and the batched emission
surface-field evaluations are offloaded (`mod_verlet.F90`). This is
implemented for the **planar** geometry (emission modes 1, 9, 10) and the
**hyperboloid tip** (mode 3) — so the `Tip-FE` and `Tip-GTF` examples also
run on the GPU with this build, unchanged. The cylindrical-tip (11) and
torus (12) modes interpolate their fields from a mesh with kd-tree
searches, which cannot run inside a GPU kernel; those modes automatically
fall back to the OpenMP CPU path at runtime, they do not fail.

## `MH_BATCH`

Field emission (mode 10) samples emission positions with Metropolis-
Hastings chains. With `MH_BATCH = .True.` all chains advance in lockstep
and the surface field for every chain is evaluated in a single batched GPU
kernel per jump instead of one small call per chain — this is where most
of the GPU speed-up of an emission-heavy run comes from (about an order of
magnitude on an A100 for this deck). The batched chains all see the
particle configuration from the start of the time step, so results are
statistically identical to the serial default but not bit-identical run
for run; the default is `.False.` so that existing CPU decks keep their
exact behaviour.

`MH_BATCH` only affects mode 10 (`mod_field_emission_v2`); the samplers of
the other emission modes are unaffected by it.
