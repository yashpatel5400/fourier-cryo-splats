# Compute plan (development estimate)

The user authorized using the current Mac's compute resources. It is an M4 Pro
with 12 logical CPUs and 24 GiB unified memory. The initial numerical, synthetic,
and reduced-resolution experiments run locally. Caffeinate prevents idle sleep.

No cloud instance has been requested or purchased. On 2026-09-29 the user asked
for a bounded estimate for possible Vast.ai use. The first proposed cloud batch
has a **$100 total ceiling**, conditional on selecting an offer whose complete
charges fit that ceiling: one ordinary RTX 4090 (24 GB), up to 100 GPU-hours at
no more than $0.65/hour, and at most $35 combined storage/transfer. This is a
planning cap, not a claim that the full research program needs only 100 hours.
Local profiling and exact baseline configurations must precede a rental request.
No additional batch is implied by this estimate.

Official pricing checked: https://vast.ai/pricing and the RTX 4090 marketplace
https://cloud.vast.ai/?gpu_name=RTX+4090. Listings are variable. Exact compute,
storage, transfer, shutdown/deletion charges must be checked on the selected
offer before provisioning. Preserve checkpoints locally before releasing storage.

Preferred initial hardware: one 24 GB CUDA GPU, at least 64 GB host RAM and
150–250 GB scratch if the final data subset fits. Consider a 48 GB A6000 only if
measured memory use requires it and its full cost fits the same approved batch.

## Local runtime notes

The small physical certificate pilot runs in seconds. A first independent-map
dictionary audit used 256 images per split, 56 Fourier samples per image and
up to 2,110 real coefficients. Fitting/rendering a noiseless generator and real
particles took approximately 24–25 seconds per dataset for the finest dictionary.
These are low-frequency development timings, not whole-dataset uncertainty or
neural-baseline timings, and cannot yet justify a precise GPU request.

The stock macOS FINUFFT 2.5.1 wheel bundles a second OpenMP runtime that conflicts
with PyTorch in the same process. The full test suite exposed the native abort.
Rebuild FINUFFT without OpenMP for this Mac, rather than enabling duplicate
OpenMP runtimes:

```sh
CMAKE_BUILD_PARALLEL_LEVEL=8 .venv/bin/pip install --force-reinstall --no-deps \
  --no-binary finufft -Ccmake.define.FINUFFT_USE_OPENMP=OFF finufft==2.5.1
```

The resulting wheel SHA256 was
`393556bf902dcfc830540fcd3c6aeba1725441d127d056120a3d004510079ed5`.
The independent-map generator requests one FINUFFT thread on macOS; BLAS and
other numerical kernels remain parallel. Linux CUDA environments need their
own dependency check. The build option is documented in FINUFFT's CMake source:
https://github.com/flatironinstitute/finufft/blob/master/CMakeLists.txt.

A stock cryoDRGN 4.3.1 `train_nn` CPU timing probe used 512 experimental 64x64
particles, three hidden layers of width 128 (74,498 parameters), batch size 8,
and one epoch. Its own log reports 1.32 seconds of training and 2.42 seconds
including reconstruction output; process startup/preparation brought the wrapper
measurement to 8.25 seconds. This is a deliberately small network and **not a
converged reconstruction or an uncertainty baseline**. The installed stock code
chooses CUDA or CPU and does not select MPS. These timings support continuing
small neural development runs locally; they do not predict high-resolution or
ab initio GPU time.

A second bounded probe executed the stock cryoDRGN 4.3.1 `abinit` command with
no supplied poses, homogeneous latent dimension zero, the default three-layer
width-256 Hartley network, and 32 old development particles. On this Mac it
completed pretraining, one hierarchical pose-search epoch, and one pose-SGD
epoch in 49.40 seconds including process startup. The training log reports
39.03 seconds for pose search and 0.40 seconds for pose SGD. Batch sizes were
two and eight, respectively, with four CPU threads. Exact arguments and source
snapshots are in `results/uncertainty/development/abinit-runtime-probe/profile.json`.
This establishes that the stock unknown-pose path executes locally. It neither
demonstrates convergence nor measures GPU throughput. Linear particle-count
extrapolation from this tiny run would ignore changing search difficulty,
batching, memory, and convergence, so it is not a cloud cost quotation.

## Larger continuous audit and transform reuse

A physical 10 Å central-target fit on 1,024 particles from 10049, Fourier
radius 12 and order-80 quadrature completed its solve/reference checks in
3,583 seconds under concurrent Mac load. Its fixed-pose sum-objective gap is
0.000992; this is not a pose-robust fit. The separate pose audit completed in
1,491 seconds at 1.08 GB peak memory, but returns the no-data interval at one
degree and 0.5 Å because its cubic remainder dominates.

Reusable CPU FINUFFT plans were checked against fresh plans on the same real
acquisition coordinates and seeded coefficients. Outputs agreed to working
precision. With 1,024 particles, radius 12 and order 80, mean warm matvec time
improved by factors 1.111 for the observation Gram and 1.055 for the pose
operator (three calls per variant, the first treated as cold). At 128 particles
the factors were 1.317 and 1.009. These concurrent-load probes are not isolated
performance benchmarks and do not justify a large speedup claim. Optional
subclasses preserve the frozen implementations; they are used in the new
two-pose development test, while existing running jobs retain their original
source snapshots. The API is documented at
https://finufft.readthedocs.io/en/latest/python.html.

Current GPU documentation includes type-3 transforms:
https://finufft.readthedocs.io/en/latest/python_gpu.html. Any GPU port still
requires double-precision conformance and measured-memory tests on actual
hardware. No CUDA speedup or GPU-hour requirement has been measured here,
and no GPU rental has been made.
