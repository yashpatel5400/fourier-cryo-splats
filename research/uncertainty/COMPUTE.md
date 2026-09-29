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
