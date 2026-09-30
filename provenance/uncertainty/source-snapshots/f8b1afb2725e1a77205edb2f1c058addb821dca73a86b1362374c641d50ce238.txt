# Native RELION runtime probe

Declared before execution, 30 September 2026. This is a feasibility measurement,
not an uncertainty comparison or converged ab initio reconstruction. It uses
256 existing EMPIAR-10028 pilot particles, selected without replacement with
seed 670001. No reserved calibration pixels or supplied orientations/shifts
are read into the reconstruction.

Use stock RELION 5.0.1, Bioconda Apple Silicon build `h7c52720_0`, in a separate
environment. Record the package environment, executable hash, version, command,
input hashes, elapsed time and peak child-process memory. Keep all failures.
The probe has five VDAM mini-batches, one C1 class, two CPU threads, initial/final
batch sizes 128/256, and a ten-minute wall-clock ceiling. Initial/final resolution
limits are 60/25 A. The diameter is 80% of the field of view, with the standard
solvent/masking and angular/shift defaults from the initial-model job.

Write a new MRC/optics STAR input, background-center and scale each image using
pixels outside radius 0.43 of the field. Preserve CTF parameters in physical
units. RELION's CTF has the opposite sign to this project's cryoDRGN convention;
therefore multiply raw pixels by minus the recorded data-sign factor. Check the
two analytic CTF conventions across the sampled frequencies, and retain the
small discrepancy from their published wavelength constants. Do not supply
poses, translations, deposited maps, or Gaussian pilot density to RELION.

The official initial-model tutorial and source define the command structure:
https://relion.readthedocs.io/en/release-5.0/SPA_tutorial/InitialModel.html
https://github.com/3dem/relion/blob/5.0.1/src/pipeline_jobs.cpp
https://github.com/3dem/relion/blob/5.0.1/src/ctf.h
https://github.com/3dem/relion/blob/5.0.1/src/ctf.cpp

A successful short run establishes executable compatibility and a local timing
only. It does not establish convergence, calibration, GPU throughput, or a
linear forecast for the full three-stack comparison. A subsequent scientific
baseline requires a separately declared fit and evaluation protocol.
