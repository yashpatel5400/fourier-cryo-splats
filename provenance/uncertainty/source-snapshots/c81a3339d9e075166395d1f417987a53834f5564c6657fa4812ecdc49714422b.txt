# One-movie feasibility pilot for an independent-frame processing graph

1 October 2026 UTC. Declared before accessing this movie's pixel array. This is acquisition/processing development, not an experimental coverage study, and does not replace any frozen cohort or analysis.

Select EMPIAR-10028 `Micrographs/Micrographs_part2/004_movie.mrcs`. The archive's 1,024-byte header was read: 4096 x 4096 x 16, float32, no extended header; total length 1,073,742,848 bytes. Download just this one movie (about 1 GiB), keeping its URL, response headers and SHA256. The existing 472 used source groups do not include `MRC_1901/004_movie.mrcs`. The part2/source-group mapping is supported by the 481/600 group counts, January-19 acquisition date and nearby particle coordinates, but the deposited autopicks are not identical to final recentered coordinates. Retain this identity limitation; do not label this a fully verified independent confirmatory cohort.

Published particle coordinates and CTF logs were inspected only for source identity. Do not use those picks, consensus poses, or full-movie-derived transformations to claim independent inference noise. Any later particle selection/motion/pose pipeline must be specified separately and use only its designated alignment frames. The external reference map also is not experimental ground truth.

Initial analysis is restricted to header and acquisition diagnostics: per-frame means, SDs and quantiles; fixed-grid frame-difference spectra and correlations; unaligned odd/even sums; and a downsampled display. Use every frame and a fixed nonoverlapping patch grid, with no outcome-based patch selection. Difference fields can contain motion/dose signal and detector effects; they are not automatically pure noise or proof of independence. No intervals, pose-radius calibration or density results are obtained in this pilot. Preserve all diagnostic outcomes, including nonstationarity or correlation.

This pilot requires no rented compute or new expenditure. Keep the raw source separate from all processed-particle records and preserve the 116-file frozen calibration lock.
