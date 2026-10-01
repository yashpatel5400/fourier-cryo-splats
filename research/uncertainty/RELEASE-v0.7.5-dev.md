# Candidate-derived uncertainty tests on three Gaussian reconstructions

The 32-page ICML-format development paper now includes scores designed from
all three actually fitted Gaussian maps, without external reference maps for
region or contrast selection. It preserves every outcome: 18,900 scalar
rejection projections and 6,300 cells of actual repeated-group results.

Scale-orthogonal scores improve sensitivity, but 10-25% regional changes remain
difficult and modest viewing uncertainty reduces power. The classical CVaR
comparison reverses the earlier ordering on some candidates. All correct-null
controls remain public, including a largest count of 7/128 with a wide
pointwise interval. These hold calibration fixed; they do not establish
unconditional error over calibration repetitions.

Independent numerical checks reproduce all score directions, counts,
critical values, projections and group vectors, plus selected physical sums
and amplitude polynomials. The archive includes all three original fitted
MRC candidate maps and saved input/source hashes, not just derived arrays.
Three new score-design tests pass. The PDF's new pages are visually checked.

This is known-simulator method development. Viewing/noise calibration,
regional attribution, practical small-error power and novelty remain open.
All three full Claude Fable 5.1 reviews remain rejections. No acceptance
verdict, experimentally calibrated uncertainty or new reconstruction is claimed.

See `research/uncertainty/REPRODUCE-REVISION8.md` for commands and dependencies
on v0.7.2-dev, v0.7.3-dev and v0.7.4-dev. Earlier public assets are unchanged.
