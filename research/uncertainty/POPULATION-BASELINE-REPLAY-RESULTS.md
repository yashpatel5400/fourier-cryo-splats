# Released population baseline: author-code replay

1 October 2026 UTC. All 23 released two-state conditions complete using the unchanged author functions at commit `b9099e1b7fb3f03207f94d89b153839bfcd6a7c4`. This is a retrospective computational check of published likelihood arrays, not a new cryo-EM reconstruction, CAHRA result or coverage experiment. The original files are unchanged.

The first attempt failed before any data analysis because the isolated environment lacked CVXPY. That log is preserved. Installing the dependency allowed the second attempt to complete in 3.646 seconds. JAX 0.11.2, JAXLIB 0.11.2 and CVXPY 1.9.3 are current dependencies, not the original authors’ environment.

The independent monotone root solver agrees with a separate tanh/log-odds bisection to 1.22e-15. All four primitive checks (interior, both boundaries, flat likelihood) pass. Deconvolution agrees with its analytic constrained solution to 1.55e-15. These checks do not validate the supplied image likelihoods.

Every author run prints an exit message. This is not proof of high-precision optimization: the float32 weights drift slightly off the simplex, by at most 1.82e-05, and the weakest synthetic case differs from the independently optimized population by 0.00785. Normalizing the returned weights changes the comparison only slightly. `independent-check.json` retains normalized gradient gaps and objective differences; the raw summary retains the original diagnostics, including tiny negative objective differences caused by the unnormalized weight sum. No run is silently replaced.

The saved synthetic images have an empirical class-0 fraction of .80074; that is distinct from the generator probability .8. The experimental .8 comparison is a constructed reference from classified subsets, not known biological truth. The two separately named noisy-image likelihood files are byte-identical and therefore do not constitute independent controls.

| Condition | Hard | Soft | Author EM | Independent EM | Deconvolution | Saved author EM |
|---|---:|---:|---:|---:|---:|---:|
| synthetic-0 | 0.800740 | 0.800740 | 0.800728 | 0.800740 | 0.800740 | 0.800740 |
| synthetic-1 | 0.800740 | 0.800740 | 0.800728 | 0.800740 | 0.800740 | 0.800740 |
| synthetic-2 | 0.799950 | 0.799528 | 0.800726 | 0.800742 | 0.800761 | 0.800742 |
| synthetic-3 | 0.773260 | 0.759718 | 0.799422 | 0.799417 | 0.800412 | 0.799417 |
| synthetic-4 | 0.693690 | 0.649571 | 0.797846 | 0.797846 | 0.801938 | 0.797846 |
| synthetic-5 | 0.611520 | 0.559338 | 0.796265 | 0.796268 | 0.797547 | 0.796267 |
| synthetic-6 | 0.562060 | 0.518987 | 0.794482 | 0.794497 | 0.804472 | 0.794494 |
| synthetic-7 | 0.535530 | 0.505497 | 0.791656 | 0.791701 | 0.827604 | 0.791687 |
| synthetic-8 | 0.524890 | 0.501523 | 0.785945 | 0.786111 | 0.934014 | 0.786049 |
| synthetic-9 | 0.524210 | 0.500408 | 0.766605 | 0.774456 | 1.000000 | 0.774290 |
| experimental-default | 0.796218 | 0.795360 | 0.798215 | 0.798228 | 0.796989 | 0.798228 |
| experimental-noisy-images-labelled-ground-truth | 0.752457 | 0.736021 | 0.791800 | 0.791924 | 0.787715 | 0.791797 |
| experimental-noisy-images | 0.752457 | 0.736021 | 0.791800 | 0.791924 | 0.787715 | 0.791797 |
| rotation-10 | 0.633076 | 0.631303 | 0.649705 | 0.649749 | 0.633441 | 0.649709 |
| rotation-2 | 0.787079 | 0.785465 | 0.792636 | 0.792627 | 0.787829 | 0.792621 |
| rotation-4 | 0.751353 | 0.749421 | 0.765196 | 0.765208 | 0.752018 | 0.765203 |
| rotation-6 | 0.712065 | 0.709589 | 0.728998 | 0.729008 | 0.712636 | 0.728992 |
| rotation-8 | 0.668498 | 0.666497 | 0.686470 | 0.686504 | 0.668957 | 0.686473 |
| shift-0 | 0.782938 | 0.781115 | 0.789880 | 0.789894 | 0.783686 | 0.789883 |
| shift-1 | 0.779735 | 0.777467 | 0.787455 | 0.787488 | 0.780477 | 0.787471 |
| shift-2 | 0.775814 | 0.773292 | 0.784807 | 0.784831 | 0.776547 | 0.784805 |
| shift-3 | 0.770265 | 0.768174 | 0.780891 | 0.780911 | 0.770985 | 0.780876 |
| shift-4 | 0.763528 | 0.761045 | 0.775221 | 0.775213 | 0.764232 | 0.775210 |

All numbers are class-0 proportions. `results/uncertainty/development/population-author-replay-v1/comparison.csv` contains full precision; the summary has input hashes, versions, runtimes and all returned weights. Per-condition stdout remains available. The source LICENSE is GPLv3 while its package metadata says MIT; no author implementation was copied into our package. See the [frozen protocol](POPULATION-BASELINE-REPLAY-PROTOCOL.md).

Reproduce by obtaining the pinned author checkout under `background/counting_particles_paper`, creating the isolated environment described in the protocol, and running `scripts/replay_population_baseline.py` with an absent output directory. The script loads only allowed array constructors from the archived pickles. Then run `scripts/check_population_baseline_replay.py`. The whole-stack likelihood/reweighting comparator is now inspected and executed; robustness to unknown poses, noise and templates is still unestablished.
