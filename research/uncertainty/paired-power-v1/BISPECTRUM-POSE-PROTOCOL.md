# Continuous-orientation falsification of moment contrasts

Post-outcome development, 1 October 2026 UTC. The completed 24-case moment
hull screen has twelve iteration-limited removal fits and twelve exact controls.
10028/10049 contrasts are nonzero; all four 10076 fits return no positive
separator. No convergence, global null maximum or density coverage is claimed.

For every nonzero saved direction on all three stacks, search continuous
SO(3) using twelve largest-score orientations from the separate 10,000-view
catalog plus eight Haar draws with seed 261011+dataset. Profile amplitude
analytically over the saved interval at every evaluation, including its
stationary point. Translation is absent because each feature is invariant.
Optimize a rotation-vector chart around each start with L-BFGS-B, component
bounds [-pi,pi], at most 150 iterations / 2,000 evaluations. Preserve all
outcomes and choose the larger of the initial and final guide scores.

The 256-padded Fourier interpolator is only a search guide. Recompute each
selected view by the physical cell Fourier operator, then reprofile amplitude.
Store rotations, scores, amplitude choices and exact complex means. The
largest verified score is a lower value for the continuous null maximum,
never an upper bound. If it exceeds the true-map alternative expectation,
the saved contrast has no positive uniform mean separation on that enlarged
set. If it stays below, separation remains unresolved without a continuous
upper certificate. A finite number of starts cannot prove absence of a
counterexample. The all-zero directions are explicitly recorded as skipped.

For each tested contrast independently verify its worst selected view using
direct sums over all 64^3 physical density cells, rather than the NUFFT or
interpolator. This verification is declared before the pose-search outcomes;
ordinary floating-point error remains distinct from validated arithmetic.
Do not run confidence/power trials or launch rentals from a sampled positive
margin alone. A more precise optimizer can be considered only after preserving
the current directions and failures; this protocol does not authorize silently
changing the previous source or summary.
