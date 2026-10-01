# Nested maxima and tail risk: targeted primary reading

1 October 2026 UTC. These are methodological prior-art notes, not additional
fully surveyed cryo-EM papers. Source versions and downloaded-file hashes are
in `provenance/uncertainty/nested-risk-reading-downloads.json`. PDFs remain local.

**Giles and Goda, arXiv:1708.05531v4 (2018).**
https://arxiv.org/abs/1708.05531v4
Read the nested formulation, section 2.2 antithetic correction, and assumptions
1-3/Theorem 3; not the full proof or application code. Their correction couples
two half-sample maxima with the full-sample maximum and cancels when the same
decision wins. Fast variance decay needs assumptions controlling probability
near decision ties and how competing conditional means separate. In our use,
amplitude cells play the role of decisions and neighboring cells can be nearly
tied. Those assumptions are not established by bounded event values alone.
We cannot import the complexity rate without additional analysis.

**Giles and Haji-Ali, arXiv:1802.05016v2 (2019).**
https://arxiv.org/abs/1802.05016v2
Read Assumption 2.1 and section 4's CVaR positive-part construction/Theorem 4.1,
not all adaptive-sampling proofs. The antithetic cancellation also applies to
positive parts of inner conditional expectations. Their rate theorem includes
regularity near the standardized decision boundary and conditional-moment/
variance conditions. It supplies relevant prior art for efficient nested risk
estimation, not a ready-made finite-sample high-confidence bound for our
viewing-cap and amplitude-profiled event probability.

**Project consequence (our inference).** The expectation of a maximum of noisy
amplitude-cell averages can exceed the maximum of their conditional expectations.
Increasing replicas may reduce this gap while reducing independent viewing
samples can enlarge the confidence radius. The equal-noise-budget experiment
checks that tradeoff directly. A future multilevel variant would need independent
calibration, a finite-sample confidence allocation over levels and a fair
comparison with the classical nested-risk methods. No such algorithm, novel
complexity theorem, or favorable result is claimed here.
