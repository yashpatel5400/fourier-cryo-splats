# Global rotation-envelope cost gate: complete outcomes

The declared curvature calculation completes for all eight nonzero directions; four zero directions are retained as skipped. No rotation grid was evaluated. The analytic covering is sufficient, not optimal, and its size is not a lower bound on certification complexity. Ordinary floating-point constants are not outward-rounded numerical certificates.

| Stack | Contrast | Curvature bound | Existing sampled margin | Proposed grid points |
|---|---|---:|---:|---:|
| 10028 | power_fixed | 1.30459e+06 | 0.14652 | 11122470703127 |
| 10028 | power_range09_11 | 2.10663e+06 | -0.00369244 | already invalid |
| 10028 | power_bispectrum_fixed | 1.34188e+08 | 0.102879 | 19718313248600912 |
| 10028 | power_bispectrum_range09_11 | 2.11478e+08 | -0.0410645 | already invalid |
| 10049 | power_fixed | 536067 | 0.0436558 | 18014276769435 |
| 10049 | power_range09_11 | 562923 | 0.017434 | 76805868171475 |
| 10049 | power_bispectrum_fixed | 1.02836e+07 | 0.0507924 | 1205927247450744 |
| 10049 | power_bispectrum_range09_11 | 1.56889e+07 | 0.0296193 | 5102943272005897 |

All six positive-margin constructions exceed the 10-million-point gate, by more than six orders of magnitude. This particular absolute-density derivative envelope is computationally impractical. Its failure does not establish that a tighter continuous certificate is impossible. The two 10028 ranged-amplitude cases already have counterexamples and receive no proposed grid. No GPU rental or enumeration follows.

The separate Monte Carlo protocol changes the model: it assumes a known bound on the viewing density relative to Haar measure. It cannot replace or be described as the missing arbitrary-view result.
