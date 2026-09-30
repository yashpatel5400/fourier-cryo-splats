# A bounded sign/total-norm sensitivity ablation

This is a development response to the request to quantify the cost of the broad
signed density class. It does not supply a physically calibrated density norm
or assert that processed cryo-EM maps are nonnegative. Keep the existing
pilot-ball study unchanged. The calculations below are standard Hilbert-space
support-function arguments; the computational point is a continuous certificate
from exact cell averages rather than imposing a voxel representation on density.

## Declared alternative class

For a real density f on the unit cube, let f_- = max(-f,0). Consider

    C(R,eta) = { f : ||f||_2 <= R, ||f_-||_2 <= eta }, 0 <= eta <= R.

Eta=0 imposes nonnegativity; eta=R is the origin-centered L2 ball. This differs
from radius B around a nonzero pilot. It is not valid to merely replace B+P
by R in an old interval while retaining its old affine center or fallback.
The appropriate estimator below is w^T y + c with c recomputed for this class.

For h=ell-A_0^*w, write a=||h_+||, b=||h_-||. Its support function is

    S_eta(a,b) = R sqrt(a^2+b^2)                     if R b <= eta sqrt(a^2+b^2),
                 sqrt(R^2-eta^2) a + eta b         otherwise.

Proof: optimal positive and negative density components align with h_+ and
h_- on their disjoint supports. Their norms u,v satisfy u^2+v^2<=R^2 and
v<=eta. Maximizing au+bv over this quarter-disk gives the displayed formula.
The lower endpoint of <h,f> is -S_eta(b,a). This is an elementary specialization
of constrained inverse-inference support functions, not a new coverage principle.

## Continuous bounds using an orthogonal cell projection

Let P_m be conditional expectation onto orthonormal constant-cell functions,
whose coefficients of h are available exactly in real arithmetic from the
existing target integrals and Fourier cell adjoint. Jensen gives

    ||(P_m h)_-|| <= ||h_-||,    ||(P_m h)_+|| <= ||h_+||.

Suppose H>=||h|| and b0<=||h_-||. Maximize the preceding support function over
a^2+b^2<=H^2 and b>=b0, discarding other information. The result is

    U = R H                                      if R b0 <= eta H,
        sqrt(R^2-eta^2) sqrt(H^2-b0^2)+eta b0     otherwise.

For the second branch the objective decreases as b exceeds eta H/R, so its
maximum occurs at b=b0. Obtain an upper bound L on S_eta(-h) by exchanging
positive and negative projections. If the numerical projection has certified
L2 error e, max(0,||(computed projection)_-||-e) is a valid b0 because the
negative-part map is nonexpansive. The analogous bound applies to h_+.
As nested projections converge in L2, e vanishes, and H converges to ||h||, both bounds converge to
the exact continuous support endpoints. The unknown f is never made voxelwise
constant. In practical floating-point probes, heuristic pads and NUFFT checks
do not become rigorous interval arithmetic; state that limitation explicitly.

For a uniform nuisance field bound ||(A_u-A_0)^*w||<=d, both endpoints enlarge
by R d. Thus c=(U-L)/2 and bias b=(U+L)/2+R d give the original folded-normal
coverage statement for w^T y+c on this alternative class, conditional on its
other assumptions. The same spectral failure event can be reused. For a
nonnegative Gaussian target ell, the no-data interval is
[-eta ||ell||, R ||ell||]; its midpoint and half-width must both replace the
fitted interval on fallback. This construction does not assert efficient weights.

## Bounded initial probe, specified before its outcomes

Use the completed 10049 central target of sigma 10 A, 1,024 particles, radius 12
and the existing one-degree / 0.5 A enclosing-cube pose audit. Keep its fixed
weights. Evaluate cell projections at m=24,48,96, with all combinations
R=1,2,3 and eta/R=0,0.1,0.25,0.5,1. Report fixed-pose and pose-sensitive
intervals, all reference sign powers, and the reference's norm and negative-part
norm. The reference was normalized to unit L2 by the original simulation; this
does not estimate the physical norm of an experimental density. Label every
reference outside a sign constraint as outside-class, never as a coverage test.
No new acquisition data or pose prior is introduced. A failure at one degree
is retained as evidence that restricting signs alone did not fix the remainder.

## Completed initial probe

All specified projections/classes complete, with the original proposal and
source bytes archived in the result's source snapshot. Eight independent tests
include a conic-program support check and a continuous linear field whose
projection bounds converge from above. The normalized 64-cell reference has
negative-part norm 0.457856: R=1 with eta=0,0.1,0.25 excludes this processed
reference. This does not establish negativity of the underlying physical density.

At m=96, R=1 and eta=0.5, the nominal support half-width is 0.10609 instead
of the unrestricted-sign 0.10902. Both have fixed-pose sign power one for the
reference. The inherited one-degree nuisance field bound is 9.16209 per unit
density norm, so every one-degree case still has zero sign power, including
all cases whose class contains the reference. At R=1/eta=0.5 the pose interval
has relative half-width 0.73850 against its own narrower no-data class. This
is a same-weight negative diagnostic, not optimization over the sign class or
an impossibility theorem. It shows why merely restricting signs and total norm
does not fix this already-computed large pose remainder.
