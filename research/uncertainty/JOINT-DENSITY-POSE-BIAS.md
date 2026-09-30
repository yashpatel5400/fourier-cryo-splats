# A residual/pose cross-term post-audit

This development is an elementary Hilbert-space refinement of the existing
triangle bound. It is not a claim of a novel general inequality, nor a
calibration result for experimental nuisance assumptions.

Let h=ell-A_0*w and D=A_t*w-A_0*w. For rho=rho_0+u with ||u||<=B and
||rho_0||<=P, the scalar bias equals <u,h-D>-<rho_0,D>. Write the existing
second-order pose expansion D=F v+e, where fixed positive block scales imply
||v||<=L=sqrt(sum d), and ||e||<=r. Suppose f bounds ||Fv|| uniformly and
c bounds ||F*h||. Then, for every allowed pose and density,

  |bias| <= B sqrt(||h||² + 2 L c + f²) + P f + (B+P)r.

Proof: expand ||h-Fv||² and use |<F*h,v>|<=L c and ||Fv||<=f. Apply
Cauchy--Schwarz separately to u and rho_0, then charge the same remainder e
to their sum. The minimum with the old bound B||h||+(B+P)(f+r) is valid.
This uses the SAME spectral event already bounding f for fixed weights; it
does not require another numerical failure budget. No claim is made that
this upper expression is convex as a function of estimator weights.

## Continuous target cross products

The operator F has degree-at-most-two polynomial Fourier columns. Gaussian
target moments on the cube are analytic and include finite-boundary terms.
For a normalized one-dimensional Gaussian g with mean mu and variance s²,
let I_j=int[-1/2,1/2] x^j g(x) exp(i omega x) dx. Integration by parts gives

  I_1=(mu+i s² omega) I_0 - s² [g(x)exp(i omega x)]boundary,
  I_2=(mu+i s² omega) I_1 + s² I_0
        - s² [x g(x)exp(i omega x)]boundary.

Products of these moments evaluate <ell,F_j> without replacing the unknown
density by a grid. For <A_0*w,F_j>, use the existing polynomial Fourier Gauss
quadrature. If T is the absolute Fourier coefficient sum of A_0*w, m_j the
scaled polynomial coefficient bound of F_j, and E the existing kernel error,
then the Euclidean error of the full cross vector is at most E T ||m||.
The implementation adds this error before using c. Ordinary special-function,
NUFFT and floating-point errors remain numerical checks, as in the original
method; they are not certified by interval arithmetic.

Independent tests compare the moments with high-order one-dimensional
quadrature including a Gaussian near the cube boundary, compare analytic
pose target pairings with direct field integration, and check the bound on
independent Hilbert-space instances. A case with orthogonal residual and
pose range illustrates a strict improvement over the triangle bound. Archive
experiments must still establish whether this refinement is useful there.

## Optional sharp cube directional moments

For X uniform on the cube, ||u||=1 and p_j=u_j², direct independence gives

  E(u.X)^4 = 1/48 - (sum p_j²)/120,
  E(u.X)^6 = 5/576 - (sum p_j²)/96 + (sum p_j³)/252.

The fourth moment is maximal at p=(1/3,1/3,1/3). For the sixth, the difference
of partial derivatives divided by p_i-p_j is -1/48+(p_i+p_j)/84<0 on the
simplex. Pairwise averaging therefore increases the expression, so the same
point maximizes it. The maxima are 13/720 and 205/36288, respectively.
The second moment is exactly 1/12. These replace the preceding Gaussian
moment majorants 3/12² and 15/12³ in the already integrated cubic remainder.
All Minkowski/Holder steps and the nonlinear rotation/translation path stay
unchanged. This elementary constant refinement is not claimed as a new
general moment theorem. The `--sharp-cubic` post-audit stores separate output
and preserves both old and new remainder values.

Independent degree-exact cube quadrature verifies the moments and their
attainment. Further random single-frequency tests evaluate the exact nonlinear
pose field and its independently assembled second-order Taylor polynomial.
A NumPy scalar-conversion error in that test fixture was corrected before
the archive audit; the failed and successful test logs remain distinct.
