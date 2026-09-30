# Prospective Fourier-cancellation remainder probe

This is a next-step derivation to implement and independently check, not an
available result or a claimed interval improvement. The completed 1,024-particle
10 Å Gaussian-width audit has cubic bias 56.4927, polynomial pose bias 6.3760
and no-data fallback. Simply using the next Taylor order would still need a
large fourth-order remainder: an exploratory calculation gives 15.5022 for
that term alone with Gaussian eighth-moment and sharp fourth/sixth constants.
That number is NOT a valid replacement bias unless third-order fields are also
included. No such higher-order implementation has been completed.

## Candidate bound retaining frequency cancellations

For particle i define the real adjoint field

    f_i(x)=Re sum_q c_iq exp(2 pi i k_iq.x),
    c_iq=(w_R+ i w_I) C_iq / noise.

The nonlinear pose path is `f_i(R_i(t)x+t b_i)`, where the embedding-dual
translation vector satisfies `k_iq.b_i=q.tau_i`. For orthonormal detector
embeddings `||b_i||<=s`; with recorded approximate embeddings use a bound
`s_i=s/sigma_min(E_i)` and check the representation residual explicitly. Do
not silently assume float32 metadata rotations are exactly orthonormal.

On the unit cube put R=sqrt(3)/2. Every transformed point lies in the ball of
radius r_i=R+s_i, and the path derivatives satisfy

    ||y'||<=L=a R+s_i, ||y''||<=H=a²R, ||y'''||<=T=a³R.

The chain rule, triangle inequality and a rigid change of variables suggest

    ||g_i'''(t)||_L2(cube)
      <= L³ S_i3 + 3 L H S_i2 + T S_i1,

where `S_im=||D^m f_i||_{F,L2(ball(r_i))}` is a Sobolev seminorm of the KNOWN
finite Fourier field, not a smoothness assumption on the unknown density.
Taylor's integral remainder would then give scalar bias at most

    (B+P)/6 sum_i [L³ S_i3 + 3 L H S_i2 + T S_i1].

The enclosing ball and Frobenius norms may be conservative. Unlike the current
frequency-wise triangle bound, these S terms retain cross-frequency
cancellations. The old and new bounds could be compared on the same class;
no density-support restriction is introduced by integrating the known field
over a larger ball.

## Exact real-arithmetic integrals to verify

For `K_r(v)=int_ball(r) exp(2 pi i v.x) dx`, with z=2 pi r||v||,

    K_r(v)=4 pi r³ j_1(z)/z,
    K_r(0)=4 pi r³/3.

The candidate seminorm identity is

    S_m² = .5 (2 pi)^(2m) Re sum_q,l (k_q.k_l)^m
      [c_q conj(c_l) K_r(k_q-k_l) + (-1)^m c_q c_l K_r(k_q+k_l)].

The (-1)^m in the sum-frequency term comes from i^(2m). Compute in particle
blocks, with small-z series for the ball transform, explicit cancellation
diagnostics and no claim of validated floating-point bounds.

Before a real-data probe, compare S_1,S_2,S_3 against independent high-order
spherical quadrature of all ordered derivative tensors. Test the resulting
nonlinear remainder against direct cube integration on unrelated random
Fourier fields and allowed joint pose directions. Derive/check the embedding
dual translation carefully, including its numerical representation error.
Only then evaluate the saved high-band weights, retaining an unfavorable result.
If the enclosing-ball bound is still too large, it cannot be called a solution
to the fine-feature precision problem.

## Implementation and second prospective enclosure (2026-09-30)

The first ball probe is now complete. Independent spherical integration of all
ordered derivative tensors through order three and direct nonlinear cube
remainders pass. On the saved 1,024-particle case, choosing the smaller original
or ball bound for each particle reduces the cubic bias from 56.4927 to 32.6233;
the interval still falls back to no data. The original result and source snapshot
are retained under `ball-remainder-probe`. This is conditional on the same
known-noise model, not experimental calibration.

The implementation uses the joint-pose speed `L=hypot(a R,s_i)`, justified by
Cauchy--Schwarz on the allowed five-vector. It also bounds the small embedding
residual explicitly. Let `b=E^+ tau`, `delta=2 pi (q-k E^+) tau`, and M(t) be
the main transformed exponential. Differentiate `(exp(i t delta)-1)M(t)`
three times. With `e >= |delta|`, phase derivative bounds v,h,j, the third
norm is at most `e(v^3+3vh+j)+3e(v^2+h)+3e^2 v+e^3`. Summing coefficient
magnitudes and dividing by six gives the recorded residual pad.

Before running the next probe, the following tighter enclosure is specified.
For every point of the unit cube, `||(R(t)-I)x|| <= 2 R sin(min(a,pi)/2)`.
Consequently the same transformed cube lies inside the coordinate cube with
half-side `r=.5+2 R sin(min(a,pi)/2)+s_i`. Its Fourier integral is exactly
`K(v)=prod_j 2r sinc(2r v_j)`. The same derivative-tensor identity and chain-rule
argument apply, because a rigid change of variables has determinant one and
the integrand is nonnegative. No density-support assumption changes. The new
probe will use this enclosure and will again retain an unfavorable result.
Independent tensor-product quadrature now also checks all ordered derivative
norms on the expanded cube; direct nonlinear pose tests cover both enclosures
and deliberately distorted (nonplanar) frequency embeddings.

A useful bound follows from Taylor's integral formula:
`||g(1)-g(0)-g'(0)-g''(0)/2|| <= sup_t ||g'''(t)|| / 6`.
This establishes the real-arithmetic claim. Floating-point summation and
special functions only receive heuristic magnitude guards; this remains a
numerical certificate, not validated interval arithmetic.
