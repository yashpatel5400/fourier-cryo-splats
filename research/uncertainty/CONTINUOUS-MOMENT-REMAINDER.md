# Integrating the nonlinear remainder over the support

Development derivation, 29 September 2026. This refines the continuous pose
audit after the wider-radius experiment exposed a dominant cubic remainder.
It does not change the locked continuous-v1 implementation or protocol.

## Proposition and proof

Let D=[-1/2,1/2]^3 with unit-volume Lebesgue measure. For any vector v and
p=2,4,6,

||v dot x||_Lp(D) <= c_p ||v||_2,
c_2=(1/12)^(1/2), c_4=(3/12^2)^(1/4), c_6=(15/12^3)^(1/6).

To prove this, expand the even moment of a linear form in independent centered
uniform coordinates. Terms with an odd exponent vanish. Every remaining term
has a nonnegative coefficient and uses moments E[X^2]=1/12, E[X^4]=1/80 and
E[X^6]=1/448, each bounded by the respective moment of N(0,1/12). Comparison
term by term with independent Gaussian coordinates gives the displayed bound.
The uniform coordinates are an integration device for the known operator
columns. No stochastic prior is placed on the unknown density.

For the pose path k(t)=k exp(t a [u_rot]_cross), with ||u||<=1, and a detector
shift t s q dot u_shift, write the Fourier phase as phi(t,x). Rotation preserves
norm, so ||k^(j)(t)|| <= a^j ||k||, j=1,2,3. The translation is linear in t.
Uniformly along the path,

||phi'||_Lp <= L_p = 2 pi (a ||k|| c_p + s ||q||),
||phi''||_Lp <= H_p = 2 pi a^2 ||k|| c_p,
||phi'''||_L2 <= T_2 = 2 pi a^3 ||k|| c_2.

The exact third derivative of exp(i phi) is
exp(i phi) [i phi''' - 3 phi' phi'' - i (phi')^3].
Minkowski's inequality and Holder's inequality imply the integrated bound

||d^3 exp(i phi)/dt^3||_L2 <= D(k,q) = L_6^3 + 3 L_4 H_4 + T_2.

Consequently the per-particle third-derivative operator is bounded by its
realified Hilbert--Schmidt norm

K_i = sqrt(sum_q |C_i(q)/sigma|^2 D(k_iq,q)^2).

For a specific estimator, its known complex weights
c_iq=(w_R,iq+i w_I,iq) C_i(q)/sigma permit the stronger scalar-field bound

|| sum_q Re[c_iq d^3 exp(i phi_iq)/dt^3] ||_L2
    <= sum_q |c_iq| D(k_iq,q) <= ||w_i|| K_i.

Taylor's integral remainder has weight (1-t)^2/2, integrating to 1/6. Thus
for P=||rho0||, ||rho-rho0||<=B, the entire scalar cubic remainder is at most

r_mom(w) = (B+P)/6 sum_iq |c_iq| D(k_iq,q).

Replace only the old cubic remainder in the continuous polynomial/spectral
audit with this expression (or retain the minimum of separately valid bounds).
The density residual and joint first/second field terms are unchanged. The
same conditional coverage theorem applies. This is an application of classical
moment comparison, Holder and Taylor inequalities to the Fourier-slice operator;
it is not a new general probability inequality.

## Verification and evaluation

Independent Gauss integration checks the even-moment envelope, exact third
derivatives computed through SciPy matrix exponentials, and exact nonlinear
Fourier remainders. This verifies implementation, not all admissible cases;
the proof supplies the uniform statement. Ordinary floating point remains an
implementation limitation.

First re-audit every existing development target at 0.1, 0.5, 1, 2 and 5
degrees, with no refitting, weight selection or new data-dependent tuning.
Preserve the older bounds and no-data fallbacks. Compare feasible coherent-pose
and separately searched nonlinear biases with both bounds. This is exploratory
method development after observing the older procedure's conservatism.
