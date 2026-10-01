# Sharp boundary of the normal bias/variance envelope

1 October 2026 UTC. Elementary supplementary derivation; no priority claim and no new cryo-EM coverage claim.

For z>0 define F_z(t)=Phi(z sqrt(1+t^2)-t)+Phi(z sqrt(1+t^2)+t)-1. The existing proof shows F_z(t)>=F_z(0) for all t>=0 when z>=sqrt(3). This sufficient threshold is also necessary. Taylor expansion at zero gives

    F_z(t) = 2 Phi(z)-1 + phi(z) z (z^2-3) t^4 / 6 + O(t^6).

For 0<z<sqrt(3), the coefficient is strictly negative. Thus arbitrarily small fixed nonzero bias violates the proposed nominal coverage. The positive-z condition excludes the trivial zero-coverage case z=0. This sharpens the range of the existing proposition; it does not make the underlying optimal-recovery or prior comparison new.

There is a different sharp statement when the bounded bias can adapt to the realized Gaussian noise but the projected noise remains fixed N(0,s^2). Write Z~N(0,1), |beta(Z)|<=b, with s=1 by scaling. For a proposed half-width q, the adversarial choice beta(Z)=b sign(Z) maximizes absolute error pointwise, giving worst coverage 2 Phi(max(q-b,0))-1. For q=c sqrt(1+b^2), c>1, the minimum over b>=0 is attained at b=1/sqrt(c^2-1) and equals

    2 Phi(sqrt(c^2-1))-1.

For 0<c<=1 the infimum is zero. To retain nominal 2 Phi(z)-1 coverage uniformly over every b and adaptive bounded bias, an RMS-form constant must therefore satisfy c>=sqrt(1+z^2). At z=1.95996398454, using c=z instead has worst coverage about .908, despite valid fixed-bias coverage. The usual b+z s bound remains pointwise no wider than sqrt(1+z^2) sqrt(s^2+b^2), by Cauchy–Schwarz.

This analysis does not handle data-dependent weights, an estimated variance, or non-Gaussian processed noise. It isolates why a fixed-bias calculation cannot be transferred to an estimator with a noise-adaptive error term. The local-pose study still needs empirical validation; no new interval is selected from its outcomes.
