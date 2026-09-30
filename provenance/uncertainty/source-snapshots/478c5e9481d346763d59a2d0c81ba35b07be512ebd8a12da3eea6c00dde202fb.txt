# Per-frequency Taylor disks: a local computational refinement

30 September 2026. Motivated by Change 1 in the unmodified focused
`mixture-audit-01` Fable review. This is classical Fenchel duality specialized
to the already derived Taylor errors, not a new statistical theorem.

For one particle let r be the center residual, J the mean Jacobian, and v a
three-dimensional Euler displacement with |v_j| <= h_j. The complex Taylor
error in frequency q has modulus at most

    eps_q = H^2/2 |C_q| (K_q^2 B2_q + K_q B1_q),  H = sum_j h_j.

Thus the realified error consists of two-dimensional blocks e_q in disks
||e_q|| <= eps_q. These are disks, not independent real/imaginary intervals:
repeating eps_q in a coordinate rectangle enlarges the aggregate error radius
by up to sqrt(2), so it would not support an unqualified dominance claim over
the earlier image-wide ball. We preserve the complex block structure.

For every realified dual vector u, weak Fenchel duality gives

    D(u) = 2 u.r - ||u||^2 - 2 sum_q eps_q ||u_q||
           - 2 sum_j h_j |(J^T u)_j|
      <= min_{|v|<=h, ||e_q||<=eps_q} ||r-Jv-e||^2.

Proof: use ||x||^2 >= 2u.x-||u||^2 and minimize the linear terms separately
over the Euler box and error disks. Every nonlinear residual is feasible in
this relaxed minimization. Consequently -max(D(u),0)/2 bounds its Gaussian
log kernel from above. The Gaussian normalization is unchanged.

Eliminating e gives the convex primal

    min_{|v|<=h} sum_q max(||(r-Jv)_q||-eps_q,0)^2.

At each feasible v, choose u_q = (1-eps_q/||(r-Jv)_q||)_+ (r-Jv)_q, with zero
at zero norm. This supplies a valid dual point and a feasible primal value.
The implementation retains all lower dual values and all upper primal values.
Eight ordinary least-squares coordinate sweeps provide starting anchors,
followed by 24 projected-gradient steps on the disk primal with step size
1/lambda_max(J^T J), interpreted as a guide in ordinary floating point.
No convergence or eigenvalue certificate is needed for dual validity; a bad
step affects only the chosen anchor. Unit tests compare with a separately
constructed conic primal, replay dual vectors, check an exact allocation
example and verify nonlinear feasible orientations.

At the exact optimization level, the product of frequency disks is contained
in the image-wide ball of radius ||eps||. Its residual minimum is therefore
at least the ball minimum, which is max(min_v ||r-Jv||-||eps||,0)^2.
Finite dual iterates need not achieve that comparison. We explicitly retain
the minimum of the earlier likelihood envelope and the new disk envelope,
ensuring that adding this computation cannot worsen the recorded upper bound
in real arithmetic. Small-box curvature remains included by the earlier
five-degree rule. All cells retain the same meaning and full-cover requirement.

This module also rejects DC, duplicate and conjugate-duplicate plane frequencies.
That guard does not establish noise independence, a valid whitening estimate,
CTF independence, an experimental viewing distribution or physical homogeneity.
No result here removes those statistical assumptions or validates rounding.
