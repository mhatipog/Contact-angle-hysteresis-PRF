# First-passage yield envelope for contact-angle hysteresis

This branch is a **science-only**, independent investigation based on the
original source of Bradley, Hewitt & Vella (2021), and on static metastability
and driving-protocol ideas of Hatipogullari et al. (2019) and Hatipogullari,
Ondarcuhu & Colinet (2021 unpublished dynamics draft).

**Key derivation.** After the first left moving-to-pinned transition and
while the right meniscus remains advancing, the original PDE vector field
`fout(t,y,N,nu,[0,2],lambda)` is *exactly independent of lambda*. Let
`w=x_plus-x_minus`, `U=h*w`, and `p_minus` the pressure diagnostic in
the original `EventsFcn`. Set

    Gamma(t) = -p_minus(t)*U_minus(t)/(nu*w(t)) - 1.

With the advancing-angle normalization `cos(theta_a)=1/(1+lambda)`,
the original receding event evaluates
`cos(theta_left)=(1+Gamma)/(1+lambda)`.
Because the original code uses tolerance `tol=0.001`, the first receding
event after pinning has the first-passage condition

    Gamma(t) = lambda + 0.001*(1+lambda),

if the crossing is in the triggering direction. Without that numerical
tolerance the physical ideal condition is `Gamma(t)=lambda`. In this
regime *one* pinned PDE integration gives the first depinning times for a
continuum of lambda values, via a first-passage inverse, and the sufficient
no-depinning inequality `lambda > sup_t Gamma(t)` (over the same horizon,
subject to all other regime switches remaining absent).

This is a structural property of a conditional hybrid dynamical branch,
not a newly discovered universal dynamical system theorem. It does NOT
classify eventual escape, which depends on the post-depinning trajectory,
nor is the steady lambda_e by itself a safe threshold if Gamma overshoots.
The stationary reference uses the independent four-equation beam BVP.
Numerical grid and time convergence are assessed in GitHub Actions.

Run using original unmodified Matlab/Octave physics:

    octave --quiet --no-gui --eval "run('research/threshold_from_trajectory.m')"

The repo's original source has its own MIT licence. No upstream PR is created.
