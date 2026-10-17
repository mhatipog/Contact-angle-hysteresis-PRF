# Phase-sensitive wetting and first-passage lifetimes — controlled pilot

Original motivation: Hatipogullari et al. 2019 static metastability and the
2021 Hatipogullari–Ondarcuhu–Colinet dynamics draft (unpublished); Bradley,
Hewitt & Vella 2021 elastocapillary hybrid PDE.

## Rigorous no-go property of the *original* uniform-angle model

On the pinned-left / advancing-right branch, original fout.m is independent
of hysteresis lambda. For identical initial state, prescribed nu/V and
contact-line release threshold, ODE uniqueness forces two different
microscopic surfaces to have identical **first** contact-line event times
within that model. Any proposed pattern-phase effect must add physically
motivated information that original model discards.

## Exactly traceable extension without altering original solver

In the original fout.m, right receding mode [0,1] uses a ghost pressure
multiplier (1+lambda). Right advancing mode [0,2] uses multiplier 1.
For the conditional branch we can therefore call the **unchanged** function

    fout(t,y,N,nu,[0,1],A*g(xplus))

to implement right boundary pressure -nu*(1+A*g(xplus))/h_right, even though
the right line is *physically advancing*. This is **algebraic reuse** of the
original discretization; no claim about contact-line regime [0,1] is made.
No right contact-line reversal is allowed in the time window studied.

The constitutive addition is kappa(x)=1+A*g(x) for the *advancing*
right meniscus, i.e. a spatially varying advancing pressure/angle.
This is a controlled candidate model, **not a derived homogenization of
the 2019 microstructure** and not experimentally validated. Real surfaces
may require microscopic contact-line mobility, consistent advancing/receding
thresholds and nonlocal energetic forces. In particular, a full physical
extension would need an energy/dissipation audit.

For two *opposite-phase* patterns g(x)=+sin(2pi(x-xminus)/ell) and -g(x),
the minimum, maximum, average and variance of the global pressure multiplier
are identical, and both patterns match at the pinned left location. However,
the advancing right meniscus samples different wetting histories.

## Sensitivity theorem (conditional branch)

Let y_A(t) solve y'=F(y,A) from the same pinned initial state.
Let Gamma(y) be the left pressure-yield observable and define
t_*(A) through Gamma(y_A(t_*(A)))=lambda+delta*(1+lambda).
For transverse crossings,

    dt_*/dA = - [d Gamma(y_A(t))/dA]_{A=0,t=t_*}
                        / [d Gamma(y_0(t))/dt]_{t=t_*}.

The numerator can be found by the tangent equation

    z'=J_0(t) z + [partial_A F]_{A=0}(y_0(t)),   z(tpin)=0,

or by centered differences on a pair of small-amplitude trajectories.
The formula is implicit-function calculus, not a new mathematical theorem.
Its useful wetting content lies in how microscopic phase changes observable
release time despite identical global extrema and local left threshold.

The benchmark uses the author's original unmodified PDE for N=12 and 16,
first-pin state at nu=4,V=.2,xplus0=.65, and two release thresholds .02/.028.
It checks first-order crossing shifts against finite-amplitude trajectories.

**Limitations**: initial pinned state is PREPARED using unpatterned advancing
boundary; pattern is applied after first pinning. This is a controlled
conditional-state experiment, not a complete patterned droplet's initial
history or full escape/trapping classification. Numerical tolerance=0.001
belongs to original event logic and is not physical hysteresis.
No upstream PR is created.

## Independent verified beam threshold and inverse sensitivity

The exact quartic clamped/free beam boundary conditions admit an implicit
function theorem. With q=(xplus,p,c,d), R(q,A)=0 and
Gamma(q)=-p*h_left/nu-1, differentiating yields

    dq/dA=-R_q^{-1} R_A,   dGamma/dA=Gamma_q*(dq/dA).

The complex-step 4x4 Jacobian gives slope dGamma/dA=-0.896031545419
at nu=4,V=.2,xminus=.4454383377, period=.25. For pinning threshold
lambda=.0293, the LINEAR predicted positive-phase transition A=0.00161456166;
the exact nonlinear equilibrium gives A=0.00161419996, a 0.0224% relative
error. The original event safeguard of 0.001 changes the apparent numerical
switch amplitude to A~=0.00046580; these are different quantities.

A=+.002 yields stationary Gamma=0.02895414 < lambda; A=-.002 yields
Gamma=0.03253826 > lambda. Original source discretization plus physically
patterned right boundary confirms the first-event classification through
observation t=0.1 across N=12,16,20. This is a conditioned, model-level
phase contrast, **not** validated experimental hysteresis and **not**
a demonstrated universal novelty result.
