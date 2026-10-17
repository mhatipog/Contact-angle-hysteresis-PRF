# Two patterned menisci from time zero — physics and scope

Builds on Hatipogullari et al. (2019) and the unpublished dynamics draft,
plus the unmodified Bradley-Hewitt-Vella (2021) flexural hydrodynamics.

On the branch where left and right menisci advance, pressure multipliers are
kappa_i=1+A*sin(2*pi*(xi-xref)/ell)-mu*s_i*v_i, with s_left=-1, s_right=+1.
For left-pinned/right-advancing branch only right multiplier matters.
The original RHS is affine in both multipliers:
F=F0+B*(kappa-1), B_i=F(kappa_i=2)-F(kappa_i=1).
Therefore the implicit velocity-dependent closure is EXACTLY
v=(I+mu*B_vel*diag(-1,+1))^{-1}*(v0+B_vel*A*g).
This only solves an algebraic 2x2 system on top of the unchanged PDE
discretization. A=0,mu=0 gives the original wetting trajectory exactly.

Local allowable pressure for left pinning is kappa_A(xl)<=kappa<=
kappa_A(xl)+lambda; the numerical guard adds tol*(1+lambda). For a
physically possible cosine interpretation choose reference cos(theta_a)=.94,
so cos(theta_i,A)=.94*kappa_A and cos(theta_i,R)=.94*(kappa_A+lambda),
both in (0,1) over the tested patterns. The additive hysteresis in
pressure multiplier is a *model assumption*; not deduced from the
microstructural free energy in the user's statics theory.

Linearized kinetic term is dissipative in the moving direction:
(c_Y-c_app)*(s_i*v_i)=cos(theta_ref)*mu*v_i^2>=0, with s_left=-1 and
s_right=+1. This is ONLY the contact-line contribution, NOT the full
fluid-beam energy audit. Bulk viscous dissipation and boundary meniscus
energy must be added before claiming the full PDE is thermodynamically
closed. Also, the original first-event label [2,2]->[0,2] is assumed to
remain reachable; the test reports a failure otherwise.

The test integrates BOTH patterned advancing menisci from t=0, detects
their first velocity event, then tracks the left-pinned branch, compares
its yield observable against the independently solved exact quartic
steady beam, and checks event detection and mass conservation at
N=12,16,20 and mu=0,.01. No permanent trapping or experimental
calibration is claimed. No upstream PR is affected.