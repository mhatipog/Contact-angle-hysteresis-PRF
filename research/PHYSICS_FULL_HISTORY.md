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
## Nonlinear contact-line closure directly from user's 2021 dynamics law

Use the *unlinearized* positive-angle cubic relation
s_i*v_i=M*(theta_app,i^3-theta_Young,i^3), with left s=-1, right s=+1.
For Young cosine c_Y(x)=0.94*(1+A*g(x)), set

   theta_Young=acos(c_Y),
   theta_app=(theta_Young^3 + s*v/M)^(1/3),
   pressure_multiplier=cos(theta_app)/0.94.

At each RHS evaluation, the original PDE's affine pressure map gives a
2-by-2 Newton solve for contact-line velocities. The right-hand-side
solver, pressure ghost-point formulas and beam equations stay unchanged.
The physical microscopic contact-line dissipation is nonnegative:

   D_CL=sum[(c_Y-c_app)*(s_i*v_i)] >= 0

because x -> cos x decreases on (0,pi) and the cubic relation is monotone.
This checks **contact-line dissipation only**, not the complete fluid-beam
energy budget (bulk viscous dissipation and mechanical work remain to be
explicitly audited). The local receding angle is an independent hysteresis
offset from advancing pressure, which is an effective model assumption,
not derived from resolved defect-scale stick-slip energy. Tests run from
initial patterned state through first release, N=12/16/20, M=40/100.

## Exact conservative stationary energy: thermodynamic connection to CAH

The fixed-volume quartic Euler-Bernoulli beam has an explicit compliance
I(a,b)=(b-a)^2*((b-a)^3+40*((a+b)/2)^3)/120.
Its reduced interfacial+elastic energy for Young pressure multiplier
kappa_Y(x)=1+A sin(k(x-xref)) is

  E(a,b)=(V-(b-a))^2/(2I(a,b)) - nu * integral_a^b kappa_Y(x)dx.

At each exact stationary quartic beam with free right contact line,
complex-step differentiation independently verifies

  dE/db = 0;
  dE/da = -nu*[Gamma - A sin(k*(a-xref))],

where Gamma=-p*h_left/nu-1 is precisely the pressure-based first-depinning
quantity. The local receding criterion Gamma-loc>=lambda is thus a
rate-independent *yield condition* expressed by this generalized energy
gradient. 25 independent parameter cases agree at <4e-13. A reduced
gradient inclusion with dry friction Y and viscous drag zeta,

  0 in grad E + Y*partial |v| + zeta*v,

satisfies the exact reduced energy identity
  dE/dt=-Y*sum|v|-zeta*sum(v_i^2)<=0.
The audit tests this identity on 1000 additional arbitrary states. It
does not establish the COMPLETE time-dependent fluid-beam energy law.

## Caveat about applying the full angle-cubed relation at startup

The author's original initial state produces O(40-70) dimensionless
initial contact-line speeds, well outside the small-speed quasi-steady
range when interpreting the cubic Cox law at the tested M=40/100.
The full two-line *linearized* kinetic closure remains in use during
that initial stiff transient, retaining each surface's own history.
The unlinearized angle-cubed law is applied **after the first left pin**
when the advancing-right speed is small, with M=100 and 300.
Do not represent this as all-time full Cox integration. At t=0 the
cubic law can cross its physically admissible angle domain during
stiff solver trials; silently clipping contact angles is prohibited.
