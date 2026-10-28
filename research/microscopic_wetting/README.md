# Microscopic-energy folds: exact reduced beam, two moving contact lines

**Research pilot, 28 October 2026.** Separate from the original Bradley–Hewitt–Vella (2021) MATLAB code and from the user's unpublished microchannel manuscript. It is not a new full lubrication solver or a validated experimental prediction.

## Aim and important distinction

Earlier branches imposed a local advancing/receding interval and demonstrated a patterned-surface phase effect in a hybrid lubrication/beam PDE. This pilot instead specifies **one conservative chemical surface-energy landscape**, and asks which *stationary* contact-line configurations become metastable and where they cease to exist. There is no externally prescribed hysteresis width.

For left/right meniscus positions \(0<a<b<1\) and nondimensional liquid volume \(0<V<b-a\), the clamped–free beam Green function is

\[
G(x,s)=\frac{\min(x,s)^2(3\max(x,s)-\min(x,s))}{6}.
\]

Define \(I=\int_a^b\int_a^b G(x,s)\,ds\,dx\). This integrates *exactly* to

\[
I(a,b)=\frac{(b-a)^5}{120}+\frac{(b-a)^2(a+b)^3}{24}.
\]

For \(\kappa_Y(x)=\kappa_0+A\sin(2\pi x/\ell+\phi)\), define the reduced conservative beam/Young energy

\[
E(a,b;V)=\frac{[V-(b-a)]^2}{2I(a,b)}
-\nu\int_a^b\kappa_Y(x)\,dx .
\]

The Young cosine itself is \(\cos\theta_Y(x)=c_{\rm ref}\kappa_Y(x)\); the default \(c_{\rm ref}=0.9\), \(A=0.1\) ensures \(0.81\leq\cos\theta_Y\leq0.99\). No arbitrary angle clipping is used. The pressure conjugate to volume is exactly \(p=\partial E/\partial V=[V-(b-a)]/I\). We check positive gap at the beam's free end, which is the extremal deflection for the suction branch considered here.

The equilibrium condition is \(\nabla_{a,b}E=0\). The coupled fold condition adds \(\det \nabla^2_{a,b}E=0\), and the second Hessian eigenvalue must be positive. This is a *local* metastability criterion, not a proof of global trapping or a fully dynamical release time.

## Concrete reproducible result

Default dimensionless inputs: \(\nu=4,\ \kappa_0=1,\ A=0.1,\ \ell=0.075,\ \phi=0\).

The two-meniscus stationary fold solves to

| Quantity | Result |
| --- | ---: |
| Left position \(a_c\) | 0.485098866979384 |
| Right position \(b_c\) | 0.715416013259478 |
| Critical volume \(V_c\) | 0.213851157109733 |
| Pressure \(p_c\) | -4.299629075650792 |
| Hessian eigenvalues | approximately \(0,\ 380.868699\) |
| Soft eigenvector (sign arbitrary) | (-0.730407, -0.683012) |

At \(V=0.214\), numerical stationarity finds a nearby local minimum at approximately \((a,b)=(0.482107,0.712778)\) and a saddle at approximately \((0.488052,0.718348)\). Their smallest Hessian eigenvalues are \(+1.955\) and \(-1.785\), respectively. Thus a genuine minimum/saddle pair is created by the coupled beam/chemical energy near the fold, rather than by a prescribed local hysteresis interval.

**Why two-line coupling matters:** at the stable configuration, the apparent left-meniscus stiffness at fixed \(b\) is \(H_{aa}\approx182.29\). Letting the right line relax gives the *exact* static Schur complement

\[
K_{\mathrm{left,relaxed}}=H_{aa}-H_{ab}^2/H_{bb}\approx3.67390 .
\]

That is a **49.6-fold** curvature overestimate if one incorrectly freezes the right contact line. The *beam deformation itself* has already been eliminated through the exact Green function; this Schur complement only eliminates the other contact-line coordinate.

This fold is a dimensionless numerical existence result for the specified hypothetical pattern, **not** a measured threshold, universal constant, complete history-dependent hysteresis loop, or original mathematical theorem.

## Thermodynamics

For an *illustrative* positive diagonal mobility \(M\) and an externally varied volume \(V(t)\), consider \(\dot{\mathbf x}=-M\nabla E\). The reduced identity is exact:

\[
\frac{dE}{dt}=\frac{\partial E}{\partial V}\dot V
-\nabla E^\top M\nabla E .
\]

The last term is nonnegative dissipation. This is **not yet** the full Bradley–Hewitt–Vella fluid/beam dissipation balance, which needs the actual hydrodynamic mobility, bulk viscous dissipation, capillary boundary work, and complete pin/depin history.

## Reproduce

From this repository's root:

    python -m pip install numpy scipy
    python research/microscopic_wetting/beam_fold.py
    python -m unittest discover -s research/microscopic_wetting -p 'test_*.py' -v

Six core tests cover: direct numerical quadrature of the Green kernel; independent complex-step energy gradient and finite-difference Hessian; physical Young cosine range; fold plus adjacent minimum/saddle; volume-work and dissipative energy identity; invalid/collapsed configurations.

## Publication gates

1. Derive **microscopic** periodic surface energy from measured geometry/chemistry, not an arbitrarily selected \(A,\ell,\phi\); include fabrication error and wetting ageing.
2. Follow an actual loading/unloading cycle, distinguishing metastability from rate-dependent delay and obtaining both contact-line histories. Verify convergence in time, spatial pattern resolution, and energy.
3. Replace the illustrative diagonal mobility with the original coupled lubrication/beam mobility. Confirm the full energy balance and handle re-pinning.
4. Blind-validate predicted event positions and critical volumes against repeated observations on patterned and unpatterned channels, with uncertainty bars.
5. Establish novelty against Bradley–Hewitt–Vella (PRF 2021), Chen et al. (Capillarity 2025, doi:10.46690/capi.2025.04.02), Millet et al. (IJNAMG 2026, doi:10.1002/nag.70363), and Zhang et al. (JCIS 2026, doi:10.1016/j.jcis.2025.139781). In particular *a variational hysteresis framework or pattern-dependent trapping alone is not novel*.

**Publication status:** a properly checked numerical intermediate result, not journal-ready. Full experimental or dynamical validation remains the decisive next step.

## Inverse-design / uncertainty sensitivity (local fold branch)

The fold equation is \(F(a,b,V;A,\phi)=(E_a,E_b,\det E'')=0\).
Provided this three-equation Jacobian remains nonsingular, implicit differentiation gives

\[
\frac{d(a,b,V)}{d\theta}=-F_{(a,b,V)}^{-1}F_\theta,\qquad \theta\in\{A,\phi\}.
\]

The new 'sensitivity.py' evaluates this on the **same** identified fold branch. At the reference point:

- \(dV_c/dA=-0.0482989900\) per unit Young-multiplier amplitude.
- \(dV_c/d\phi=+0.0008219084\) per radian of surface phase.
- These agree with independent centered re-solving of the fold to better than \(10^{-5}\) on the tested perturbations.
- Nearby folds also survive \(A=0.04,0.06,0.08,0.10\), although a *global* phase diagram has not been constructed.

This is a local uncertainty/design calculation, not an extrapolation to arbitrary large pattern shifts or uncontrolled real surfaces. The two added tests bring the local suite to **8 tests**.
