# Finite-wavelength wettability in elastic channels: a stationary pilot

**Scientific prototype; not a validated transport or trapping model.** This is
deliberately separate from the published Figure 11 reproducibility PR.

The goal is to connect the sinusoidal local Young-angle field and finite-wavelength
stick-slip physics of Hatipogullari et al., *Phys. Rev. Fluids* 4, 044008
(2019), to the stationary elastic beam of Bradley, Hewitt & Vella,
*Phys. Rev. Fluids* 6, 114003 (2021).

A directly relevant 2025 study already exists: Chen et al., *Capillarity*
15, 4–11, DOI: 10.46690/capi.2025.04.02, applies heterogeneous wettability
in a phase-field model of deformable nonparallel channels. Therefore
"patterned wetting on flexible channels" alone is **not a novelty claim**.

## Model and audit

Take c_Y(x) = c0 + A sin(2 pi x / ell + phi). Impose the local Young-angle
law **at both menisci**, rather than simultaneously imposing an unrelated
empirical advancing/receding threshold. With fixed volume V, pressure p
and channel half-gap h, the beam obeys h''''=p on [x_minus,x_plus], zero
otherwise, clamped at x=0 and free at x=1. The exact cantilever Green
kernel is G(x,s)=min(x,s)^2 (3 max(x,s)-min(x,s))/6, so
h(x)=1+p integral_[x_minus,x_plus] G(x,s) ds.

Volume gives p directly in terms of the unknown contact-line positions.
The two residual equations are
p h(x_minus)+nu c_Y(x_minus)/c_ref=0,
p h(x_plus)+nu c_Y(x_plus)/c_ref=0.
The code only accepts interior, noncontacting, wetting solutions with
small (<1e-8) pressure residual and positive gap.

At nu=4, V=.2, c0=.95, A=.025 and phase=0, the multistart search finds:
ell=.25: 2 candidates; ell=.15: 6; ell=.10: 0. The .15/.10 contrast
survives two numerical seed resolutions; ell=.10 is also zero for four
phases. **Counts are observed candidates, not exhaustive mathematical
counts**, and do not imply stability, accessibility, or trapping.

A possible geometric mechanism is that two menisci separated by d
sample equal chemical phases when d/ell is approximately an integer.
For sinusoidal wettability their local cosine difference equals
2 A sin(pi*(x_minus-x_plus)/ell)
    cos(pi*(x_minus+x_plus)/ell+phi).
The local contrast can cancel despite identical mean and amplitude.
This is a *hypothesis* for the stationary result, not a demonstrated
dynamic effect or a novelty claim.

## Run

From this folder:
```sh
python -m pip install numpy scipy pytest
python -m pytest -q test_green_solver.py
python green_solver.py
```

## What remains to qualify as research

1. Compute constrained-energy stability and check the full wetting
   free energy, including contact-line dissipation and external work.
   Results for one rigid-channel meniscus must not be assumed valid
   for two menisci interacting with elastic walls.
2. Derive an energy-consistent dynamic law with the fluid mobility and
   reproduce the original constant-hysteresis limit before replacing it.
3. Verify physical basins and finite-time escape versus stationary
   arrest, grid convergence and sensitivity to pinning phase.
4. Contrast wavelength/correlation effects and the force/velocity
   protocol with the latest literature, particularly Chen et al. 2025.

This code does not revise Bradley et al.'s published conclusions.

## Reduced fixed-volume energy and local stability

An **additional model-level result** follows from eliminating the quasi-static
beam shape by its Green function and enforcing fixed liquid volume V.
With d=x_plus-x_minus and I the double integral of G on the wetted interval,

    E(x_minus,x_plus) = [V-d]^2/(2I)
                       - (nu/c_ref) integral_(x_minus)^(x_plus) c_Y(s) ds.

The first term is the Euler–Bernoulli bending energy; the second is the
position-dependent solid wetting energy in the same dimensionless convention
as the Laplace endpoint conditions. Direct differentiation gives
dE/dx_minus = p h_minus + nu*c_Y(x_minus)/c_ref and
dE/dx_plus = -p h_plus - nu*c_Y(x_plus)/c_ref, so stationarity is exactly
equivalent to the two local Young–Laplace conditions imposed by the solver.

At nu=4, V=.2, c0=.95, A=.025, phase=0, the computed *unconstrained interior
Hessian* eigenvalues classify the located candidates as:

- wavelength .25: 1 strict local minimum and 1 saddle;
- wavelength .15: 3 strict local minima and 3 saddles;
- wavelength .10: none found by the multistart search.

The model therefore exhibits several local resting configurations with the
*same* mean wettability and pattern amplitude but different wavelength.
This is a numerical local-minimum result, not proof that all roots have been
found or that the original hydrodynamic model reaches these states.

For illustration only, a diagonal positive mobility gives
dx/dt = -M grad(E) and dE/dt = -grad(E)^T M grad(E) <= 0.
The three minima at wavelength .15 attract small nearby perturbations in this
**synthetic gradient flow**. This is NOT a derived Cox–Voinov law or a
reproduction of transient capillary dynamics. Thus no physical escape rate,
basin probability, or hysteresis exponent is inferred.

Run energy checks with
`python -m pytest -q test_energy_stability.py`.
