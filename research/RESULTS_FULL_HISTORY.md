# Phase-sensitive wetting, full contact-line histories and energetic yield

**Model-level independent research note — 17 October 2026**

## Question

Can two chemically patterned surfaces with the same Young-cosine
distribution and advancing/receding interval cause different *first*
contact-line release outcomes in an elastic fluid-filled channel, even
when both contact lines experience their particular surface from the
initial state?

## Validated answer (conditional on constitutive assumptions)

Yes. We solve the original Bradley-Hewitt-Vella fluid–beam discretization
with a spatially varying advancing pressure multiplier on BOTH moving
contact lines and a local receding yield interval. We use its original
unmodified RHS, only replacing how boundary pressures are selected.
Spatial pattern: kappa_A(x)=1+A sin(2pi*(x-.45)/.25), where A=±.002.
Local receding multiplier kappa_R=kappa_A+.0293. Choosing a reference
advancing-angle cosine .94 keeps the implied local advancing and
receding cosines in the physically admissible interval (0,1).

At nu=4, V=.2, initial x_plus=.65 and x_minus=.45, the N=20 solver
at mu=0 reports:

| Signed pattern amplitude | First left pin t | First receding event t <= .1 |
|---|---:|---:|
| +.002 | .00349275 | none by t=.1 |
| 0 | .00330396 | .0141715 |
| -.002 | .00329147 | .0101992 |

At mu=.01: A=+.002 remains pinned through t=.1; A=-.002
releases near .0108028. The same classification passes N=12/16/20.
Individual first pinning *times* are still sensitive to the initial
event-solver step and are not claimed grid-converged; pinned positions,
late equilibrium and qualitative first-release classification are much
more consistent. Unreleased through t=.1 is NOT proof of permanent
trapping. Subsequent release-to-escape dynamics is not included.

## Exact energy bridge to the user's equilibrium theory

The reduced clamped Euler–Bernoulli beam compliance is exactly

    I(a,b)=(b-a)^2*((b-a)^3+40*((a+b)/2)^3)/120.

With conserved droplet volume V and local Young multiplier kappa_Y,

    E(a,b)=(V-(b-a))^2/[2I(a,b)]-nu*integral_a^b kappa_Y(x) dx.

For the separately solved stationary beam with right meniscus free,
the generalized force derivatives obey EXACTLY

    dE/db = 0,
    dE/da = -nu * [Gamma - A sin(2pi*(a-xref)/ell)]

where Gamma=-p*h(a)/nu-1 is the same pressure-based yield observable
in the time-dependent solver. No fitting: 25 independently computed
quartic beam equilibria agree to <4e-13.

This is the correct route to energetic dry-friction/contact-line
complementarity (Moreau-type rate-independent inclusions). In a
separate *reduced* gradient-flow benchmark,

   0 in grad(E)+Y*partial|velocity|+zeta*velocity

one obtains identically

   dE/dt = -Y*sum|velocity|-zeta*sum(velocity^2) <= 0.

The code checks 1,000 additional random admissible configurations.
**This is not a demonstrated full energy balance for the moving
lubrication PDE**, which also requires bulk viscous dissipation,
meniscus surface energy and mechanical boundary work.

## Contact-line kinetics

The full-history benchmark employs the linearized local advancing
law kappa_app=kappa_Y-mu*s_i*velocity_i. With s_left=-1 and
s_right=+1, this has nonnegative local contact-line dissipation.
It is solved through an exact 2x2 Schur complement of the ORIGINAL
fluid–beam discretization, avoiding a new nonlinear PDE system.

An attempt to apply the full angle-cubed Cox/Voinov law throughout
the stiff startup *failed*: trial states exceeded its positive-angle
validity domain. We did not clip negative angle cubes. The exact
angle-cubed constitutive law is separately checked pointwise on
slow, near-equilibrium original-PDE states, where it reproduces the
linearized tangent law and gives nonnegative contact-line dissipation.
This does NOT verify an integrated all-time Cox evolution.

## Scientific significance and next decisive test

Compared with previous isolated phase-pattern or prepared-pinning
demonstrations, the release/no-release distinction survives separate
surface-dependent starting histories. Its energetic driving quantity
is independently related to the variation of the actual reduced
clamped-beam energy.

However, a *pattern-dependent threshold* by itself is not established
as novel: heterogeneous wettability and flexible-channel droplet
transport have substantial existing literature (e.g., Chen et al.,
Capillarity 2025, DOI 10.46690/capi.2025.04.02; 2026 phase-field
studies of multistable wetting). Novelty would require a distinct,
quantitatively validated prediction, ideally compared with
measured/local microscopic hysteresis and fully converged trajectories.

The next decisive physics task is to derive the advancing and receding
pinning laws from a single specified microscopic interfacial energy
landscape (instead of imposing an additive phenomenological interval),
then verify full moving-fluid energy closure and solve the subsequent
receding and possible repeated-pinning stages.

## Reproduce

- `octave --quiet --no-gui --eval "run('research/test_fullhistory_linear.m')"`
- `python research/energy_yield_identity.py`
- `octave --quiet --no-gui --eval "run('research/test_fullhistory_cox.m')"` (near-equilibrium *pointwise* only)

This is work on a separate fork research branch, not an upstream
scientific correction claim. The author's original source remains
unmodified, and the pre-existing upstream Figure 11 PR remains open.
