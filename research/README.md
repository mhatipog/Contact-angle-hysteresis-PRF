# Independent pinned-meniscus equilibrium check

Scientific reference: Bradley, Hewitt & Vella, *Physical Review Fluids* **6**, 114003 (2021), DOI https://doi.org/10.1103/PhysRevFluids.6.114003, Eqs. (33)–(36), Figure 5.

This folder independently reconstructs the **stationary** wetting-droplet configuration: pinned left contact line, advancing right contact line. It imposes clamped-left beam continuity, free-right beam conditions, right-meniscus Laplace pressure and volume conservation. It does not invoke the original time-dependent solver or modify any of the published source files.

Run from the repository root with GNU Octave or MATLAB:

```sh
octave --quiet --no-gui --eval "run('research/test_steady_pinned_equilibrium.m')"
```

For `nu=4`, `V=0.2`, `x_minus=0.45`, the stationary calculation gives `x_plus≈0.660902319`, `p≈-4.28495113` and `lambda_e≈0.03149606`. The paper reports a trapped-state asymmetry around 0.03. This is a *steady-state*, fixed-contact-line benchmark, **not** a reproduced trajectory or proof of dynamical stability.

## Why the separately identified contact-line bugs may not affect the paper

Section II.D (page 7) of the paper says that its wetting-droplet initial conditions always keep the **right** meniscus advancing. The repository includes code paths for right pinned/receding cases but those are outside the main wetting trajectories; defects in those paths should not be advertised as affecting Figure 5 or Figure 11 without a counterexample.

In a separate preliminary Python port, coarse `N=8/12` time evolution falsely leaves the pinned-left state shortly after reaching it, whereas `N=16` remains pinned through the tested interval `t=0.02`. At the first pinning event, the estimated left cosine was `0.95780185, 0.95998635, 0.96069454` for `N=8,12,16`, versus the advancing cosine `1/1.04=0.96153846`. The first two exceed the MATLAB event tolerance (`1e-3`) on the wrong side due to finite-difference error. This is an **unvalidated-port resolution diagnostic**, not proof of any fault in the authors' default `N=25` simulations.

Keep this investigation separate from the Figure 11 regeneration PR; no scientific correction is proposed yet.

## Archived Figure 11 data provenance

The paper's Section V describes the initial bisection bounds for V=0.2 as
[V+0.03, 0.97], while the deposited Generate_numerics_data_Fig11.m script
uses [V+0.01, 1]. Running the independent script
python -m research.fig11_provenance
against Data/numerics_comparison_data.mat establishes that all 41 finite
archived thresholds fit the eight-step bisection lattice of the deposited
code's interval, and none fit the article's stated interval at machine
precision. This is a provenance difference, NOT a demonstrated problem
with the trapping model or figures. Original generation history is unknown.

## Pilot physics continuation: independent threshold discrepancy

The archived `Data/numerics_comparison_data.mat` stores **initial positions**
above which a full dynamic integration escaped during the authors' finite
simulation horizon. They are not equilibrium meniscus positions.

Run `python -m research.steady_threshold_comparison` to independently solve
Eqs. (33)–(42) and compare its stationary `X_plus` against the deposited
dynamic escape thresholds. Across 22 valid (lambda, nu) pairs, the discrepancy
changes sign with lambda. For nu=4, it is +0.04156 for lambda=0.02, +0.00093
for 0.04, and -0.01993 for 0.10. This does not demonstrate an error: the paper
states that the equilibrium-versus-initial-position comparison is approximate.
However, it motivates checking initial squeezing, subsequent depinning and
finite stopping time.

`research/test_original_first_pinning.m` is a narrow integration of the
original (unmodified) `fout.m` with `ode15s` and `EventsFcn.m`, **without**
PMAD and without changed event equations. It compares the first pinning
time and contact-line positions to an independent SciPy finite-difference
reconstruction. It is not a full trapping/escape test. All research scripts
remain on this independent branch and are not part of the upstream PR.

## Original-model hysteresis transition check

The Octave test research/test_original_hysteresis_divergence.m runs the authors'
UNMODIFIED fout.m, EventsFcn.m and update_pinning_conditions.m up to time
0.03 for N=16, nu=4, V=0.2, xplus0=0.65. It checks that lambda_max=0.02
reaches the left receding state after first pinning, while lambda_max=0.04
remains pinned. This is a *short-time* event check, not a proof of full-time
escape versus indefinitely stable trapping or numerical convergence.
