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
