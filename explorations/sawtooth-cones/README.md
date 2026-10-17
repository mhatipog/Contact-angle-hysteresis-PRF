# Exploratory Young–Laplace contact-line mechanics on sawtooth cones

Independent, exploratory extension of the contact-angle-hysteresis theme. **Not part of Bradley et al. (2021), and not a reproduction of their bendotaxis model.** This self-contained subdirectory retains the original repository's code and attribution unchanged.

## Goal

Investigate whether constant-mean-curvature liquid barrels, Gibbs-style pinning and partial prewetting explain contact-line motion on toothed conical fibres. The new numerical check uses the **actual Young–Laplace BVP**, the true polygonal solid area, and independently differentiates its capillary energy with respect to volume and contact-line position.

At pinned endpoints \(z_-\) and \(z_+\), using the axisymmetric radial surface \(p(z)\),

\[
V=\pi\int_{z_-}^{z_+}[p(z)^2-r(z)^2]\,dz, \qquad
E=\gamma[A_{lv}-A_{\mathrm{reference}}].
\]

The reference area uses \(\cos\theta_Y\,dA_s\) for dry faces and \(dA_s\) for ideal prewetted faces. With fixed wetted pattern and away from corners, the pressure and endpoint first variations should obey

\[
\frac{\partial E}{\partial V}=\gamma\kappa,\quad
\frac{\partial E}{\partial z_+}=2\pi\gamma r\left[\frac{1+p'r'}{\sqrt{1+p'^2}}-c\sqrt{1+r'^2}\right],\quad
\frac{\partial}{\partial V}\frac{\partial E}{\partial z_+}
=\frac{\partial}{\partial z_+}(\gamma\kappa).
\]

These are envelope-theorem/Maxwell-type **consistency identities**, not newly discovered physical laws. The `thermodynamic_checks.py` script numerically tests them at both dry and prewetted ramp locations, independently of the BVP's differential residual.

## Reproduce

Requires Python 3.10+, NumPy, and SciPy. In this directory:

```bash
python -m pip install -r requirements.txt
python thermodynamic_checks.py
```

Units: coordinates in **mm**, volume in **mm³**, and energy in **pJ**. The console output reports dimensionless relative errors; assertions fail if the independent identities disagree.

## Experimental basis versus assumptions

Published SI for Omer, Raza and Zhang (2025) reports teeth pitches of 10/20/40 μm, facet angle \(79\pm5^\circ\), prewetted gaps and apparent static/advancing/receding angles on planar resin of about 35°/52°/0°. **None of these is an independent measurement of the local equilibrium Young angle for all toothed faces.** The model's 20 μm pitch, 3° half-angle, 95% gentle facet, 90% ideal-film coverage, 3 nL, and 35° dry-angle proxy are **illustrative assumptions**. The experimental snap trajectory is not fitted or predicted here; no video measurement has been made.

The BVP holds its rear boundary fixed, uses an axisymmetric zero-thickness wet-film reference and checks nonpenetration on a discrete evaluation mesh. It does **not** solve dynamic Navier–Stokes, a moving film equation, or two freely evolving contact lines. Contact-line jumps cannot be inferred merely from equal-energy states or endpoint-angle agreement.

A useful analytical continuation is the distinction between a **finite-slope loss of a remaining dry facet** (barrier \(\Delta E\sim C\ell\) when the endpoint force has a nonzero finite limit) and a **smooth saddle-node** (barrier \(\Delta E\sim\varepsilon^{3/2}\) under nondegenerate normal-form conditions). Neither exponent is claimed novel. A subsequent study should resolve rounding of the edges, film advance, and admissibility of both contact lines, then preregister geometric parameters before comparison with front/rear high-speed trajectories.

## References and credit

- Bradley et al. (2021), *Droplet trapping in bendotaxis caused by contact angle hysteresis*, *Physical Review Fluids* 6, 114003, https://doi.org/10.1103/PhysRevFluids.6.114003 — original code elsewhere in this repository; **independent topic and authorship**.
- Fuentes, Hatipogullari et al. (2019), *Contact line stick-slip motion and meniscus evolution on micrometer-size wavy fibres*, *J. Colloid Interface Sci.* 540, 544–553, https://doi.org/10.1016/j.jcis.2019.01.045.
- Omer, Raza and Zhang (2025), *Self-Propelling Water Droplets on Conical Spikes with Sawtooth Surface Structure*, *ACS Appl. Mater. Interfaces* 17, 31604–31613, https://doi.org/10.1021/acsami.5c03846.

The numerical code here is an independent exploratory implementation. No supplemental publisher PDFs or video files are redistributed.