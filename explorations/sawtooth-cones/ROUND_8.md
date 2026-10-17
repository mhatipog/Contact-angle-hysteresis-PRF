# Round 8 — SI-consistent geometry and obstacle contact (17 October 2026)

**Exploratory research / tested numerical diagnostics — NOT a validated experimental reproduction or a manuscript-ready prediction.**
This is an independent extension within a public archive of Bradley et al. (2021)'s bendotaxis code. It does not alter or represent the original bendotaxis simulations.

## The actual experimental constraint

Figure S4(b) of Omer, Raza and Zhang (2025), Supporting Information to doi:10.1021/acsami.5c03846, defines **beta = 79 ± 5 degrees from the vertical** to the *long, descending* tooth face. This must NOT be interpreted as a 79-degree inclination from the horizontal. The nearly vertical rising reset follows that facet when moving toward the broader end of a cone.

Let the cone half-angle be alpha, tooth pitch lambda, and beta measured from vertical. In the ideal two-facet approximation, the descending facet has
```
dr/dz = -cot(beta).
```
If the steep reset occupies an *unmeasured* fraction epsilon of one pitch, set d=1-epsilon and use
```
r(z) = r0 + tan(alpha)*z + a * tooth_profile(z/lambda)
a = d*lambda*[tan(alpha)+cot(beta)],
tooth_profile = 1 - f/d                 (0 <= f < d)
                (f-d)/(1-d)            (d <= f < 1)
```
so that the long facet has the measured slope. The shape is a polygonal **hypothesis**, not reconstructed SEM metrology. epsilon=0.01 or 0.03 is a numerical regularization, NOT measured 0.1-1 micrometre print resolution.

At pitch 20 µm, alpha=3 degrees, beta=79 degrees, the model predicts:
- epsilon=.01: reset amplitude 4.886 µm, planar path-length ratio 1.252.
- epsilon=.03: reset amplitude 4.788 µm, planar path-length ratio 1.229.

For a vertical ideal reset, SI Eq. S1 predicts ~1.265 with alpha=3 degrees, beta=79 degrees; the paper reports ~1.25. These are **not independent measurements** of amplitude and reset width: Eq. S1 uses the observed beta and geometric construction. A full axisymmetric *wetted area* ratio need not equal this planar path-length ratio, because the tooth radius varies over each facet. No independent tooth height or reset width has been measured here.

The earlier arbitrary 20 µm pitch / 3 µm amplitude / 80%-20% facet split gave a planar path length ratio ~1.064 and therefore should NOT be used as experimental metrology.

## Elegant geometric feasibility condition

The axisymmetric Young-Laplace barrel has radial free surface p(z) constrained by a tooth surface r(z), with p(z_i)=r(z_i) and **nonpenetration** p(z)>=r(z). This immediately implies necessary one-sided endpoint inequalities:
```
p'(z_rear+) >= r'(z_rear+),
p'(z_front-) <= r'(z_front-).
```
At a rear valley the solid rises along an almost vertical reset. In the sharp-reset limit its right derivative diverges as epsilon -> 0, and a regular finite-slope free surface CANNOT depart from that valley without intersecting the solid. Its actual apparent meniscus must instead contact the next facet/tooth, be supported by a wetting film, or depart with near-vertical tangent. This is a local feasibility condition, not a new contact-angle law.

This local test is useful because it can rule out an inadmissible Young-Laplace solution before globally resolving all intermediate tooth contacts. Passing the local inequality does **not** prove global nonpenetration.

## Parameter-screen definition and findings

Run `python si_round8.py`. The scan uses:
- pitches 10, 20, 40 µm;
- beta 74, 79, 84 degrees (reported mean ± 5, *scenario bounds*, not a confidence interval);
- unmeasured reset fractions 1% and 3%;
- illustrative droplet volumes 2.5, 3.0, 3.5 nL;
- rear/front corners in specified windows near 0.8 and 1.0 mm (axial coordinates illustrative);
- alpha = 3 degrees, r0 = 30 µm (both illustrative rather than specimen-specific SEM fits);
- two **proxy** local dry Young angles (35 and 52 degrees); patterned wet fractions 0, .5, .9, whose actual experimental values are unknown.

**Numerical outcomes, for this defined screen only:**

| Pitch | Candidate BVPs | Geometrically admissible |
|---|---:|---:|
| S-10 | 648 | 324 |
| S-20 | 648 | 325 |
| S-40 | 450 | 266 |
| TOTAL | 1,746 | 915 |

1,731 BVPs converged. Of the 816 converged solutions penetrating the solid, 811 also failed the **local** endpoint-slope test; all 915 globally nonpenetrating solutions passed it. Five nonpenetrating-screen failures passed the local gate, which proves it is necessary, not sufficient. The implementation checks **all polygon vertices** as well as a dense profile grid. A converged BVP is *never* by itself labeled physical.

The 915 geometrically admissible endpoint pairs enter 5,490 proxy wetting/angle combinations for one-sided **first-order Gibbs** corner signs. None of the tested pairs passes the rear AND front necessary corner-minimum signs simultaneously. This is **NOT a demonstration that the real experiment has no metastable droplets**. The screen samples corner-bound lines only, does not free the film edge or interior obstacle contacts, and lacks full second-variation stability.

A smaller independent taper check at cone half-angles 2, 3, 4 degrees for beta=79 degrees, epsilon=.01 and volume=3 nL also found zero simultaneous Gibbs passes across the checked corner pairs. This is an illustrative local robustness check, not a full uncertainty propagation.

## A constructive replacement: the obstacle problem

Run `python obstacle_case.py`. Consider beta=79 degrees, S-20, epsilon=.01, z_rear=(40+d)*20µm (a rear valley), z_front=51*20µm (crest), and V=3 nL.

The unconstrained *converged* CMC boundary-value solution intersects the sawtooth by **4.712111 µm**. We then solve the discrete constrained variational problem

```
min A[p] = 2*pi*integral p(z)*sqrt(1+p_z^2) dz,
subject to p(z) >= r(z), p(z_rear/front)=r(z_rear/front),
pi*integral[p(z)^2-r(z)^2] dz = 3 nL,
```
on a mesh containing **every polygon vertex**. The area functional is discretized as exact conical frusta, the volume as exact frusta differences, and the gradient supplied analytically to SLSQP.

| Maximum mesh spacing | Nodes | Area (mm2) |
|---|---:|---:|
| 4 µm | 62 | 0.144057677140 |
| 2 µm | 112 | 0.144054992270 |
| 1 µm | 212 | 0.144054321376 |

Each computed solution has at least one additional interior contact, nonnegative gap to optimization tolerance, and negligible excess-volume error. The last-two-mesh relative area change is ~4.7e-6.

**Limitations:** SLSQP convergence and mesh refinement do not prove a global continuum minimum. Moreover, this is a *fully wetted/reference-interface geometric limit*: dry solid–liquid interfacial-energy accounting, moving film edge, line friction and 3-D hydrodynamics are **not** solved. A contact with the obstacle changes the physical partition of interfaces; do not interpret the resulting profile as a two-free-dry-line metastable branch.

## Why this is worth pursuing — but not yet publishable

The robust outcome is that **a metrology-aware polygonal model makes a single spanning Young-Laplace barrel insufficient for some nominally pinned pairs**. Local geometric tests explain almost all intersections, and an explicit obstacle solve repairs an example. This provides a sound numerical pathway toward a contact-set and film-front model.

This is not a universal mathematical novelty claim; Gibbs pinning and sharp-edge saddle-node instabilities are established, including Graña-Otero and Parra Fabián, *Physical Review Fluids* 4, 114001 (2019). A strong paper would require:
1. independent **actual** facet/film geometry and reset widths, not just beta and Eq. S1;
2. Movie S4/S5 trajectories: z_front(t), z_rear(t), film-edge positions, tooth-jump count, radii, time calibration;
3. mixed wet/dry energy with moving film and contact-set complementarity;
4. matched forward/reversed-tooth and S-10/S-20/S-40 predictions **without retuning**;
5. continuum convergence and measurement uncertainty intervals.

**Suggested collaboration target:** Omer / Raza / Zhang experimental group to request their calibrated trajectories and geometry. No outreach has been sent.

References:
- Omer, Raza and Zhang (2025), *ACS Applied Materials & Interfaces*, doi:10.1021/acsami.5c03846, Supporting Information Fig. S4 and Eq. S1.
- Fuentes, Hatipogullari et al. (2019), *Journal of Colloid and Interface Science* 540, 544–553, doi:10.1016/j.jcis.2019.01.045.
- Bradley et al. (2021), *Physical Review Fluids* 6, 114003, doi:10.1103/PhysRevFluids.6.114003 (original independent code in this repository).
- Graña-Otero and Parra Fabián (2019), *Physical Review Fluids* 4, 114001, doi:10.1103/PhysRevFluids.4.114001.
