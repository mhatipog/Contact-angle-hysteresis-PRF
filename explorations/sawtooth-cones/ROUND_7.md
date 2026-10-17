# Round 7 — symmetry control and two-contact-line Gibbs screening (17 October 2026)

**Status: exploratory numerical research, not an experimentally validated explanation.**
Original MATLAB code by Bradley et al. remains untouched. This extension
concerns Omer et al.'s 2025 conical sawtooth specimens, not bendotaxis.

## 1. Periodicity gives a no-go control

For a cylinder with periodic radius r(z+lambda)=r(z) and a periodic
wetting-energy reference c(z+lambda)=c(z), shift BOTH contact lines and
the entire admissible Young-Laplace barrel by exactly lambda. The
axisymmetric volume and capillary functional are unchanged:

    V = pi integral [p(z)^2 - r(z)^2] dz
    E = gamma [2pi integral p sqrt(1+p_z^2) dz
               - 2pi integral c(z) r(z) sqrt(1+r_z^2) dz]

Consequently E(x_minus+lambda,x_plus+lambda,V)=
E(x_minus,x_plus,V) **for corresponding equilibrium branches**.
Tooth chirality can create asymmetric barriers but not a net conservative
energy slope across complete periods on a cylinder. Cone taper, volume
growth, film evolution, or another drive breaks these hypotheses.
This is a standard symmetry argument, **not a new physical theorem**.

The reproducible run `python periodic_no_go.py` tests each tooth
orientation at two periodic wetting fractions (0 and 0.9), using
21 translations per pitch and the resolved, fixed-volume Young-Laplace
surface. Example parameters: r0=80 um, lambda=20 um, tooth height=3 um,
duty=0.8, volume=3 nL, footprint=200 um, and a 35-degree *proxy*
for local dry reference energy.

For each cylindrical orientation: the intra-pitch variation is about
139.272 pJ, while the one-pitch energy drift is numerically zero to
the displayed precision. The periodic wetting fraction does not change
this particular full-footprint translation trace, since the footprint
covers an integer number of complete periods; it **does** affect
individual moving-endpoint pinning criteria.

For comparison, a smooth cone of slope dr/dz=0.05 and r0=80 um exhibits
an endpoint-to-endpoint drop of about 12.832 pJ for one 20 um shift; the
sawtooth cone gives about 18.577 pJ for the same illustrative endpoint
states. An energy difference between endpoints is NOT proof that every
intervening Young-Laplace surface is geometrically accessible.

## 2. Both contact-line coordinates need stability checks

The two endpoint first variations, evaluated from their Young-Laplace
surface slopes and the facet-local wetting references, must obey the
one-sided necessary condition at a corner:

    dE/dx_i from left <= 0 <= dE/dx_i from right, i in {rear,front}.

In an illustrative nonpenetrating cylindrical specimen (rear=816 um at
a crest, front=1020 um at a valley, V=3 nL), the dry 52-degree proxy
passes the **first-order** tests for both endpoints. A dry 35-degree
proxy fails at the front. Allowing the ideal wet region to extend onto
the following gentle facet also fails the front-corner test even with
the 52-degree proxy. This does NOT establish loss of all metastability:
the contact line could arrest at the new film edge, and second
variations, obstacle contacts, or kinetics can change the result.

The measured 35 +/- 3-degree *static apparent* angle and 52 +/- 5-degree
advancing angle reported for flat printed resin are **not equivalent
independent measurements of the equilibrium Young angle** on every
microscopic tooth facet. Hence both are deliberately labelled proxies.

## 3. What is already known and what remains genuinely open

This round is a **falsification/control**, not a journal manuscript.
In particular, the previous ideal-patch linear barrier law is
conditional on the endpoint force retaining a nonzero finite limit as
the dry interval shrinks. It is NOT a universal law for sharp corners:
Graña Otero and Parra Fabián (2019), *Contact line depinning from
sharp edges*, already analyze corner stability and saddle-node
depinning [arXiv:1907.03499](https://arxiv.org/abs/1907.03499).

A consequential paper would require: experimentally digitized tooth
facets and film coverage; source-video contact-line trajectories,
including rear displacement and jump count; a nonpenetrating
two-line (obstacle) solver with a moving wetting front; convergence
and uncertainty intervals; and a preregistered comparison across
S-10/S-20/S-40 and reversed imbricated teeth.

The near-term collaboration target is the experimental group of
Omer, Raza and Zhang, not journal submission. Ask for Movie S4/S5
coordinate time series, per-event geometry and wetting protocol before
claiming a reproduction.

## Run

From `explorations/sawtooth-cones`:

```bash
python -m pip install -r requirements.txt
python periodic_no_go.py
```

Two CSV files are created locally in `results/` and each run prints
asserted symmetry closure and two-endpoint one-sided sign checks.

## Literature

- Fuentes, Hatipogullari et al. (2019), J. Colloid Interface Sci. 540, 544–553, doi:10.1016/j.jcis.2019.01.045.
- Bradley et al. (2021), Phys. Rev. Fluids 6, 114003, doi:10.1103/PhysRevFluids.6.114003 — original independent bendotaxis code in this repository.
- Graña Otero and Parra Fabián (2019), *Contact line depinning from sharp edges*, arXiv:1907.03499.
- Omer, Raza and Zhang (2025), ACS Appl. Mater. Interfaces 17, 31604–31613, doi:10.1021/acsami.5c03846.
