# Volume-cycle metastability and sequential repinning (28 October 2026)

**Independent dimensionless model study. Not a reproduction of an experiment, a derivation of a microscopic defect energy, or a full moving-fluid simulation.**

## Mechanism and no-go control

This study uses the single conservative chemical energy and exact clamped/free
Euler–Bernoulli beam compliance in [beam_fold.py](beam_fold.py). Two contact
lines move on **the same prescribed Young-cosine landscape**, with no
independently imposed advancing or receding-angle threshold.

Write \(q=(a+b)/2\), \(d=b-a\), \(R=V-d<0\). Exact compliance is

\[
I(d,q)=d^5/120+d^2q^3/3.
\]

For a **homogeneous** Young cosine (\(A=0\)), at fixed \(d\)

\[
\partial_q E = -R^2 d^2q^2/(2I^2)<0.
\]

Thus this particular suction-driven, clamped/free beam has no *interior*
center equilibrium when the substrate is homogeneous. The channel's elastic
asymmetry biases the droplet towards the free end; chemical patterning can
create local minima that oppose that bias.

For a patterned surface the exact center derivative adds

\[
\partial_q E_{\rm chem}=-\nu[\kappa_Y(q+d/2)-\kappa_Y(q-d/2)].
\]

This is a force contrast across the two lines, not a made-up dry-friction
window. The interplay of this contrast and the **fully coupled** elastic
energy yields distinct metastable contact-line pairs.

## Protocol

- Parameters: \(\nu=4\), \(c_{\rm ref}=0.9\), \(\kappa_Y(x)=1+0.1\sin(2\pi x/0.075)\), with \(0.81\leq\cos\theta_Y(x)\leq0.99\).
- Begin in a specific stable interior well at \(V=0.22\): \((a,b)=(0.468814814433,0.706609934261)\).
- Reduce \(V\) quasi-statically to \(0.203\) and raise it to \(0.22\). Continue only *nearby stable roots*; when the branch disappears, integrate the local reduced gradient kinetics \(\dot z=-M\nabla E\) to the next well.
- Default illustrative mobility \(M=I_2\). This positive mobility ensures \(\dot E=-\nabla E^\top M\nabla E\leq0\) during a jump at fixed volume. **It is not the original lubrication/beam mobility.**
- Check each autonomous jump's energy decrease against time-integrated dissipation, positive wall gap, and outlet non-crossing. Avoid unconstrained global energy minimization, which would skip metastable wells.

## Independent saddle-node calculations

Solve \(E_a=E_b=\det E''=0\), and verify exactly one zero Hessian eigenvalue,
the other positive.

| Event | \(V_{\rm fold}\) | \(a_{\rm fold}\) | \(b_{\rm fold}\) |
| --- | ---: | ---: | ---: |
| Release 1 | 0.213851157109733 | 0.485098866979384 | 0.715416013259478 |
| Release 2 | 0.207647266128320 | 0.558542832798029 | 0.791063295145950 |
| Release 3, **outside this volume cycle** | 0.198698333854667 | 0.631382720272622 | 0.867426759371862 |

These are *local loss-of-stability volumes*, not fitted empirical CAH angles.
The first two are crossed during the protocol. The third would be crossed
only if the volume were reduced further; outlet escape must then be handled
as a distinct event, **not mislabeled as another repinning well**.

## Numerical result

At 136 steps per half-cycle, the two actual repinning steps occur at grid
volumes 0.21375 and 0.207625; the analytical fold volumes are independent of
that discretized control increment.

| Quantity | Initial | After complete down/up volume cycle |
| --- | ---: | ---: |
| Volume \(V\) | 0.22 | 0.22 |
| Left contact \(a\) | 0.4688148144 | 0.6056612548 |
| Right contact \(b\) | 0.7066099343 | 0.8707248986 |
| Center \(q\) | 0.5877123743 | 0.7381930767 |
| Reduced energy \(E\) | -0.9140478159 | -0.9605332614 |

The net center shift is **0.15048070, or 2.0064 chemical pattern periods**.
This is not an exactly quantized displacement (the width changes too);
it is a reproducible, parameter-specific *one-time relocation*.

For a complete path \(dE=p\,dV-D\,dt\), with pressure \(p=\partial E/\partial V\),
the volume-work integral and final minus initial energy give the total
irreversible energy loss \(D_{\rm cycle}=\int p\,dV-\Delta E\).

| Steps per half-cycle | Net \(D_{\rm cycle}\) |
| ---: | ---: |
| 68 | 0.043658760845 |
| 136 | 0.043656710568 |
| 272 | 0.043656385871 |
| 544 | 0.043655674653 |

At 136 steps: volume work \(-0.002828734922\), final minus initial energy
\(-0.046485445490\), and net irreversible loss \(0.043656710568\).
Individual repinning jumps lose approximately \(0.0179141\) and \(0.0256915\)
at that resolution. The small difference between their sum and cycle
dissipation is attributable to relaxation and finite volume increments along
the smooth portions of the discrete trajectory.

## Controls and limitations

1. **Locality:** Replacing the local continuation with global minimization
   can jump all the way to the outlet and yields a spurious apparent
   catastrophic escape. We follow local gradient flow at each fold instead.
2. **Mobility sensitivity:** The two tested jump endpoints agree within
   \(3\times10^{-7}\) when positive diagonal left/right mobilities are varied
   through \((0.2,1),(1,1),(1,0.2),(5,1),(1,5)\). This does **not**
   establish independence from the *true* coupled hydrodynamic mobility.
3. **Not an autonomous pump:** A second identical volume cycle from the
   relocated final state has no jumps or net motion in the tested interval.
   Calling this persistent/continuous ratcheting would be incorrect.
4. **Energy honesty:** The pressure-volume work and jump dissipation are
   thermodynamic properties **within the reduced energy/kinetics only**.
   Neither bulk viscous losses nor moving-meniscus boundary work of the
   full lubrication PDE has been explicitly closed.
5. **Finite geometry:** Contact lines must satisfy \(0<a<b<1\),
   \(0<V<b-a\) and noninterpenetration. A further reduction toward the
   third fold may produce an outlet event, not a sustainable interior well.
6. **Numerical precision:** Fold locations are solved with analytic
   gradients/Hessians; cycle energy loss is converging with volume step.
   A temporal convergence campaign of the full lubrication/beam solver and
   an experimental comparison are still missing.

## Reproduce

From the repository root:

    python -m pip install numpy scipy
    python research/microscopic_wetting/quasistatic_cycle.py --steps 136 --output cycle.csv
    python -m unittest discover -s research/microscopic_wetting -p 'test_*.py' -v

For a refined energy estimate:

    python research/microscopic_wetting/quasistatic_cycle.py --steps 544

## Literature priority and novelty check

Do **not** claim that pattern-induced hysteresis, metastable droplet states,
or directional microfluidic transport are new in isolation:

- Bradley, Hewitt & Vella (2021), *Phys. Rev. Fluids* **6**, 114003,
  doi:10.1103/PhysRevFluids.6.114003 (original bendotaxis setting).
- *Width effect on contact angle hysteresis in a patterned heterogeneous
  microchannel*, *Journal of Fluid Mechanics* (2022), includes evaporation /
  condensation and jumping contact lines:
  https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/width-effect-on-contact-angle-hysteresis-in-a-patterned-heterogeneous-microchannel/0032450EBB1389DD630AC72CFFDE6E6C
- Chen et al. (2025), *Capillarity* **15**, 4–11, doi:10.46690/capi.2025.04.02,
  already couple patterned wetting, hysteresis and deformable channels.
- Millet et al. (2026), doi:10.1002/nag.70363, derive a variational
  formulation including contact-line hysteresis energy.
- *Spontaneous droplet transport on shape-evolving microfiber rails*,
  *Nature Communications* (2025), https://www.nature.com/articles/s41467-025-65884-9,
  demonstrates ratchet-like directional contact-line transport.

The presently defensible contribution is a **reproducible, reduced
counterexample/control and quantitative two-line energy landscape** that can
form the basis of new experimentally testable predictions. Manuscript-level
novelty would need the coupled full-PDE and measured-pattern validation,
including whether the predicted sequence and energy accounting survive it.


## Minimum chemical contrast for a local pinning well

At **fixed** volume \(V=0.22\) and the same wavelength and beam properties,
there is a second distinct simple fold when **pattern amplitude** \(A\) is
used as the control parameter. Solving \(E_a=E_b=\det E''=0\) gives

\[
A_c=0.0470208189241591,\quad
(a_c,b_c)\simeq(0.48196748,\,0.71958230).
\]

This is a local contrast threshold for the *reference* well, not a global
guarantee that no other traps exist anywhere in the channel below \(A_c\).

Set \(\mu=A-A_c>0\), and let \(u\) be the normalized soft Hessian eigenvector.
Expanding about the fold gives a standard saddle-node normal form,

\[
E_{\rm soft}(\xi;\mu)
  =E_c+s\mu\xi+\tfrac{t}{6}\xi^3+\cdots,\quad
s=u\cdot\partial_A\nabla E,\quad
t=\partial^3_{\xi}E.
\]

Numerical independent checks give
\(s\approx-2.8511776210\),
\(t\approx472.57426\), and the *derived*, not fitted, barrier coefficient

\[
\Delta E \sim
\frac{4|s|}{3}\sqrt{\frac{2|s|}{|t|}}\,\mu^{3/2}
\approx 0.41759484\,\mu^{3/2}.
\]

Solve both nearby equilibria for seven contrasts between
\(\mu=10^{-4}\) and \(10^{-6}\); their computed barrier values have
log–log slope approximately **1.499819**, in agreement with the generic
saddle-node exponent \(3/2\). At \(\mu=10^{-4}\), the measured model barrier
is approximately \(4.17194\times10^{-7}\), versus the leading prediction
\(4.17595\times10^{-7}\).

The exponent is **not** a novel universal physical discovery. It checks
that the energy landscape has the required non-degenerate bifurcation
structure. It creates a computable contrast threshold and near-threshold
barrier for subsequent geometry/chemistry matching and experiments.

Reproduce:

    python research/microscopic_wetting/contrast_threshold.py

The code independently finds the minimum, saddle, curvature and
energy barriers, and tests the asymptotic coefficient/exponent.
