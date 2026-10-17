# Spectral commensurability and compliant pinning — research only

Analytically soluble and numerically benchmarked **reduced** model extending
the exact two-meniscus wettability contrast with a *surrogate* longitudinal
elastic stiffness. No original MATLAB physics files are changed.

For a sinusoidal chemical Young cosine cY(x)=c0+A sin(kx+phi) and fixed
meniscus separation d, the chemical force on droplet center q equals

    f(q)=2 gamma [cY(q+d/2)-cY(q-d/2)]
        =4 gamma A sin(kd/2) cos(kq+phi).

Thus Fc=4 gamma |A sin(pi d/wavelength)|, with exact rigid cancellation
at integer-wavelength lengths. This is a size-selective Fourier filter.

Allowing the length to relax in

    E(q,d)=K(d-d0)^2/2-(4 gamma A/k)sin(kd/2)sin(kq+phi)

restores pinning at d0=m*wavelength at order Fc~2 gamma^2 A^2 k/K.
The first residual mode has twice the substrate wavenumber. This is a
surrogate model: K is NOT yet derived from the elastic cantilever and
lubrication flow in Bradley, Hewitt & Vella (2021).

For a tuned two-harmonic chemical pattern giving an opposing force
B cos(kq)-(B/4)cos(2kq), the pinning maximum is quartic and the
force-controlled overdamped velocity scales as excess_force^(3/4) rather
than the ordinary 1/2. This is engineered degeneracy, NOT a generic
new critical exponent. Small detuning restores the 1/2 exponent
arbitrarily close to threshold.

For a general stationary Young-cosine field the 2-point relation is
Var(f)=8 gamma^2[Var(cY)-Cov(cY(x),cY(x+d))].
This is about force fluctuation magnitude, not automatically escape.

These phenomena require validation with the fixed-volume Green-function
cantilever, a physically justified two-line mobility, and matched
measured advancing/receding thresholds before claiming a new
prediction about bendotaxis. In particular, patterned wetting and
flexible channels have been considered previously by Chen et al.
(Capillarity 2025, DOI 10.46690/capi.2025.04.02).

Run:

    pip install numpy scipy pytest
    python -m pytest -q research/commensurability

Related context: Hatipogullari et al. (2019), Phys. Rev. Fluids 4,
044008; Bradley et al. (2021), Phys. Rev. Fluids 6, 114003.

No unpublished manuscripts are included here.
