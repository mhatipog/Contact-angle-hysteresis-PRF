"""Reduced fixed-volume Euler-Bernoulli beam: wavelength commensurability.

Minimize droplet length at FIXED center q. This is a local energetic
validation, not full transport, true trapping or hydrodynamic dissipation.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar

@dataclass(frozen=True)
class Beam:
    volume: float=.2
    nu: float=4.
    cosine_reference: float=.95
    cosine_mean: float=.95
    amplitude: float=0.
    wavelength: float=.2
    phase: float=0.

def green_integral(xm,xp):
    """Integral of clamped-free Green kernel on wetted rectangle."""
    if not(0<xm<xp<1):raise ValueError("Invalid contact-line locations")
    return ((xp**5-xm**5)/20-xm**3*(xp*xp-xm*xm)/6
            +xm**4*(xp-xm)/12)

def beam_energy(q,d,b=Beam()):
    a=q-d/2;z=q+d/2
    I=green_integral(a,z)
    k=2*np.pi/b.wavelength
    wet=b.cosine_mean*d+2*b.amplitude/k*np.sin(k*q+b.phase)*np.sin(k*d/2)
    return (b.volume-d)**2/(2*I)-b.nu*wet/b.cosine_reference

def stationary_length(q,b=Beam(),lo=None,hi=None):
    lo=b.volume+1e-5 if lo is None else lo
    hi=2*min(q,1-q)-1e-4 if hi is None else hi
    if hi<=lo:raise ValueError("Center too close to end")
    result=minimize_scalar(lambda d:beam_energy(q,d,b),bounds=(lo,hi),
                           method="bounded",options={"xatol":1e-14})
    if not result.success or min(result.x-lo,hi-result.x)<2e-5:
        raise RuntimeError("No verified interior length minimum")
    return float(result.x),float(result.fun)

def effective_stiffness(q,d,b=Beam(),step=7e-5):
    return (beam_energy(q,d+step,b)-2*beam_energy(q,d,b)
           +beam_energy(q,d-step,b))/step**2

def commensurate_experiment(q=.6, amplitudes=(.004,.008,.012,.016)):
    from dataclasses import replace
    base=Beam()
    d0,E0=stationary_length(q,base,lo=.20001,hi=.45)
    K=effective_stiffness(q,d0,base)
    phi=np.pi/2-2*np.pi*q/d0
    records=[]
    for A in amplitudes:
        perturbed=replace(base,amplitude=float(A),wavelength=d0,phase=phi)
        d,E=stationary_length(q,perturbed,lo=.20001,hi=.45)
        predicted=-(base.nu*A/base.cosine_reference)**2/(2*K)
        pressure=(base.volume-d)/green_integral(q-d/2,q+d/2)
        records.append((A,d0,d,K,E-E0,predicted,pressure))
    return records
