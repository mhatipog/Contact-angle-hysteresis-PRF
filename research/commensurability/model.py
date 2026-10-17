"""Two-line wettability filtering and elastic commensurability, reduced physics.

Not a validated hydrodynamic replacement for the Bradley-Hewitt-Vella PDE.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

@dataclass(frozen=True)
class Surface:
    wavelength: float = 0.2
    reference_length: float = 0.4
    amplitude: float = 0.025
    gamma: float = 1.0
    stiffness: float = 50.0
    phase: float = 0.0

    @property
    def k(self):
        return 2*np.pi/self.wavelength

def line_force(q, length, s=Surface()):
    """2gamma*(cosYoung(right)-cosYoung(left)), sinusoidal heterogeneity."""
    return 4*s.gamma*s.amplitude*np.sin(s.k*length/2)*np.cos(s.k*q+s.phase)

def rigid_threshold(length, s=Surface()):
    return 4*s.gamma*abs(s.amplitude*np.sin(np.pi*length/s.wavelength))

def rigid_velocity(excess_force, length, s=Surface()):
    fc=rigid_threshold(length,s)
    if excess_force <= fc:
        return 0.0
    return float(np.sqrt(excess_force**2-fc**2))

def relaxed_length(q, s=Surface()):
    """Strictly convex in length when stiffness > gamma*abs(A)*k."""
    if s.stiffness <= s.gamma*abs(s.amplitude)*s.k:
        raise ValueError("Unique relaxed length not guaranteed")
    q=np.asarray(q,dtype=float)
    angle=s.k*q+s.phase
    d=np.full_like(q,s.reference_length)
    for _ in range(25):
        z=s.k*d/2
        residual=s.stiffness*(d-s.reference_length)-2*s.gamma*s.amplitude*np.cos(z)*np.sin(angle)
        deriv=s.stiffness+s.gamma*s.amplitude*s.k*np.sin(z)*np.sin(angle)
        step=residual/deriv
        d=d-step
        if np.max(np.abs(step))<3e-14:
            break
    if np.max(np.abs(s.stiffness*(d-s.reference_length)-
                     2*s.gamma*s.amplitude*np.cos(s.k*d/2)*np.sin(angle)))>3e-11:
        raise RuntimeError("Relaxed length did not converge")
    return d

def relaxed_energy(q, s=Surface()):
    d=relaxed_length(q,s)
    return .5*s.stiffness*(d-s.reference_length)**2 - (
        4*s.gamma*s.amplitude/s.k*np.sin(s.k*d/2)*np.sin(s.k*np.asarray(q)+s.phase))

def relaxed_force(q, s=Surface()):
    return line_force(q,relaxed_length(q,s),s)

def relaxed_threshold(s=Surface(), points=4096):
    q=np.arange(points)*s.wavelength/points
    v=relaxed_force(q,s)
    idx=int(np.argmin(v))
    step=s.wavelength/points
    opt=minimize_scalar(lambda y:float(relaxed_force(y,s)),bounds=((idx-1)*step,(idx+1)*step),
                        method='bounded',options={'xatol':1e-14})
    return -float(opt.fun)

def asymptotic_threshold(s=Surface()):
    if abs(np.sin(np.pi*s.reference_length/s.wavelength))>1e-10:
        raise ValueError("Need integer-wavelength reference length")
    return 2*s.gamma**2*s.amplitude**2*s.k/s.stiffness

def spectral_force(q,length,wavelength,amplitudes,phases=None,gamma=1.):
    phases=[0.]*len(amplitudes) if phases is None else phases
    if len(amplitudes)!=len(phases):raise ValueError("Unmatched modes")
    return sum(4*gamma*a*np.sin((n+1)*np.pi*length/wavelength)*
               np.cos((n+1)*2*np.pi*q/wavelength+phase)
               for n,(a,phase) in enumerate(zip(amplitudes,phases)))

def force_rms_from_covariance(samples,shift,gamma=1.):
    """Exact for equally spaced periodic samples; variance != max pinning force."""
    x=np.asarray(samples,float)
    anomaly=x-np.mean(x)
    var=np.mean(anomaly**2)
    cov=np.mean(anomaly*np.roll(anomaly,-shift))
    predicted=2*gamma*np.sqrt(max(0.,2*(var-cov)))
    direct=2*gamma*np.sqrt(np.mean((np.roll(x,-shift)-x)**2))
    return float(predicted),float(direct)

def two_harmonic_design(wavelength=.2,length_ratio=.25,A1=.025):
    a=np.pi*length_ratio
    A2=-A1/(8*np.cos(a))
    b1=4*A1*np.sin(a)
    return (A1,A2),b1

def quartic_velocity(excess, B):
    """Exact period average for tuned B*cos(y)-(B/4)*cos(2y)."""
    if excess<=0 or B<=0:raise ValueError("Positive excess and B required")
    integral=2*quad(lambda y:1/(excess+.5*B*(1-np.cos(y))**2),
                    0,np.pi,epsrel=5e-10,limit=650)[0]
    return float(2*np.pi/integral)

def generic_velocity(excess,B):
    return float(np.sqrt(excess*(2*B+excess)))
