"""Zero-film-thickness heterogeneous wet-coverage energy model.
A_wet coverage is spatial: fully prewetted steep faces, followed by a
specified initial proportion of each gentle tooth facet; outside remains dry.
This is a limiting MODEL, not measured wet coverage or a dynamic film PDE.
"""
import numpy as np
from math import floor,ceil,cos,pi

def energy_weighted(geo,z1,z2,thetaY,wet_ramp_fraction=0):
    """Return weighted reference area integral c(z) dA(z), mm2."""
    pitch=geo.pitch;duty=geo.duty;phase=geo.phase
    a,b=float(z1),float(z2)
    nmin=int(floor(a/pitch+phase))-2;nmax=int(ceil(b/pitch+phase))+2
    breaks=[a,b]
    for n in range(nmin,nmax+1):
        breaks.extend([(n-phase)*pitch,(n+duty-phase)*pitch,
                       (n+wet_ramp_fraction*duty-phase)*pitch])
    breaks=np.unique(np.asarray([v for v in breaks if a-1e-12<=v<=b+1e-12]));breaks[0]=a;breaks[-1]=b
    val=0.
    dry=cos(np.deg2rad(thetaY))
    for la,lb in zip(breaks[:-1],breaks[1:]):
        if lb-la<1e-11: continue
        mid=(la+lb)*.5
        f=((mid/pitch+phase)%1.)
        wet=(f>=duty) or (f<duty*wet_ramp_fraction)
        mm=float(geo.dr([mid])[0]);rm=float(geo.radius([mid])[0])
        seg=2*np.pi*rm*np.hypot(1,mm)*(lb-la)
        val+=seg*(1 if wet else dry)
    return float(val)

def local_c(geo,z,thetaY,wet_ramp_fraction,side):
    f=((z+side*1e-8)/geo.pitch+geo.phase)%1.
    return 1. if f>=geo.duty or f<wet_ramp_fraction*geo.duty else np.cos(np.deg2rad(thetaY))

def one_sided_force(geo,z,q,endpoint,thetaY,wet_ramp_fraction):
    sign=1 if endpoint=='front' else -1
    r=float(geo.radius([z])[0]);g=[]
    for side in (-1,1):
        m=float(geo.dr([z+side*1e-8])[0]);c=local_c(geo,z,thetaY,wet_ramp_fraction,side)
        g.append(float(sign*2*np.pi*r*((1+q*m)/np.hypot(1,q)-c*np.hypot(1,m))))
    return g
