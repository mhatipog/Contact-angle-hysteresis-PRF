"""Exploratory *stationary* elastocapillary equilibria with local Young angles.

Based on Hatipogullari et al. PRFluids (2019) and Bradley et al. (2021).
No dynamics, energy stability, or trapping probability is implied.
"""
from dataclasses import dataclass
from math import pi
import numpy as np
from scipy.optimize import root

@dataclass(frozen=True)
class Pattern:
    mean: float = .95
    amplitude: float = .025
    period: float = .15
    phase: float = 0.
    def __call__(self,x):
        return self.mean + self.amplitude*np.sin(2*pi*x/self.period + self.phase)

def coefficients(xm,xp):
    # Exact integrals of clamped-free cantilever G(x,s):
    # G = min(x,s)^2 * (3*max(x,s)-min(x,s)) / 6.
    km=xm*xm/6*(1.5*(xp*xp-xm*xm)-xm*(xp-xm))
    kp=(.75*xp**4-xp*xm**3+.25*xm**4)/6
    integral=(xp**5-xm**5)/20-xm**3*(xp*xp-xm*xm)/6+xm**4*(xp-xm)/12
    kend=(xp**3-xm**3)/6-(xp**4-xm**4)/24
    return km,kp,integral,kend

def candidate(xm,xp,pattern,V=.2,nu=4.,cos_ref=.95):
    if not(.03<xm<xp<.99 and xm+V<1):return None
    km,kp,I,kend=coefficients(xm,xp)
    if I<1e-10:return None
    p=(V-(xp-xm))/I
    hm,hp,h1=[1+p*x for x in (km,kp,kend)]
    if not(p<0 and min(hm,hp,h1)>0):return None
    cm,cp=pattern(xm),pattern(xp)
    if not(0<cm<=1 and 0<cp<=1):return None
    em=(p*hm+nu*cm/cos_ref)/nu
    ep=(p*hp+nu*cp/cos_ref)/nu
    return dict(xm=float(xm),xp=float(xp),p=float(p),hm=float(hm),
                hp=float(hp),h1=float(h1),eminus=float(em),eplus=float(ep),
                mismatch=float(max(abs(em),abs(ep))),
                lambda_geometry=float(hm/hp-1),
                lambda_chemistry=float(cm/cp-1))

def residual(q,pattern,V,nu,cos_ref):
    xm,xp=q
    c=candidate(xm,xp,pattern,V,nu,cos_ref)
    if c is not None:return [c['eplus'],c['eminus']]
    if not(.01<xm<xp<1):return [10*(xm-.5),10*(xp-xm-.25)]
    km,kp,I,kend=coefficients(xm,xp)
    if I<=1e-9:return [10.,10.]
    p=(V-(xp-xm))/I
    return [p*(1+p*kp)/nu+pattern(xp)/cos_ref,
            p*(1+p*km)/nu+pattern(xm)/cos_ref]

def equilibria(pattern,V=.2,nu=4.,cos_ref=.95,nseed=24):
    """Multistart numerical search; NOT proven exhaustive."""
    results=[]
    for xm0 in np.linspace(.05,.75,nseed):
        for width in np.linspace(V+.005,.35,nseed):
            if xm0+width>=.98:continue
            solution=root(residual,[xm0,xm0+width],args=(pattern,V,nu,cos_ref),tol=1e-9)
            xm,xp=solution.x
            c=candidate(xm,xp,pattern,V,nu,cos_ref)
            if c is not None and c['mismatch']<1e-8 and not any(
                    abs(c['xm']-z['xm'])<1e-5 and abs(c['xp']-z['xp'])<1e-5
                    for z in results):
                results.append(c)
    return sorted(results,key=lambda z:z['xm'])

if __name__=="__main__":
    for L in (.25,.15,.10):
        roots=equilibria(Pattern(period=L),nseed=42)
        print(f"wavelength={L:.2f}: {len(roots)} candidate stationary states")
