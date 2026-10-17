"""Fixed-volume reduced elastocapillary energy for two patterned contact lines.

Static energy and an illustrative isotropic gradient relaxation only.
Neither the fluid mobility nor Cox--Voinov dynamics is derived here.
"""
import numpy as np
from scipy.integrate import solve_ivp
from green_solver import Pattern, equilibria, coefficients

def energy(q, pattern, V=0.2, nu=4., cos_ref=.95):
    xm,xp=np.asarray(q,dtype=float)
    if not (0<xm<xp<1):
        raise ValueError("Contact lines must lie inside the channel")
    I=coefficients(xm,xp)[2]
    k=2*np.pi/pattern.period
    wetting=(pattern.mean*(xp-xm)
             +pattern.amplitude*(np.cos(k*xm+pattern.phase)
                                 -np.cos(k*xp+pattern.phase))/k)
    return (V-(xp-xm))**2/(2*I)-nu*wetting/cos_ref

def gradient(q, pattern, V=0.2, nu=4., cos_ref=.95):
    xm,xp=np.asarray(q,dtype=float)
    km,kp,I,_=coefficients(xm,xp)
    p=(V-(xp-xm))/I
    return np.array([p*(1+p*km)+nu*pattern(xm)/cos_ref,
                     -p*(1+p*kp)-nu*pattern(xp)/cos_ref])

def curvature(q, pattern, h=1e-5, **kwargs):
    q=np.asarray(q,dtype=float)
    J=np.column_stack([(gradient(q+h*np.eye(2)[i],pattern,**kwargs)
                        -gradient(q-h*np.eye(2)[i],pattern,**kwargs))/(2*h)
                       for i in range(2)])
    return (J+J.T)/2,abs(J[0,1]-J[1,0])

def classify(pattern,nseed=42):
    out=[]
    for s in equilibria(pattern,nseed=nseed):
        q=[s['xm'],s['xp']]
        H,skew=curvature(q,pattern)
        eig=np.linalg.eigvalsh(H)
        label='minimum' if min(eig)>1e-5 else 'saddle' if np.prod(eig)<0 else 'unresolved'
        out.append(dict(xm=s['xm'],xp=s['xp'],energy=energy(q,pattern),
                        eigenvalues=eig.tolist(),symmetric_error=skew,
                        classification=label))
    return out

def relax(initial,pattern,final_time=5.,mobility=1.):
    """Illustrative gradient flow: qdot=-mobility*grad(E); NOT original hydrodynamics."""
    return solve_ivp(lambda t,z:-mobility*gradient(z,pattern),[0,final_time],
                     initial,method='BDF',rtol=1e-9,atol=1e-11)

if __name__=='__main__':
    for L in (.25,.15,.10):
        roots=classify(Pattern(period=L))
        print(f'period={L}: {sum(r["classification"]=="minimum" for r in roots)} minima, '
              f'{sum(r["classification"]=="saddle" for r in roots)} saddles')
