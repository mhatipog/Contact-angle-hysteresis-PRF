"""Steady beam approximation compared with deposited Figure 11 dynamic thresholds.

Stationary X_plus is NOT the same physical observable as the initial escape
threshold x_plus^0,escape. Numerical differences are expected; these are
diagnostics, not an accusation of erroneous published data.
"""
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.optimize import root, brentq

def equilibrate(xminus, V, nu):
    def residual(a):
        xp,p,c,d=a
        w=xminus-xp
        hl=p*w**4/24+c*w+d
        hp=p*w**3/6+c
        return [p*w*w/2-2/xminus**2*(2*xminus*hp-3*hl+3),
                p*w-6/xminus**3*(xminus*hp-2*hl+2),
                p*d+nu,
                -p*w**5/120-c*w*w/2-d*w-V]
    q=root(residual,[xminus+V+.02,-nu,-.04*nu,1-.02*nu],tol=1e-11)
    xp,p,c,d=q.x
    if not(max(abs(np.asarray(residual(q.x))))<1e-8 and
           xminus<xp<1 and p<0 and d>0 and
           c*(1-xp)+d>0):raise ValueError('No admissible equilibrium')
    w=xminus-xp
    hl=p*w**4/24+c*w+d
    if min(p*(x-xp)**4/24+c*(x-xp)+d
           for x in np.linspace(xminus,xp,20))<=0:
        raise ValueError('Nonpositive wet gap')
    return xp,hl/d-1.

def predicted_xplus(V,nu,lam):
    scan=[]
    for xminus in np.linspace(.04,.99-V,80):
        try:scan.append((xminus,equilibrate(xminus,V,nu)[1]))
        except ValueError:pass
    for (a,fa),(b,fb) in zip(scan,scan[1:]):
        if (fa-lam)*(fb-lam)<=0:
            x=brentq(lambda z:equilibrate(z,V,nu)[1]-lam,a,b,xtol=1e-12)
            return equilibrate(x,V,nu)[0]
    return None

def audit(path=Path('Data/numerics_comparison_data.mat')):
    d=loadmat(path); V=float(d['V'].item()); rows=[]
    for i,lam in enumerate(d['hysparmax_out'].ravel()):
        for j,nu in enumerate(d['nu_out'].ravel()):
            if nu not in (2,3,4,5,6):continue
            xp=predicted_xplus(V,float(nu),float(lam))
            if xp is None:continue
            ref=float(d['xu0_escape'][i,j])
            rows.append((round(float(lam),2),float(nu),ref,xp,ref-xp))
    return rows

if __name__=='__main__':
    rows=audit()
    for lam,nu,observed,stationary,delta in rows:
        print(f'lambda={lam:.2f} nu={nu:.1f} archive={observed:.8f} stationary={stationary:.8f} archive-minus-steady={delta:+.8f}')
    assert len(rows)==22
    assert all(delta>.025 for lam,nu,_,_,delta in rows if lam==.02)
    assert all(abs(delta)<.003 for lam,nu,_,_,delta in rows if lam==.04)
    assert all(delta<-.006 for lam,nu,_,_,delta in rows if lam>=.06)
    print('PASS: 22 physical equilibrium comparisons; sign change independently reproduced.')
