"""Exact static potential vs contact pressure yield in the original quartic beam.

Tests at 25 stationary patterned solutions with complex-step gradients.
Separately tests nonnegative energy rate for a reduced dry-friction
gradient inclusion. This is NOT the energy law of the full dynamic PDE.
"""
import numpy as np
from scipy.optimize import root
V=.2;NU=4.;ELL=.25;XREF=.45

def compliance(a,b):
    w=b-a;q=(a+b)/2
    return w*w*(w**3+40*q**3)/120

def energy(a,b,A):
    k=2*np.pi/ELL
    w=b-a
    W=w+A*(np.cos(k*(a-XREF))-np.cos(k*(b-XREF)))/k
    return (V-w)**2/(2*compliance(a,b))-NU*W

def gradient(a,b,A):
    z=1e-26
    return np.array([np.imag(energy(a+1j*z,b,A))/z,
                     np.imag(energy(a,b+1j*z,A))/z])

def equilibrium(a,A):
    k=2*np.pi/ELL
    def R(y):
        b,p,c,d=y;w=a-b
        hl=p*w**4/24+c*w+d;hp=p*w**3/6+c
        return np.asarray([p*w*w/2-2/a**2*(2*a*hp-3*hl+3),
           p*w-6/a**3*(a*hp-2*hl+2),
           p*d+NU*(1+A*np.sin(k*(b-XREF))),
           -p*w**5/120-c*w*w/2-d*w-V])
    sol=root(R,[a+V+.01,-NU,-.04*NU,1-.02*NU],tol=1e-12)
    b,p,c,d=sol.x
    w=a-b;hl=p*w**4/24+c*w+d
    assert max(abs(R(sol.x)))<1e-9 and a<b<1 and p<0 and d>0 and hl>0
    return b,-p*hl/NU-1

def soft_yield(g,friction=.04,mobility_inv=10):
    v=-np.sign(g)*np.maximum(np.abs(g)-friction,0)/mobility_inv
    diss=friction*np.abs(v).sum()+mobility_inv*(v@v)
    return v,float(g@v),float(diss)

def audit():
    max_error=0.
    for a in [.35,.4,.4454383377,.49,.55]:
      for A in [-.002,-.001,0,.001,.002]:
        b,G=equilibrium(a,A)
        ga,gb=gradient(a,b,A)
        local=A*np.sin(2*np.pi*(a-XREF)/ELL)
        residual=max(abs(gb),abs(ga+NU*(G-local)))
        max_error=max(max_error,residual)
        assert residual<3e-9,(a,A,residual)
        v,power,diss=soft_yield(np.array([ga,gb]))
        assert abs(power+diss)<1e-12 and diss>=0
    rng=np.random.default_rng(271828)
    for _ in range(1000):
      a=rng.uniform(.2,.55);b=rng.uniform(max(a+.205,.48),.85)
      A=rng.uniform(-.005,.005);g=gradient(a,b,A)
      v,power,diss=soft_yield(g,rng.uniform(0,.3),rng.uniform(.1,30))
      assert np.isclose(power,-diss,atol=1e-10,rtol=1e-11)
    print(f'PASS 25 stationary beam yield identities, max gradient residual {max_error:.3g}, plus 1000 exact dissipative gradient flow checks')

if __name__=='__main__':audit()
