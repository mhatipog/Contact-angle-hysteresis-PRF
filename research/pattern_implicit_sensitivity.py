"""Complex-step implicit derivative of quartic clamped/free equilibrium.

No time-dependent integrations needed: a single 4x4 Jacobian inversion
predicts the critical patterned-wetting amplitude.
"""
import numpy as np
from scipy.optimize import root,brentq
XM=.4454383377;V=.2;NU=4.;ELL=.25;PIN=.0293
def R(q,A):
    xu,p,c,d=q;w=XM-xu;hl=p*w**4/24+c*w+d;hp=p*w**3/6+c
    return np.asarray([p*w*w/2-2/XM**2*(2*XM*hp-3*hl+3),
                       p*w-6/XM**3*(XM*hp-2*hl+2),
                       p*d+NU*(1+A*np.sin(2*np.pi*(xu-XM)/ELL)),
                       -p*w**5/120-c*w*w/2-d*w-V])
def G(q):
    xu,p,c,d=q;w=XM-xu
    return -p*(p*w**4/24+c*w+d)/NU-1
def solve(A):
    z=root(lambda q:R(q,A),[XM+V+.01,-NU,-.04*NU,1-.02*NU],tol=1e-12)
    assert np.max(np.abs(R(z.x,A)))<1e-9
    return z.x
def result():
    q=solve(0); h=1e-28
    J=np.column_stack([np.imag(R(q.astype(complex)+1j*h*np.eye(4)[j],0))/h for j in range(4)])
    RA=np.imag(R(q,1j*h))/h
    dq=-np.linalg.solve(J,RA)
    Gq=np.array([np.imag(G(q.astype(complex)+1j*h*np.eye(4)[j]))/h for j in range(4)])
    slope=float(Gq@dq)
    linear=(PIN-G(q))/slope
    exact=brentq(lambda a:G(solve(a))-PIN,0,.004,xtol=1e-14)
    return slope,float(linear),float(exact)
if __name__=='__main__':
    slope,linear,exact=result()
    print(f'G_A={slope:.12g} linear_Acrit={linear:.12g} nonlinear_Acrit={exact:.12g} relative_error={abs(linear-exact)/exact:.9g}')
    assert abs(slope+.896031545419)<1e-6
    assert abs(linear-exact)/exact<.0003
    assert G(solve(.002))<PIN<G(solve(-.002))
    print('PASS complex-step 4x4 implicit gradient predicts phase transition')
