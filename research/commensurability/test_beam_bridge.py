import numpy as np
from scipy.integrate import dblquad
from beam_bridge import (Beam,green_integral,beam_energy,stationary_length,
                         effective_stiffness,commensurate_experiment)

def test_green_integral_with_independent_quadrature():
    for a,b in [(.15,.37),(.39,.61),(.51,.73)]:
        I=green_integral(a,b)
        numeric=dblquad(lambda s,x:min(x,s)**2*(3*max(x,s)-min(x,s))/6,
                        a,b,lambda x:a,lambda x:b,epsabs=1e-10)[0]
        assert abs(I-numeric)<1e-10

def test_stiffness_from_unpatterned_beam():
    b=Beam();q=.6
    d,E=stationary_length(q,b,lo=.20001,hi=.45)
    K=effective_stiffness(q,d,b)
    assert .2<d<.3 and K>0
    assert beam_energy(q,d+1e-4,b)>E
    assert beam_energy(q,d-1e-4,b)>E

def test_second_order_survives_real_beam_reduction():
    samples=commensurate_experiment()
    for A,d0,d,K,observed,predicted,pressure in samples:
        assert 0<d<1 and pressure<0 and observed<0 and K>0
        assert abs(observed/predicted-1)<.005
    a=np.array([r[0] for r in samples])
    shifts=-np.array([r[4] for r in samples])
    assert abs(np.polyfit(np.log(a),np.log(shifts),1)[0]-2)<.005
