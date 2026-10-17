"""Tests of a reduced stationary energy, not published trapping predictions."""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).parent))
from green_solver import Pattern
from energy_stability import energy,gradient,curvature,classify,relax

def test_variational_contact_line_equations():
    for L in (.15,.25):
        pat=Pattern(period=L)
        for row in classify(pat,nseed=24):
            q=np.array([row['xm'],row['xp']]);h=2e-5
            num=np.array([(energy(q+h*np.eye(2)[i],pat)
                           -energy(q-h*np.eye(2)[i],pat))/(2*h) for i in range(2)])
            np.testing.assert_allclose(num,gradient(q,pat),atol=2e-5,rtol=1e-6)
            assert np.linalg.norm(gradient(q,pat))<1e-7

def test_hessian_stability_and_symmetry():
    for L,expected_min,expected_saddle in [(.25,1,1),(.15,3,3),(.10,0,0)]:
        rows=classify(Pattern(period=L),nseed=42)
        assert sum(r['classification']=='minimum' for r in rows)==expected_min
        assert sum(r['classification']=='saddle' for r in rows)==expected_saddle
        for r in rows:
            assert r['symmetric_error']<1e-3
            assert min(abs(x) for x in r['eigenvalues'])>1

def test_gradient_flow_returns_to_each_local_minimum():
    pat=Pattern(period=.15)
    minima=[r for r in classify(pat,nseed=42) if r['classification']=='minimum']
    assert len(minima)==3
    for row in minima:
        q0=np.array([row['xm'],row['xp']]);qinit=q0+np.array([.004,.004])
        sol=relax(qinit,pat,final_time=5.)
        assert sol.success and np.linalg.norm(sol.y[:,-1]-q0)<1e-6
        values=np.array([energy(y,pat) for y in sol.y.T])
        assert np.all(np.diff(values)<1e-10)
        assert values[-1]<values[0]
