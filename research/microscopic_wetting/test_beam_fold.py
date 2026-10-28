"""Independent numerical checks for energy derivatives and coupled folds."""
import unittest
import numpy as np
from scipy.integrate import dblquad
from scipy.optimize import root
from beam_fold import (Pattern, beam_compliance, state, find_fold,
                       effective_left_curvature,gradient_flow_energy_balance)


def raw_energy(a,b,V,p):
    d=b-a
    I=d**5/120+d*d*((a+b)/2)**3/3
    return (V-d)**2/(2*I)-p.nu*p.integral(a,b)


class TestBeamFold(unittest.TestCase):
    def setUp(self):
        self.p = Pattern()

    def test_cantilever_compliance_against_independent_quadrature(self):
        for a,b in [(0.15,0.37),(0.39,0.61),(0.51,0.73)]:
            I=beam_compliance(a,b)[0]
            val=dblquad(lambda y,x:min(x,y)**2*(3*max(x,y)-min(x,y))/6,
                        a,b,lambda x:a,lambda x:b,epsabs=1e-10)[0]
            self.assertAlmostEqual(I,val,delta=2e-12)

    def test_exact_gradient_and_hessian_against_separate_differentiation(self):
        for a,b,V in [(.39,.63,.2),(.45,.68,.2),(.51,.74,.21),(.3,.58,.23)]:
            _,grad,H,pressure=state(a,b,V,self.p)
            h=1e-27
            independent=np.array([
                np.imag(raw_energy(a+1j*h,b,V,self.p))/h,
                np.imag(raw_energy(a,b+1j*h,V,self.p))/h])
            self.assertLess(np.linalg.norm(grad-independent),1e-10)
            e=2e-6
            for j in range(2):
                v=np.eye(2)[j]*e
                d=(state(*(np.array([a,b])+v),V,self.p)[1]
                   -state(*(np.array([a,b])-v),V,self.p)[1])/(2*e)
                self.assertLess(np.max(abs(H[:,j]-d)),2e-5)
            self.assertAlmostEqual(pressure,(V-b+a)/beam_compliance(a,b)[0])
            self.assertLess(abs(H[0,1]-H[1,0]),1e-12)

    def test_young_cosine_admissibility(self):
        Pattern()
        with self.assertRaises(ValueError):
            Pattern(amplitude=.25,cosine_reference=.9)
        with self.assertRaises(ValueError):
            Pattern(wavelength=0)
        with self.assertRaises(ValueError):
            Pattern(amplitude=1.1,cosine_reference=.1)

    def test_coupled_fold_and_adjacent_minimum_saddle(self):
        fold=find_fold((.485,.715,.214),self.p)
        self.assertAlmostEqual(fold['volume'],.213851157109733,delta=2e-9)
        self.assertAlmostEqual(fold['a'],.485098866979384,delta=2e-9)
        self.assertAlmostEqual(fold['b'],.715416013259478,delta=2e-9)
        self.assertLess(abs(fold['eigenvalues'][0]),1e-6)
        self.assertGreater(fold['eigenvalues'][1],300)
        self.assertLess(np.linalg.norm(fold['residual']),1e-7)
        # A nearby volume has two coexisting stationary configurations.
        minima=[]
        for seed in [(.482,.713),(.488,.719)]:
            sol=root(lambda z:state(*z,.214,self.p)[1],seed,tol=1e-11)
            self.assertTrue(sol.success)
            self.assertLess(np.linalg.norm(sol.fun),1e-8)
            H=state(*sol.x,.214,self.p)[2]
            minima.append((sol.x,np.linalg.eigvalsh(H),H))
        self.assertGreater(minima[0][1][0],0)
        self.assertLess(minima[1][1][0],0)
        # Holding right meniscus fixed overestimates soft curvature ~50x.
        H=minima[0][2]
        self.assertGreater(H[0,0]/effective_left_curvature(H),40)

    def test_pressure_work_and_exact_dissipation_identity(self):
        for a,b,V in [(.4,.64,.2),(.46,.71,.22),(.55,.8,.23)]:
            for rate in [-.003,0.,.011]:
                v, power, diss, work=gradient_flow_energy_balance(
                    a,b,V,self.p,(.2,2.),rate)
                self.assertGreaterEqual(diss,0)
                self.assertLess(abs(power-work+diss),1e-12)
                # Check the pressure conjugate to V independently.
                z=1e-26
                EV=np.imag(raw_energy(a,b,V+1j*z,self.p))/z
                self.assertAlmostEqual(work,EV*rate,delta=1e-12)
                self.assertEqual(v.shape,(2,))

    def test_no_extrapolation_into_collapsed_channel(self):
        with self.assertRaises(ValueError):
            state(.3,.52,.25,self.p)
        with self.assertRaises(ValueError):
            state(.94,1.04,.09,self.p)


if __name__ == '__main__':
    unittest.main()
