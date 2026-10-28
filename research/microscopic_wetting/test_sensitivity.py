"""Independent continuation checks of the fold's parameter derivatives."""
import unittest
from dataclasses import replace
import numpy as np
from beam_fold import Pattern, find_fold, state
from sensitivity import fold_sensitivity

class TestSensitivity(unittest.TestCase):
    def test_implicit_derivatives_against_nonlocal_refolding(self):
        p=Pattern()
        z,dz=fold_sensitivity(pattern=p)
        for parameter,epsilon in [('amplitude',1e-4),('phase',1e-3)]:
            plus=replace(p,**{parameter:getattr(p,parameter)+epsilon})
            minus=replace(p,**{parameter:getattr(p,parameter)-epsilon})
            fplus=find_fold(z,plus)
            fminus=find_fold(z,minus)
            fd=(np.array([fplus[k] for k in ('a','b','volume')])
                -np.array([fminus[k] for k in ('a','b','volume')]))/(2*epsilon)
            self.assertLess(np.max(abs(fd-dz[parameter])),1e-5)
        self.assertAlmostEqual(dz['amplitude'][2],-.04829899003444,delta=1e-7)
        self.assertAlmostEqual(dz['phase'][2],.000821908396257,delta=1e-8)

    def test_same_fold_survives_wetting_perturbations(self):
        for A in [.04,.06,.08,.1]:
            p=Pattern(amplitude=A)
            f=find_fold((.48,.72,.22),p)
            self.assertTrue(0<f['volume']<f['b']-f['a'])
            self.assertLess(abs(f['eigenvalues'][0]),1e-5)
            self.assertGreater(f['eigenvalues'][1],0)
            self.assertAlmostEqual(float(np.linalg.norm(state(
                f['a'],f['b'],f['volume'],p)[1])),0,delta=2e-7)

if __name__=='__main__':unittest.main()
