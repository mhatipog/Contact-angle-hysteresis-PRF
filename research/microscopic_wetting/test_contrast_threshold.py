"""Local contrast threshold and analytical saddle-node barrier tests."""
import unittest
import numpy as np
from beam_fold import Pattern, state
from contrast_threshold import (critical_contrast, normal_form,
                                barrier_above_contrast)


class TestContrastThreshold(unittest.TestCase):
    def test_simple_critical_amplitude_fold(self):
        location, A, soft=critical_contrast()
        self.assertAlmostEqual(A, .047020818924159, delta=1e-9)
        self.assertAlmostEqual(location[0], .48196748, delta=2e-7)
        self.assertAlmostEqual(location[1], .71958230, delta=2e-7)
        _,g,H,_=state(*location,.22,Pattern(amplitude=A))
        self.assertLess(np.linalg.norm(g),1e-8)
        eig=np.linalg.eigvalsh(H)
        self.assertLess(abs(eig[0]),1e-6)
        self.assertGreater(eig[1],300)

    def test_normal_form_matches_independent_minimum_saddle_barriers(self):
        nf=normal_form()
        self.assertAlmostEqual(nf["barrier_coefficient"],.41759484,delta=2e-5)
        offsets=np.array([1e-4,5e-5,2e-5,1e-5,5e-6,2e-6,1e-6])
        barriers=[]
        for delta in offsets:
            result=barrier_above_contrast(float(delta))
            self.assertGreater(result["barrier"],0)
            self.assertGreater(np.linalg.norm(result["stable"]-result["saddle"]),1e-5)
            barriers.append(result["barrier"])
        slope=np.polyfit(np.log(offsets),np.log(barriers),1)[0]
        self.assertGreater(slope,1.498)
        self.assertLess(slope,1.502)
        self.assertLess(abs(barriers[-1]/(nf["barrier_coefficient"]*offsets[-1]**1.5)-1),.002)

    def test_negative_offset_rejected(self):
        with self.assertRaises(ValueError):
            barrier_above_contrast(-.0001)


if __name__=="__main__":
    unittest.main()
