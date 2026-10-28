"""Reproducibility gates for quasi-static release and repinning."""
import unittest
import numpy as np
from beam_fold import Pattern, state, find_fold
from quasistatic_cycle import volume_cycle, local_relaxation


class TestVolumeCycle(unittest.TestCase):
    def test_three_successive_stability_thresholds(self):
        p = Pattern()
        seeds = [
            ((.485, .715, .214), .213851157109733),
            ((.559, .792, .2075), .2076472661283195),
            ((.63, .868, .1987), .19869833385466745),
        ]
        for seed, expected in seeds:
            fold = find_fold(seed, p)
            self.assertAlmostEqual(fold["volume"], expected, delta=1e-9)
            self.assertLess(abs(fold["eigenvalues"][0]), 1e-6)
            self.assertGreater(fold["eigenvalues"][1], 0)
            self.assertLess(np.linalg.norm(fold["residual"]), 1e-7)

    def test_local_cycle_has_two_repin_events_and_dissipates(self):
        S, jumps, work, dE, dissipated = volume_cycle(steps_per_leg=68)
        self.assertEqual(len(jumps), 2)
        self.assertTrue(.2135 < jumps[0].volume < .21385115711)
        self.assertTrue(.2074 < jumps[1].volume < .20764726613)
        self.assertGreater(S[-1, 1]-S[0, 1], .13)
        self.assertGreater(S[-1, 2]-S[0, 2], .16)
        self.assertAlmostEqual(S[0, 0], S[-1, 0])
        self.assertAlmostEqual(dissipated, .043658760845, delta=1e-5)
        self.assertGreater(dissipated, 0)
        self.assertAlmostEqual(work, dE+dissipated, delta=1e-12)
        for jump in jumps:
            self.assertGreater(jump.energy_drop, 0)
            self.assertGreater(jump.smallest_free_tip_gap, .7)
            self.assertAlmostEqual(jump.energy_drop, jump.integrated_dissipation,
                                   delta=3e-5)
        # Compare aggregate jump loss to the volume-cycled balance; slow
        # relaxation at discrete volume steps makes an O(delta V) correction.
        self.assertAlmostEqual(sum(x.energy_drop for x in jumps), dissipated,
                               delta=1.2e-4)

    def test_resolution_and_one_time_not_continuous_pump(self):
        low, events, _, _, D = volume_cycle(steps_per_leg=68)
        hi, finer, _, _, D2 = volume_cycle(steps_per_leg=136)
        self.assertEqual(len(finer), 2)
        self.assertLess(abs(D-D2), 1e-5)
        self.assertLess(np.max(abs(low[-1, 1:3]-hi[-1, 1:3])), 1e-7)
        S2, events2, _, change2, D3 = volume_cycle(
            steps_per_leg=68, start=tuple(low[-1, 1:3]))
        self.assertEqual(len(events2), 0)
        self.assertLess(np.linalg.norm(S2[-1, 1:3]-S2[0, 1:3]), 1e-8)
        self.assertLess(abs(change2), 1e-8)
        self.assertLess(abs(D3), 1e-7)

    def test_repins_are_robust_to_positive_diagonal_mobility(self):
        # Mobility is NOT derived from the fluid equations; this tests
        # qualitative basin selection under plausible diagonal variations.
        for V, old, expected in [
            (.21375, (.4838941378587488, .7143155265340685),
             (.5433844964578693, .7843039740239497)),
            (.207625, (.5563545036821276, .7892556676551393),
             (.6144045628031181, .8640811319087841)),
        ]:
            for M in ((.2, 1.), (1., 5.)):
                end, jump = local_relaxation(old, V, mobility=M,
                                             max_time=100)
                self.assertLess(np.linalg.norm(end-expected), 3e-7)
                self.assertGreater(jump.energy_drop, 0)

    def test_homogeneous_no_go_for_interior_center_stationarity(self):
        # Fixed width d, varying q=(a+b)/2: for A=0, the Young term is
        # translation invariant, while I_q=d^2*q^2>0. Thus E_q<0
        # for a suction-driven deformable beam. Patterning is necessary
        # for the interior metastability seen in this particular model.
        p=Pattern(amplitude=0.)
        for V in (.19, .22):
            for q in (.36, .52, .69):
                d=.24
                a,b=q-d/2,q+d/2
                E,g,H,pressure=state(a,b,V,p)
                directional=float(np.dot(g,[1,1]))
                self.assertLess(directional, 0)
                delta=V-d
                I=d**5/120+d*d*q**3/3
                predicted=-delta**2*d*d*q*q/(2*I**2)
                self.assertAlmostEqual(directional,predicted,delta=2e-12)


if __name__=="__main__":
    unittest.main()
