"""Conditional linear barrier-collapse law on an ideal wet/dry sawtooth patch.

The prescribed wet fraction is NOT measured. The barrier is a constrained
energy difference from wet-film edge to cusp, not a spontaneous jump prediction.
"""
import numpy as np
from thermodynamic_checks import Setup, illustrative_cone, state, GAMMA_PJ_PER_MM2
from yl_solver import solve_pinned

def barrier_rows(theta_y=35.0):
    geo = illustrative_cone()
    valley = 1.019  # mm, cusp after one 20 um pitch
    rear, volume = 0.820, 0.003
    sol = solve_pinned(geo, rear, valley, volume, n=110, tol=1e-7)
    assert sol['success'] and sol['min_gap'] > -1e-6
    q = float(sol['q'][-1])
    r = float(geo.radius(valley))
    # Gentle facet immediately BEFORE cusp, not the steep wet reset face.
    slope = float(geo.dr(valley - 1e-9))
    c = np.cos(np.deg2rad(theta_y))
    slope_pJ_per_um = (
        GAMMA_PJ_PER_MM2 * 2 * np.pi * r
        * ((1+q*slope)/np.hypot(1,q) - c*np.hypot(1,slope))
        / 1000
    )
    for wet in [.90, .95, .99, .995, .999]:
        setup = Setup(theta_dry_deg=theta_y, wet_ramp_fraction=wet)
        x_edge = 1 + .019 * wet
        E_edge = state(geo, setup, x_edge, volume)['energy_pJ']
        E_valley = state(geo, setup, valley, volume)['energy_pJ']
        width_um = (valley-x_edge)*1000
        actual = E_valley-E_edge
        linear = width_um*slope_pJ_per_um
        yield dict(wet_fraction=wet, remaining_dry_um=width_um,
                   numerical_barrier_pJ=actual, linear_pJ=linear,
                   ratio=actual/linear)

def run():
    for theta in [35, 52]:
        rows=list(barrier_rows(theta))
        for r in rows:
            print(f"theta={theta}, f={r['wet_fraction']:.3f}, "
                  f"barrier={r['numerical_barrier_pJ']:.6f} pJ, "
                  f"linear={r['linear_pJ']:.6f} pJ, ratio={r['ratio']:.6f}")
        assert abs(rows[-1]['ratio']-1) < .001, rows[-1]
        # Numerical noise dominates the last 19 nm: no strict monotonicity claim.
        assert abs(rows[2]['ratio']-1) < abs(rows[0]['ratio']-1)
    print('PASS: sharp-edge linear barrier scaling in two dry-angle proxy cases')

if __name__=='__main__':
    run()
