"""Translation symmetry / obstacle diagnostics in pinned Young--Laplace barrels.

At fixed liquid volume and footprint length, a dry circular cylinder is
translation invariant: E(z_-,z_+)=E(z_-+t,z_++t). Taper breaks that symmetry.
A formal CMC BVP can still cross a sawtooth obstacle: such a state is rejected,
not assigned a physical energy. This is a diagnostic, NOT a motion law.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from thermodynamic_checks import illustrative_cone, GAMMA_PJ_PER_MM2
from yl_solver import Sawtooth, solve_pinned

RHO_VOLUME_MM3 = .003
REAR_MM, FRONT_MM = .816, 1.016
THETA_DRY_DEG = 35.0

def pinned_translation_energy(geo, shift_mm, quadrature=140):
    rear, front = REAR_MM + shift_mm, FRONT_MM + shift_mm
    sol = solve_pinned(geo, rear, front, RHO_VOLUME_MM3, n=75,
                       tol=3e-6, max_nodes=6000)
    # Collocation convergence is NOT geometric admissibility.
    if not sol['success'] or sol['min_gap'] < -1e-6:
        return None
    nodes,weights=leggauss(quadrature)
    uu=(nodes+1)/2
    y=sol['sol'].sol(uu)
    area=2*np.pi*(front-rear)*np.dot(weights/2, y[0]*np.hypot(1,y[1]))
    ref=np.cos(np.deg2rad(THETA_DRY_DEG))*geo.wetted_area(rear,front)
    return GAMMA_PJ_PER_MM2*(area-ref)

def run():
    teeth=illustrative_cone()
    cases={
       'cylinder':Sawtooth(r0=.08,slope=0,amplitude=0),
       'smooth cone':Sawtooth(r0=.03,slope=teeth.slope,amplitude=0),
       'sawtooth cone':teeth,
    }
    shifts=np.linspace(0,.019,20)
    rows={name:[pinned_translation_energy(geo,t) for t in shifts]
          for name,geo in cases.items()}
    cyl=np.asarray(rows['cylinder'])
    cone=np.asarray(rows['smooth cone'])
    tooth=rows['sawtooth cone']
    cylinder_spread=float(np.ptp(cyl))
    smooth_delta=float(cone[-1]-cone[0])
    rejected=sum(v is None for v in tooth)
    print(f'cylinder translation energy spread: {cylinder_spread:.3e} pJ')
    print(f'smooth-cone translation energy change over 19 um: {smooth_delta:.3f} pJ')
    print(f'sawtooth states rejected for geometric penetration/nonconvergence: {rejected}/{len(shifts)}')
    assert cylinder_spread < 1e-4
    assert np.all(np.diff(cone)<0)
    assert 0<rejected<len(shifts)
    print('PASS: Noether-type translation check and nonpenetration screening')

if __name__=='__main__':
    run()
