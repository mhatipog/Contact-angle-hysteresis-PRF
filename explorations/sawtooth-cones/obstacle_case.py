"""Reproduce one SI-oriented CMC penetration and its obstacle correction.

This is fixed-volume capillary geometry, NOT a prediction of snapping,
finite film kinetics, or a globally certified energy minimizer.
Run: python obstacle_case.py
"""
import numpy as np
from si_round8 import geometry
from yl_solver import solve_pinned
from obstacle_solver import solve_obstacle


def run():
    g=geometry(20,79,.01)
    z1=(40+g.duty)*g.pitch   # tooth valley
    z2=51*g.pitch           # tooth crest
    V=.003                  # 3 nL
    cmc=solve_pinned(g,z1,z2,V,n=100,tol=2e-6)
    assert cmc['success']
    verts=[z1,z2]
    for a,b,_,_,_ in g._facets(z1,z2):verts.extend([a,b])
    zz=np.unique(np.r_[cmc['z'],verts])
    pen_um=float(np.min(cmc['sol'].sol((zz-z1)/(z2-z1))[0]-g.radius(zz))*1000)
    assert pen_um < -4.
    results=[solve_obstacle(g,z1,z2,V,max_dx=dx)
             for dx in (.004,.002,.001)]
    for s in results:
        assert s['success'],s['message']
        assert s['min_gap_um']>=-1e-5
        assert abs(s['volume_error_nL'])<1e-7
        assert s['contacts']>=1
    assert abs(results[-1]['area']-results[-2]['area'])/results[-1]['area']<1e-5
    print('Unconstrained CMC solid penetration: %.6f um' % pen_um)
    for dx,s in zip((.004,.002,.001),results):
        print('obstacle dx=%g mm: N=%d, area=%.12f mm2, volume_error=%.2e nL, contacts=%d'
             %(dx,s['nodes'],s['area'],s['volume_error_nL'],s['contacts']))
    print('PASS: nonpenetrating, volume conserving, mesh-convergent example')


if __name__=='__main__':
    run()
