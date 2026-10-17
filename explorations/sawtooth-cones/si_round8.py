"""S4-informed two-contact-line screening; *necessary* conditions only.

Omer et al. (2025) SI Fig. S4: beta=79+-5 degrees relative to VERTICAL.
Finite reset width is UNMEASURED and scanned (not inferred from Eq. S1).
Uses exact fixed-endpoint CMC BVP; no actual film kinetics or snap prediction.
Lengths in mm; volumes in mm^3; 0.003 mm^3 = 3 nL.
Run: python si_round8.py
"""
from itertools import product
from pathlib import Path
import csv
import numpy as np
from yl_solver import Sawtooth, solve_pinned
from patch_model import one_sided_force


def geometry(pitch_um, beta_deg, reset_fraction, alpha_deg=3):
    pitch=pitch_um/1000
    duty=1-reset_fraction
    slope=np.tan(np.deg2rad(alpha_deg))
    amplitude=duty*pitch*(slope+1/np.tan(np.deg2rad(beta_deg)))
    return Sawtooth(r0=.030,slope=slope,pitch=pitch,
                    amplitude=amplitude,duty=duty,direction=-1)


def planar_ratio(g):
    l=sum((b-a)*np.hypot(1,m) for a,b,ra,rb,m in g._facets(.8,1.2))
    return l/(.4*np.hypot(1,g.slope))


def si_s1(beta_deg,alpha_deg=3):
    a,b=np.deg2rad([alpha_deg,beta_deg])
    return np.cos(a)*(1+np.cos(b)+np.tan(a))/np.sin(b)


def corners(g,target,tolerance):
    c=round(target/g.pitch)
    p=[(n*g.pitch,'crest') for n in range(c-1,c+2)]
    p += [((n+g.duty)*g.pitch,'valley') for n in range(c-1,c+2)]
    return [(x,name) for x,name in p if abs(x-target)<=tolerance]


def test_surface(g,z1,z2,V):
    if not .11<=z2-z1<=.30:return None
    s=solve_pinned(g,z1,z2,V,n=65,tol=3e-5,max_nodes=5500)
    if not s['success']:
        return {'converged':False,'admissible':False,'reason':'bvp_fail'}
    v=[z1,z2]
    for a,b,_,_,_ in g._facets(z1,z2):v.extend([a,b])
    z=np.unique(np.r_[s['z'],v])
    gap=float(np.min(s['sol'].sol((z-z1)/(z2-z1))[0]-g.radius(z)))
    eps=1e-9
    # A free surface above an obstacle must depart the left endpoint
    # no slower than the solid, and arrive at the right no faster.
    rl=bool(float(s['q'][0])>=float(g.dr(z1+eps))-1e-9)
    fr=bool(float(s['q'][-1])<=float(g.dr(z2-eps))+1e-9)
    admissible=gap>=-1e-6
    return dict(converged=True,admissible=admissible,
        reason='valid' if admissible else 'tooth_penetration',
        min_gap_um=1000*gap,local_viable=rl and fr,
        rear_local_viable=rl,front_local_viable=fr,
        q_rear=float(s['q'][0]),q_front=float(s['q'][-1]))


def run():
    shapes=[];states=[];signs=[]
    for pitch,beta,reset in product((10,20,40),(74,79,84),(.01,.03)):
        g=geometry(pitch,beta,reset)
        shapes.append(dict(pitch_um=pitch,beta_deg=beta,reset_fraction=reset,
            height_um=1000*g.amplitude,planar_ratio=planar_ratio(g),
            ideal_SI_eq_S1=si_s1(beta)))
        for (z1,k1),(z2,k2),nl in product(corners(g,.8,.055),corners(g,1.,.065),
                                          (2.5,3.,3.5)):
            if not .11<=z2-z1<=.30:continue
            key=dict(pitch_um=pitch,beta_deg=beta,reset_fraction=reset,
                rear_um=round(1000*z1,6),front_um=round(1000*z2,6),
                rear_type=k1,front_type=k2,volume_nL=nl)
            try:s=test_surface(g,z1,z2,nl/1000)
            except (ValueError,OverflowError,ZeroDivisionError,FloatingPointError):
                s={'converged':False,'admissible':False,'reason':'solver_exception'}
            if s is None:continue
            states.append(dict(**key,**s))
            if not s['admissible']:continue
            for theta,wet in product((35.,52.),(0.,.5,.9)):
                rear=one_sided_force(g,z1,s['q_rear'],'rear',theta,wet)
                front=one_sided_force(g,z2,s['q_front'],'front',theta,wet)
                rear_pass=rear[0]<=1e-9<=rear[1]
                front_pass=front[0]<=1e-9<=front[1]
                signs.append(dict(**key,theta_proxy_deg=theta,wet_ramp_fraction=wet,
                    rear_sign_pass=rear_pass,front_sign_pass=front_pass,
                    both_sign_pass=rear_pass and front_pass))
    return shapes,states,signs


def save(path,rows):
    path.parent.mkdir(exist_ok=True)
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(dict.fromkeys(k for row in rows for k in row)))
        w.writeheader();w.writerows(rows)


def main():
    shapes,states,signs=run()
    root=Path(__file__).resolve().parent/'results'
    for name,rows in [('si_shapes.csv',shapes),('si_states.csv',states),('si_gibbs.csv',signs)]:
        save(root/name,rows)
    converged=[s for s in states if s['converged']]
    admitted=[s for s in converged if s['admissible']]
    penetrated=[s for s in converged if not s['admissible']]
    print('geometries',len(shapes),'candidate_states',len(states),
          'converged',len(converged),'nonpenetrating',len(admitted),
          'penetrating',len(penetrated))
    print('local slope rejects',sum(not s['local_viable'] for s in penetrated),
          'of',len(penetrated),'penetrations')
    print('joint Gibbs passes',sum(s['both_sign_pass'] for s in signs),
          'of',len(signs),'tested parameter combinations')
    assert not any((not s['local_viable']) for s in admitted)
    assert sum(s['both_sign_pass'] for s in signs)==0
    assert len(admitted)>900
    return shapes,states,signs

if __name__=='__main__':main()
