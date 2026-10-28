"""Implicit sensitivity of a simple two-line fold to pattern amplitude and phase."""
from dataclasses import replace
import numpy as np
from beam_fold import Pattern,find_fold,fold_residual


def fold_sensitivity(seed=(.485,.715,.214),pattern=Pattern(),step=2e-6):
    """Return fold and local d(a,b,V)/d(amplitude,phase).

    Differentiate the three stationary/fold equations:
    F_z dz/dtheta + F_theta = 0. This is strictly a local sensitivity
    of one identified fold branch, not a guarantee of global continuity.
    """
    fold=find_fold(seed,pattern)
    z=np.array([fold['a'],fold['b'],fold['volume']]);eye=np.eye(3)
    Jz=np.column_stack([(fold_residual(z+step*e,pattern)-fold_residual(z-step*e,pattern))/(2*step)
                         for e in eye])
    out={}
    for field in ('amplitude','phase'):
        value=getattr(pattern,field)
        plus=replace(pattern,**{field:value+step})
        minus=replace(pattern,**{field:value-step})
        Jtheta=(fold_residual(z,plus)-fold_residual(z,minus))/(2*step)
        derivative=-np.linalg.solve(Jz,Jtheta)
        out[field]=derivative
    return z,out

if __name__=='__main__':
    z,s=fold_sensitivity()
    print('fold (a,b,V):',z)
    for k,v in s.items():print('d(a,b,V)/d'+k,v)
