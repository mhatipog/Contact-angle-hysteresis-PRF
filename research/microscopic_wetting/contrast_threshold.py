"""Critical chemical contrast and 3/2 energy-barrier normal form.

One conservative Young-cosine surface landscape and the exact reduced
clamped/free beam energy. This is LOCAL deterministic bifurcation theory,
not a universal measured contact-angle hysteresis law.
"""
from dataclasses import replace
import numpy as np
from scipy.optimize import root
from beam_fold import Pattern, state


def critical_contrast(volume=.22, pattern=Pattern(),
                      guess=(.482, .72, .047)):
    """Find A where both stationary contact lines lose local stability."""
    def equations(z):
        a, b, amplitude = z
        try:
            local = replace(pattern, amplitude=float(amplitude))
            _, grad, H, _ = state(a, b, volume, local)
            return np.array([grad[0], grad[1], np.linalg.det(H)/100])
        except ValueError:
            return np.ones(3)*100
    solution = root(equations, guess, tol=1e-11)
    if np.linalg.norm(equations(solution.x))>1e-8:
        raise RuntimeError("No verified stable/saddle contrast fold")
    a,b,A = map(float, solution.x)
    local = replace(pattern, amplitude=A)
    _, _, H, _ = state(a,b,volume,local)
    eigenvalues, modes = np.linalg.eigh(H)
    if not(abs(eigenvalues[0])<1e-6 and eigenvalues[1]>0):
        raise RuntimeError("The contrast root is not a simple minimum/saddle fold")
    return np.array([a,b]), A, modes[:,0]


def normal_form(volume=.22, pattern=Pattern(), h=1e-5):
    """Predict barrier coefficient from soft-mode force and third derivative.

    At the fold, E_soft ~ s*dA*x + t*x^3/6 and E_soft'=0.
    Thus Delta E ~ C*dA^(3/2), C=(4/3)|s|sqrt(2|s|/|t|).
    """
    z, Ac, e = critical_contrast(volume,pattern)
    p = replace(pattern, amplitude=Ac)
    a,b=z; k=p.wave_number
    s = p.nu*(e[0]*np.sin(k*a+p.phase)
              -e[1]*np.sin(k*b+p.phase))
    Hp = state(*(z+h*e),volume,p)[2]
    Hm = state(*(z-h*e),volume,p)[2]
    t = float(e@((Hp-Hm)/(2*h))@e)
    if s*t >= 0:
        raise RuntimeError("Expected two stationary solutions for amplitude above fold")
    L = float(np.sqrt(-2*s/t))
    coefficient = 4*abs(s)*L/3
    return dict(position=z, amplitude=Ac, soft_mode=e,
                force_derivative=float(s), cubic_derivative=t,
                displacement_coefficient=L, barrier_coefficient=coefficient)


def barrier_above_contrast(delta, volume=.22, pattern=Pattern()):
    """Find local minimum/saddle pair and their energy barrier at Ac+delta."""
    if not delta>0:raise ValueError("delta must be strictly positive")
    nf = normal_form(volume,pattern)
    A = nf["amplitude"]+delta
    local=replace(pattern,amplitude=A)
    z=nf["position"]; e=nf["soft_mode"]
    shift=nf["displacement_coefficient"]*np.sqrt(delta)*e

    def grad(y):
        try:return state(*y,volume,local)[1]
        except ValueError:return np.ones(2)*100
    solutions=[]
    for sign in (-1,1):
        res=root(grad,z+sign*shift,tol=1e-11)
        if np.linalg.norm(res.fun)>1e-8:
            raise RuntimeError("Minimum or saddle did not converge")
        energy,_,H,_=state(*res.x,volume,local)
        solutions.append((res.x,float(energy),float(np.linalg.eigvalsh(H)[0])))
    stable=[x for x in solutions if x[2]>0]
    saddles=[x for x in solutions if x[2]<0]
    if len(stable)!=1 or len(saddles)!=1:
        raise RuntimeError("No distinct stable/saddle pair")
    return dict(delta=float(delta),contrast=float(A),stable=stable[0][0],
                saddle=saddles[0][0],barrier=float(saddles[0][1]-stable[0][1]),
                leading_prediction=float(nf["barrier_coefficient"]*delta**1.5))


if __name__=="__main__":
    nf=normal_form()
    print("critical contrast:",nf["amplitude"],"location:",nf["position"])
    print("normal form: s =",nf["force_derivative"],
          "third derivative =",nf["cubic_derivative"],
          "barrier coefficient =",nf["barrier_coefficient"])
    offsets=np.array([1e-4,5e-5,2e-5,1e-5,5e-6,2e-6,1e-6])
    rows=[barrier_above_contrast(float(delta)) for delta in offsets]
    for row in rows:
        print("offset=",row["delta"],"barrier=",row["barrier"],
              "leading order=",row["leading_prediction"])
    slope=np.polyfit(np.log(offsets),np.log([x["barrier"] for x in rows]),1)[0]
    print("observed exponent:",slope,"(generic saddle-node target 1.5)")
