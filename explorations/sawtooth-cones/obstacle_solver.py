"""Discrete axisymmetric capillary obstacle minimizer (exploratory).

At prescribed endpoint positions and excess volume, minimize the frustum
surface area with p >= r(z), including polygon tooth vertices in the mesh.
A fully wetted/film-reference geometric limit, NOT a dry contact-line
Young-angle model, contact dynamics, or proof of global minimization.
Units: mm, mm2, mm3.
"""
import numpy as np
from scipy.optimize import minimize, brentq


def mesh(geo,z1,z2,max_dx=.002):
    points=[float(z1),float(z2)]
    for a,b,_,_,_ in geo._facets(z1,z2):points.extend([a,b])
    initial=np.unique(np.array(points))
    res=[initial[0]]
    for a,b in zip(initial[:-1],initial[1:]):
        res.extend(np.linspace(a,b,max(2,int(np.ceil((b-a)/max_dx))+1))[1:])
    z=np.array(res)
    return z,geo.radius(z)


def volume_and_grad(p,z,r):
    dx=np.diff(z)
    a,b=p[:-1],p[1:]
    ra,rb=r[:-1],r[1:]
    v=np.pi/3*np.sum(dx*(a*a+a*b+b*b-ra*ra-ra*rb-rb*rb))
    g=np.zeros_like(p)
    g[:-1]+=np.pi/3*dx*(2*a+b)
    g[1:]+=np.pi/3*dx*(a+2*b)
    return v,g


def area_and_grad(p,z):
    dx=np.diff(z)
    d=np.diff(p)
    s=np.hypot(dx,d)
    pr=p[:-1]+p[1:]
    A=np.pi*np.sum(pr*s)
    g=np.zeros_like(p)
    g[:-1]+=np.pi*(s-pr*d/s)
    g[1:]+=np.pi*(s+pr*d/s)
    return A,g


def solve_obstacle(geo,z1,z2,V,max_dx=.002,maxiter=600,ftol=1e-11):
    z,r=mesh(geo,z1,z2,max_dx)
    u=(z-z1)/(z2-z1)
    def vol(a):
        return volume_and_grad(r+a*np.sin(np.pi*u),z,r)[0]-V
    amp=brentq(vol,0,.5)
    initial=r+amp*np.sin(np.pi*u)
    def profile(x):
        return np.r_[r[0],x,r[-1]]
    def objective(x):
        a,grad=area_and_grad(profile(x),z)
        return 100*a,100*grad[1:-1]
    def constraint(x):
        return volume_and_grad(profile(x),z,r)[0]/V-1
    def jac(x):
        return volume_and_grad(profile(x),z,r)[1][1:-1]/V
    sol=minimize(objective,initial[1:-1],jac=True,
        bounds=[(float(x),None) for x in r[1:-1]],
        constraints={'type':'eq','fun':constraint,'jac':jac},
        method='SLSQP',
        options={'maxiter':maxiter,'ftol':ftol,'disp':False})
    p=profile(sol.x)
    vol,_=volume_and_grad(p,z,r)
    area,_=area_and_grad(p,z)
    gap=p-r
    touches=(gap[1:-1]<.000025)
    return dict(success=bool(sol.success),message=sol.message,
                iterations=sol.nit,nodes=len(z),area=area,volume=vol,
                volume_error_nL=(vol-V)*1000,
                min_gap_um=float(np.min(gap))*1000,
                contacts=int(np.sum(touches)),z=z,r=r,p=p)
