"""Axisymmetric constant-mean-curvature drop with prescribed contact line positions.

Dimensions mm. Surface p(z) is solved from Young-Laplace ODE,
    kappa = 1/(p sqrt(1+p_z**2)) - p_zz/(1+p_z**2)**(3/2).
The volume constraint is pi integral[p(z)**2 - solid_radius(z)**2] dz = V.
Both contact lines are PINNED: their apparent angles are OUTPUTS, not imposed.
A pinned state is not automatically stable or a prediction of depinning.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_bvp, quad, trapezoid

@dataclass(frozen=True)
class Sawtooth:
    r0: float=.030
    slope: float=.05
    pitch: float=.020
    amplitude: float=.003
    duty: float=.80
    phase: float=.0
    direction:int=1
    def radius(self,z):
        x=np.asarray(z,dtype=float)
        if self.amplitude==0: return self.r0+self.slope*x
        f=np.mod(x/self.pitch+self.phase,1.)
        v=np.where(f<self.duty, f/self.duty, (1-f)/(1-self.duty))
        if self.direction<0:
            v=1-v
        return self.r0+self.slope*x+self.amplitude*v
    def dr(self,z):
        x=np.asarray(z,dtype=float)
        if self.amplitude==0:return 0*x+self.slope
        f=np.mod(x/self.pitch+self.phase,1.)
        v=np.where(f<self.duty, 1/self.duty, -1/(1-self.duty))
        return self.slope+self.direction*self.amplitude/self.pitch*v
    def rim_position(self,n):
        return (n+self.duty-self.phase)*self.pitch
    def _facets(self,z1,z2):
        a,b=float(z1),float(z2)
        if b<=a: return []
        x=(np.arange(np.floor(a/self.pitch+self.phase)-2, np.ceil(b/self.pitch+self.phase)+3))
        breaks=np.sort(np.unique(np.r_[a,b,(x-self.phase)*self.pitch,(x+self.duty-self.phase)*self.pitch]))
        breaks=breaks[(breaks>=a-1e-12)&(breaks<=b+1e-12)]
        breaks[0]=a;breaks[-1]=b
        segs=[]
        for la,lb in zip(breaks[:-1],breaks[1:]):
            if lb-la<=1e-12:continue
            mid=(la+lb)/2
            slope=float(self.dr(mid))
            rm=float(self.radius(mid))
            ra=rm-slope*(mid-la);rb=rm+slope*(lb-mid)
            segs.append((la,lb,ra,rb,slope))
        return segs
    def sol_volume(self,z1,z2):
        return float(np.pi*sum((b-a)*(ra*ra+ra*rb+rb*rb)/3 for a,b,ra,rb,m in self._facets(z1,z2)))
    def wetted_area(self,z1,z2):
        return float(2*np.pi*sum((b-a)*.5*(ra+rb)*np.hypot(1,m) for a,b,ra,rb,m in self._facets(z1,z2)))

def solve_pinned(geo,z1,z2,V,n=100,tol=1e-5,init=None,max_nodes=15000):
    """Return solution of exact fixed-endpoint axisymmetric CMC BVP."""
    L=z2-z1
    if L<=0 or V<=0:raise ValueError('Bad domain/volume')
    u=np.linspace(0,1,n)
    r1,r2=map(float,geo.radius([z1,z2]))
    rbase=r1+(r2-r1)*u
    H=max(.002,2*V/(np.pi*L*(r1+r2)*1.273))
    pinit=rbase+H*np.sin(np.pi*u)
    qinit=(r2-r1)/L+H*np.pi/L*np.cos(np.pi*u)
    z=z1+L*u
    vini=np.r_[0, np.cumsum(.5*np.diff(u)*(np.pi*L*(pinit[1:]**2-geo.radius(z[1:])**2 + pinit[:-1]**2-geo.radius(z[:-1])**2)))]
    if vini[-1]>0:vini*=V/vini[-1]
    if init is not None:
        try:
            pinit=np.interp(u,init.x,init.y[0]);qinit=np.interp(u,init.x,init.y[1]);vini=np.interp(u,init.x,init.y[2])*V/max(init.y[2,-1],1e-10)
        except Exception:pass
    def fun(u,y,k):
        p,q,v=y
        pp=np.maximum(p,1e-9)
        w=np.sqrt(1+q*q)
        return np.vstack((L*q, L*((1+q*q)/pp - k[0]*w**3), np.pi*L*(p*p-geo.radius(z1+L*u)**2)))
    def bc(a,b,k):
        return np.array([a[0]-r1,b[0]-r2,a[2],b[2]-V])
    guess=np.vstack((pinit,qinit,vini))
    k0=np.array([1/np.max(pinit)])
    sol=solve_bvp(fun,bc,u,guess,p=k0,tol=tol,max_nodes=max_nodes)
    uu=np.linspace(0,1,1001);yy=sol.sol(uu)
    p,q=yy[:2];zz=z1+L*uu;solid=geo.radius(zz)
    theta_left=np.degrees(np.arctan2(q[0],1)-np.arctan(geo.dr(z1-1e-8)))
    theta_right=np.degrees(np.arctan(geo.dr(z2+1e-8))-np.arctan2(q[-1],1))
    A=2*np.pi*L*trapezoid(p*np.sqrt(1+q*q),uu)
    q_u=sol.sol(uu,1)[1]
    kappa=1/(p*np.sqrt(1+q*q))-(q_u/L)/(1+q*q)**1.5
    return dict(success=bool(sol.success),message=sol.message,z1=z1,z2=z2,V=V,p=p,q=q,z=zz,
        theta_left=theta_left,theta_right=theta_right,kappa=float(sol.p[0]),
        kappa_spread=float(np.max(np.abs(kappa-sol.p[0]))),
        min_gap=float(np.min(p-solid)),min_p=float(np.min(p)),
        volume=float(sol.y[2,-1]),area_lv=float(A), area_sl=geo.wetted_area(z1,z2),
        sol=sol,shape_r=solid)

if __name__=='__main__':
    geo=Sawtooth()
    for V in [2e-4,5e-4,1e-3,2e-3]:
        try:
            s=solve_pinned(geo,geo.rim_position(40),geo.rim_position(52),V,n=80)
            print(V,s['success'],s['theta_left'],s['theta_right'],s['min_gap'],s['kappa_spread'],s['message'][:55])
        except Exception as e: print('error',V,str(e))
