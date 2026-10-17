from dataclasses import replace
import numpy as np
from scipy.integrate import quad
from model import (Surface,line_force,rigid_threshold,rigid_velocity,
                   relaxed_length,relaxed_energy,relaxed_force,relaxed_threshold,
                   asymptotic_threshold,spectral_force,force_rms_from_covariance,
                   two_harmonic_design,quartic_velocity,generic_velocity)

def test_endpoint_force_and_exact_cancellation():
    s=Surface();q=.31;d=.27
    direct=2*s.gamma*s.amplitude*(np.sin(s.k*(q+d/2))-np.sin(s.k*(q-d/2)))
    assert abs(line_force(q,d,s)-direct)<1e-12
    assert rigid_threshold(.4,s)<1e-14
    assert abs(rigid_threshold(.3,s)-.1)<1e-13

def test_rigid_average_velocity_against_quadrature():
    s=Surface();d=.3;fc=rigid_threshold(d,s)
    for F in [fc*1.01,fc*1.25,fc*2.]:
        T=quad(lambda q:1/(F+line_force(q,d,s)),0,s.wavelength,
               epsrel=1e-10,limit=300)[0]
        assert np.isclose(rigid_velocity(F,d,s),s.wavelength/T,rtol=2e-9)

def test_relaxed_energy_gradient():
    s=Surface()
    for q in np.linspace(.11,.51,13):
        h=1e-5
        numerical=-(relaxed_energy(q+h,s)-relaxed_energy(q-h,s))/(2*h)
        assert abs(numerical-relaxed_force(q,s))<1e-8

def test_residual_second_order_pinning():
    s=Surface()
    assert rigid_threshold(s.reference_length,s)<1e-13
    actual=relaxed_threshold(s)
    assert actual>0
    assert abs(actual/asymptotic_threshold(s)-1)<.001

def test_stiffness_and_amplitude_scalings():
    s=Surface()
    K=np.array([20.,40.,80.,160.])
    vals=np.array([relaxed_threshold(replace(s,stiffness=float(k))) for k in K])
    assert abs(np.polyfit(np.log(K),np.log(vals),1)[0]+1)<.02
    A=np.array([.006,.012,.024])
    vals=np.array([relaxed_threshold(replace(s,amplitude=float(a))) for a in A])
    assert abs(np.polyfit(np.log(A),np.log(vals),1)[0]-2)<.02

def test_general_fourier_filter():
    amps=[.025,-.006,.004];phases=[.12,.67,-.3]
    q=.4;ell=.2;d=.17
    direct=sum(2*a*(np.sin((n+1)*2*np.pi*(q+d/2)/ell+p)-
                    np.sin((n+1)*2*np.pi*(q-d/2)/ell+p))
               for n,(a,p) in enumerate(zip(amps,phases)))
    assert abs(spectral_force(q,d,ell,amps,phases)-direct)<1e-13
    assert abs(spectral_force(q,ell,ell,amps,phases))<1e-13

def test_wettability_correlation_identity():
    x=np.linspace(0,2*np.pi,4096,endpoint=False)
    c=.85+.02*np.sin(x+.2)+.005*np.cos(3*x+.4)
    for n in [0,15,421,1024,2048,4096]:
        pred,actual=force_rms_from_covariance(c,n)
        assert abs(pred-actual)<1e-12

def test_quartic_wettability_design_and_exponent():
    (a1,a2),B=two_harmonic_design()
    assert abs(a2+ a1/(8*np.cos(np.pi/4)))<1e-14
    assert np.isclose(4*a2*np.sin(np.pi/2),-B/4)
    eps=np.array([1e-8,3e-8,1e-7,3e-7])
    u=np.array([quartic_velocity(e,B) for e in eps])
    generic=np.array([generic_velocity(e,B) for e in eps])
    assert abs(np.polyfit(np.log(eps),np.log(u),1)[0]-.75)<.01
    assert abs(np.polyfit(np.log(eps),np.log(generic),1)[0]-.5)<.001
