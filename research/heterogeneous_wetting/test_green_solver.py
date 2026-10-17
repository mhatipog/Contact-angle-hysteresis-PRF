"""Sanity checks for static patterned bendotaxis, not a physical dynamical test."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).parent))
from green_solver import Pattern, coefficients, candidate, equilibria

def test_original_stationary_pressure_and_asymmetry():
    xm=.45;xp=.6609023191263425;V=.2;nu=4.
    km,kp,I,ke=coefficients(xm,xp)
    p=(V-(xp-xm))/I
    assert np.isclose(p,-4.284951129848285,atol=1e-8)
    assert np.isclose((1+p*km)/(1+p*kp)-1,.0314960597357,atol=1e-9)
    assert abs(p*(1+p*kp)+nu)<1e-8

def test_beam_green_volume_integral():
    xm=.45;xp=.66
    km,kp,I,ke=coefficients(xm,xp)
    x=np.linspace(xm,xp,501)
    X,S=np.meshgrid(x,x,indexing='ij')
    G=np.minimum(X,S)**2*(3*np.maximum(X,S)-np.minimum(X,S))/6
    check=np.trapezoid(np.trapezoid(G,x,axis=0),x,axis=0)
    assert np.isclose(check,I,rtol=2e-6)
    assert 0<km<kp<ke

def test_local_young_volume_pressure():
    results=equilibria(Pattern(period=.15),nseed=42)
    assert len(results)==6
    for s in results:
        assert s['mismatch']<1e-8 and s['h1']>0
        assert np.isclose(s['lambda_geometry'],s['lambda_chemistry'],atol=2e-8)
        I=coefficients(s['xm'],s['xp'])[2]
        assert np.isclose(s['xp']-s['xm']+s['p']*I,.2,atol=2e-8)

def test_search_resolution_wavelength_contrast():
    for n in (24,42):
        assert len(equilibria(Pattern(period=.15),nseed=n))==6
        assert len(equilibria(Pattern(period=.10),nseed=n))==0

def test_phase_checks_and_homogeneous_case():
    for phase in (0,np.pi/2,np.pi,3*np.pi/2):
        assert len(equilibria(Pattern(period=.1,phase=phase),nseed=45))==0
    assert len(equilibria(Pattern(period=.15,amplitude=0),nseed=24))==0

def test_two_meniscus_phase_identity():
    A=.025;L=.15;phase=.7;xm=.31;xp=.53
    pat=Pattern(amplitude=A,period=L,phase=phase)
    difference=2*A*np.sin(np.pi*(xm-xp)/L)*np.cos(np.pi*(xm+xp)/L+phase)
    assert np.isclose(pat(xm)-pat(xp),difference,atol=1e-14)
