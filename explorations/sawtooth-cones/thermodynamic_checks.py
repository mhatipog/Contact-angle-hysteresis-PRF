"""Pressure-force reciprocity for pinned Young--Laplace barrels.

Energy uses actual meniscus area and piecewise polygonal solid area.
Wet-patch geometry, contact angle and cone dimensions are illustrative.
No experimental video events are fitted or predicted.
Units: z/r in mm, V in mm^3, E in pJ, pressure pJ/mm^3,
and moving-front endpoint energy derivative in pJ/mm.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from yl_solver import Sawtooth, solve_pinned
from patch_model import energy_weighted, local_c

GAMMA_N_PER_M = 0.072
GAMMA_PJ_PER_MM2 = GAMMA_N_PER_M * 1e6

@dataclass(frozen=True)
class Setup:
    rear_mm: float = 0.820
    volume_mm3: float = 0.003  # 3 nL
    theta_dry_deg: float = 35.0  # proxy, NOT independently measured Young angle
    wet_ramp_fraction: float = 0.90

def illustrative_cone() -> Sawtooth:
    """20-micrometre-pitch example, not a metrology reconstruction."""
    pitch = 0.020
    taper = np.tan(np.deg2rad(3.0))
    return Sawtooth(
        r0=0.030, slope=taper, pitch=pitch,
        amplitude=pitch * .95 * (taper + np.tan(np.deg2rad(11.0))),
        duty=.95, direction=-1,
    )

def state(geo: Sawtooth, setup: Setup, front_mm: float, volume_mm3: float,
          quadrature: int = 220) -> dict:
    """Return exact-area energy, YL pressure, analytical endpoint energy derivative."""
    sol = solve_pinned(geo, setup.rear_mm, front_mm, volume_mm3,
                       n=95, tol=2e-7, max_nodes=10000)
    if not sol['success'] or sol['min_gap'] < -1e-6:
        raise RuntimeError(f"inadmissible BVP: {sol['message']}; gap={sol['min_gap']}")
    t, w = leggauss(quadrature)
    u = (t + 1.0) / 2.0
    yy = sol['sol'].sol(u)
    area_lv = (2 * np.pi * (front_mm - setup.rear_mm)
               * np.dot(w / 2, yy[0] * np.hypot(1.0, yy[1])))
    weighted = energy_weighted(geo, setup.rear_mm, front_mm,
                               setup.theta_dry_deg, setup.wet_ramp_fraction)
    E_pJ = GAMMA_PJ_PER_MM2 * (area_lv - weighted)
    r, dr = float(geo.radius(front_mm)), float(geo.dr(front_mm + 1e-9))
    q = float(sol['q'][-1])
    c = local_c(geo, front_mm, setup.theta_dry_deg,
                setup.wet_ramp_fraction, side=1)
    endpoint = 2 * np.pi * GAMMA_PJ_PER_MM2 * r * (
        (1 + q * dr) / np.hypot(1, q) - c * np.hypot(1, dr))
    pressure = GAMMA_PJ_PER_MM2 * sol['kappa']
    return dict(energy_pJ=E_pJ, pressure_pJ_mm3=pressure,
                endpoint_pJ_mm=endpoint, min_gap_mm=sol['min_gap'],
                kappa_spread_mm_inverse=sol['kappa_spread'])

def derivatives(geo: Sawtooth, setup: Setup, x: float,
                hx: float = 0.000075, hv: float = 0.00002) -> dict:
    """Finite differences independent of endpoint force and YL pressure.

    Maxwell identity: d/dV (dE/dx) = d/dx (dE/dV).
    Positions must lie INSIDE a smooth facet and one wet/dry segment.
    """
    base = state(geo, setup, x, setup.volume_mm3)
    xp = state(geo, setup, x + hx, setup.volume_mm3)
    xm = state(geo, setup, x - hx, setup.volume_mm3)
    vp = state(geo, setup, x, setup.volume_mm3 + hv)
    vm = state(geo, setup, x, setup.volume_mm3 - hv)
    dxE = (xp['energy_pJ'] - xm['energy_pJ']) / (2 * hx)
    dvE = (vp['energy_pJ'] - vm['energy_pJ']) / (2 * hv)
    dxP = (xp['pressure_pJ_mm3'] - xm['pressure_pJ_mm3']) / (2 * hx)
    dvF = (vp['endpoint_pJ_mm'] - vm['endpoint_pJ_mm']) / (2 * hv)
    return dict(front_mm=x,
                local_surface='wet' if local_c(geo,x,setup.theta_dry_deg,setup.wet_ramp_fraction,1)==1 else 'dry',
                E_pJ=base['energy_pJ'],
                pressure_numeric=dvE, pressure_YL=base['pressure_pJ_mm3'],
                force_numeric=dxE, force_endpoint=base['endpoint_pJ_mm'],
                mixed_dF_dV=dvF, mixed_dP_dx=dxP,
                pressure_relerr=abs(dvE-base['pressure_pJ_mm3'])/abs(base['pressure_pJ_mm3']),
                force_relerr=abs(dxE-base['endpoint_pJ_mm'])/abs(base['endpoint_pJ_mm']),
                maxwell_relerr=abs(dvF-dxP)/max(abs(dxP),abs(dvF)))

def run():
    geo = illustrative_cone()
    setup = Setup()
    for x in [1.0270, 1.0375]:
        d = derivatives(geo, setup, x)
        print(f"{d['local_surface']:3} x={x:.4f} mm "
              f"pressure error={d['pressure_relerr']:.2e}, "
              f"endpoint error={d['force_relerr']:.2e}, "
              f"Maxwell error={d['maxwell_relerr']:.2e}")
        assert d['pressure_relerr'] < 3e-5, d
        assert d['force_relerr'] < 3e-4, d
        assert d['maxwell_relerr'] < 1e-3, d
    print('PASS: conjugacy, endpoint force, and mixed derivatives, both film regimes')

if __name__ == '__main__':
    run()
