"""One-pitch thermodynamic control and two-contact-line Gibbs sign tests.

Independent exploratory extension (not part of Bradley et al. 2021).
This is NOT a predictive wet-film dynamics model. Units: mm, mm^3, pJ.
Run: python periodic_no_go.py
"""
from pathlib import Path
import csv
import numpy as np
from yl_solver import Sawtooth, solve_pinned
from patch_model import energy_weighted, one_sided_force

GAMMA_PJ_MM2 = 72000.0
VOLUME_MM3 = .003
PITCH_MM = .020
REAR, FRONT = .816, 1.016


def admissible_energy(geo, shift, wet=.9, theta=35.):
    rear, front = REAR + shift, FRONT + shift
    s = solve_pinned(geo, rear, front, VOLUME_MM3, n=80, tol=3e-6)
    if not s["success"]:
        raise RuntimeError("Young-Laplace collocation did not converge")
    # Check every polygon vertex in addition to the uniform surface grid.
    vertices = [rear, front]
    for a, b, *_ in geo._facets(rear, front):
        vertices += [a, b]
    z = np.unique(np.r_[s["z"], vertices])
    gap = s["sol"].sol((z - rear) / (front - rear))[0] - geo.radius(z)
    if np.min(gap) < -1e-6:
        return None  # Inadmissible: DO NOT assign a physical energy.
    ref = energy_weighted(geo, rear, front, theta, wet)
    return GAMMA_PJ_MM2 * (s["area_lv"] - ref)


def periodic_control(samples=21):
    rows = []
    for direction in (1, -1):
        for wet in (0., .9):
            geo = Sawtooth(r0=.08, slope=0, pitch=PITCH_MM,
                           amplitude=.003, duty=.8, direction=direction)
            phases = np.linspace(0, PITCH_MM, samples)
            energies = [admissible_energy(geo, float(t), wet)
                        for t in phases]
            assert all(e is not None for e in energies)
            drift = energies[-1] - energies[0]
            span = max(energies) - min(energies)
            assert abs(drift) < .015, ("period closure", drift)
            assert span > 1., ("no pinning landscape?", span)
            print("cylinder: direction=%+d wet=%g span=%.4f pJ, drift=%+.8f pJ"
                  % (direction, wet, span, drift))
            rows.extend(dict(direction=direction, wet_fraction=wet,
                             phase_um=1000*t, energy_pJ=e)
                        for t, e in zip(phases, energies))
    return rows


def joint_corner_signs():
    # Independent illustrative state: rear at crest, front at valley.
    # Volume, footprint, geometry, and wet reference are NOT fitted data.
    geo = Sawtooth(r0=.08, slope=0, pitch=PITCH_MM,
                   amplitude=.003, duty=.8)
    rear, front = .816, 1.020
    s = solve_pinned(geo, rear, front, VOLUME_MM3, n=110, tol=1e-7)
    assert s["success"] and s["min_gap"] >= -1e-8
    rows = []
    for theta in (35., 52.):
        for wet in (0., .5, .9, 1.):
            f_rear = one_sided_force(geo, rear, float(s["q"][0]),
                                     "rear", theta, wet)
            f_front = one_sided_force(geo, front, float(s["q"][-1]),
                                      "front", theta, wet)
            # A necessary first-order corner-minimum condition for BOTH
            # independent endpoint coordinates is E'_- <= 0 <= E'_+.
            ok_rear = bool(f_rear[0] < 0 < f_rear[1])
            ok_front = bool(f_front[0] < 0 < f_front[1])
            rows.append(dict(theta_proxy_deg=theta, wet_fraction=wet,
                             rear_left=f_rear[0], rear_right=f_rear[1],
                             front_left=f_front[0], front_right=f_front[1],
                             rear_sign_pass=ok_rear,
                             front_sign_pass=ok_front,
                             both_sign_pass=ok_rear and ok_front))
            print("joint corner: theta=%g wet=%g rear=%s front=%s"
                  % (theta, wet, ok_rear, ok_front))
    assert next(r for r in rows if r["theta_proxy_deg"] == 52.
                and r["wet_fraction"] == 0.)["both_sign_pass"]
    assert not next(r for r in rows if r["theta_proxy_deg"] == 52.
                    and r["wet_fraction"] == .5)["both_sign_pass"]
    return rows


def main():
    periodic = periodic_control()
    gibbs = joint_corner_signs()
    out = Path(__file__).parent / "results"
    out.mkdir(exist_ok=True)
    for file, rows in (("periodic_energy.csv", periodic),
                       ("gibbs_corners.csv", gibbs)):
        with (out / file).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    print("PASS: periodic closure and two-line necessary Gibbs sign tests")


if __name__ == "__main__":
    main()
