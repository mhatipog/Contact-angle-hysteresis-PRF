"""Quasi-static volume cycling of two contact lines in a patterned elastic channel.

Reduced, energy-consistent benchmark. Beam shape follows the exact clamped/
free Green function, but moving-fluid hydrodynamics and the microscopic
contact-line mobility are NOT solved. Quasi-static metastability is followed
locally; globally re-minimizing at every volume would erase pinning.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from beam_fold import Pattern, state


@dataclass
class Jump:
    volume: float
    before: tuple
    after: tuple
    energy_drop: float
    integrated_dissipation: float
    relaxation_time: float
    smallest_free_tip_gap: float


def stable_local_equilibrium(z, volume, pattern=Pattern(), max_displacement=.03):
    """Accept only a nearby stable root, never a remote/global minimizer."""
    def gradient(y):
        try:
            return state(y[0], y[1], volume, pattern)[1]
        except ValueError:
            return np.array([1e2, 1e2])
    result = root(gradient, z, tol=1e-11)
    if np.linalg.norm(result.fun) > 3e-9:
        return None
    try:
        _, _, H, _ = state(*result.x, volume, pattern)
    except ValueError:
        return None
    if np.linalg.eigvalsh(H)[0] < 1e-6 or np.linalg.norm(result.x-z) > max_displacement:
        return None
    return result.x


def local_relaxation(z, volume, pattern=Pattern(), mobility=(1., 1.),
                     gradient_tolerance=2e-9, max_time=80., max_step=.06):
    """Energy-decreasing jump with an explicit wall-escape and energy check."""
    M = np.asarray(mobility, float)
    if M.shape != (2,) or not np.all(np.isfinite(M)) or np.any(M <= 0):
        raise ValueError("Mobility must contain two finite positive values")

    def rhs(t, y):
        # Invalid trial geometry raises, rather than silently clipping a wall.
        return -M*state(y[0], y[1], volume, pattern)[1]

    def settled(t, y):
        return np.linalg.norm(state(y[0], y[1], volume, pattern)[1])-gradient_tolerance
    settled.terminal = True
    settled.direction = -1

    def outlet(t, y):
        return .9999-y[1]
    outlet.terminal = True
    outlet.direction = -1

    sol = solve_ivp(rhs, (0., max_time), np.asarray(z, float),
                    method="BDF", rtol=2e-10, atol=1e-12, max_step=max_step,
                    events=[settled, outlet])
    if sol.status != 1 or len(sol.t_events[0]) != 1:
        raise RuntimeError("Outlet escape or unresolved transient: not a repinning event")

    energies = []
    powers = []
    tip_gaps = []
    for time, y in zip(sol.t, sol.y.T):
        E, g, _, p = state(*y, volume, pattern)
        energies.append(E)
        powers.append(np.dot(M*g, g))
        a, b = y
        K_tip = (b**3-a**3)/6-(b**4-a**4)/24
        tip_gaps.append(1+p*K_tip)
    energies = np.asarray(energies)
    dissipated = float(np.trapezoid(powers, sol.t))
    lost = float(energies[0]-energies[-1])
    if np.max(np.diff(energies)) > 1e-9:
        raise AssertionError("Autonomous gradient flow increased energy")
    if abs(lost-dissipated) > max(3e-5, dissipated*.0005):
        raise AssertionError("Jump energy loss is not accounted for by dissipation")
    if min(tip_gaps) <= 0:
        raise AssertionError("Elastic gap collapsed")
    return sol.y[:, -1], Jump(float(volume), tuple(map(float, z)),
                              tuple(map(float, sol.y[:, -1])),
                              lost, dissipated, float(sol.t[-1]),
                              float(min(tip_gaps)))


def volume_cycle(volume_high=.22, volume_low=.203, steps_per_leg=136,
                 pattern=Pattern(), start=(.46881481443305373,
                                          .7066099342611455),
                 mobility=(1., 1.)):
    """Adiabatically change V down and back, with local relaxation at folds."""
    if not 0 < volume_low < volume_high < 1 or steps_per_leg < 2:
        raise ValueError("Invalid cycle")
    volumes = np.r_[np.linspace(volume_high, volume_low, steps_per_leg+1),
                    np.linspace(volume_low, volume_high, steps_per_leg+1)[1:]]
    z = np.array(start, float)
    records, jumps = [], []
    for volume in volumes:
        candidate = stable_local_equilibrium(z, volume, pattern)
        if candidate is None:
            z, jump = local_relaxation(z, volume, pattern, mobility=mobility)
            jumps.append(jump)
        else:
            z = candidate
        E, g, H, pressure = state(*z, volume, pattern)
        records.append((volume, z[0], z[1], E, pressure,
                        np.linalg.eigvalsh(H)[0]))

    records = np.asarray(records)
    volume_work = float(np.trapezoid(records[:, 4], records[:, 0]))
    delta_energy = float(records[-1, 3]-records[0, 3])
    dissipated = volume_work-delta_energy
    return records, jumps, volume_work, delta_energy, dissipated


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=136)
    parser.add_argument("--output", default=None, help="Optional CSV destination")
    args = parser.parse_args()
    rows, jumps, work, change, diss = volume_cycle(steps_per_leg=args.steps)
    for idx, jump in enumerate(jumps, 1):
        print(f"release/repin {idx}: {jump}")
    print("initial:", rows[0])
    print("final:", rows[-1])
    print(f"volume work={work:.12g}, delta energy={change:.12g}, "
          f"irreversible loss={diss:.12g}")
    if args.output:
        np.savetxt(args.output, rows, delimiter=",",
                   header="volume,left,right,energy,pressure,smallest_hessian_eigenvalue",
                   comments="")
