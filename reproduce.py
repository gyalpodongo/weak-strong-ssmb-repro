"""Reproduce every number in the note and write results.json.

    python3 reproduce.py                              # all rows, 5e4 particles, 8 processes (about 1 min)
    python3 reproduce.py --only paper_assumptions --particles 100000
"""
from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field, replace
from itertools import repeat

from model import (D_delta, Params, bisect, bunching_gain, emittance_eq20, energy_spread_after_mod1, radiated_power_W,
                   radiator)
from touschek import acceptance_for, lifetime_h, microbunch_enhancement
from tracking import Options, track

DAMPING_TIME_TURNS = 3093
REALISTIC = {"backaction": 0.641, "ibs": True, "radiator_qe": 4.05e-7, "turns": DAMPING_TIME_TURNS}
REALISTIC_SIGMA_DELTA_R = 3.37e-3              # tracked rms spread at the radiator in the realistic rows


def heating_ratio(p: Params) -> float:
    """Rev Mod 1 heating from a path error, kL^2 (A1^2 + 9 A3^2) sigma_Dz^2 / 2, over the (sigma_Dz / r56)^2 / 2
    that Eq. 21 budgets. Above 1 because the beam at Mod 1 spans several laser periods (the sine is not linear)."""
    return (p.kL * p.r56 * p.A1) ** 2 * (1 + 9 * p.A3_over_A1 ** 2)


SIGMA_DZ_MATCHED = Params().sigma_Dz / math.sqrt(heating_ratio(Params()))   # heating equal to the Eq. 21 budget
ENERGY_CONSERVING = D_delta(radiator(Params()).N_u, REALISTIC_SIGMA_DELTA_R)  # beam's coherent loss = Eq. 23 power


@dataclass(frozen=True)
class Row:
    note: str                                     # value quoted in the note
    params: dict = field(default_factory=dict)    # changes to Params
    options: dict = field(default_factory=dict)   # changes to Options


ROWS = {
    "paper_assumptions":       Row("645 W"),
    "paper_sigmaDz_0":         Row("3.2 kW", {"sigma_Dz": 0.0}),
    "paper_sigmaDz_0.091nm":   Row("1.25 kW", {"sigma_Dz": 0.091e-9}),
    "paper_sigmaDz_0.18nm":    Row("0.41 kW", {"sigma_Dz": 0.18e-9}),
    "paper_heating_matched":   Row("0.91 kW", {"sigma_Dz": SIGMA_DZ_MATCHED}),
    "A1_3.6e-3":               Row("(not in note)", {"A1": 3.6e-3}),
    "A1_3.8e-3":               Row("6 W", {"A1": 3.8e-3}),
    "A1_4.1e-3":               Row("37 W", {"A1": 4.1e-3}),
    "A1_4.3e-3":               Row("(not in note)", {"A1": 4.3e-3}),
    "paper_plus_backaction":   Row("35 W", {}, {"backaction": 0.641}),
    "realistic_no_csr_1tau":   Row("36 W", {"sigma_Dz": 0.0}, REALISTIC),
    "realistic_no_ibs_1tau":   Row("(IBS negligible)", {"sigma_Dz": 0.0}, {**REALISTIC, "ibs": False}),
    "realistic_energy_conserving_1tau": Row("70 W", {"sigma_Dz": 0.0}, {**REALISTIC, "backaction": ENERGY_CONSERVING}),
    "realistic_with_csr_1tau": Row("3 W or less", {"sigma_Dz": 0.0}, {**REALISTIC, "csr": True}),
}


def run_row(name: str, particles: int) -> dict:
    row = ROWS[name]
    params = replace(Params(), **row.params)
    options = replace(Options(particles=particles), **row.options)
    return {**track(params, options), "note": row.note}


def gaussian_b1(p: Params, eps_z: float) -> float:
    """|b1| at the FF entrance of a Gaussian beam with emittance eps_z."""
    return math.exp(-p.kL ** 2 * eps_z * p.beta_z0 / 2)


def formulas() -> dict:
    """Paper formulas only (no tracking)."""
    p = Params()
    A1_eq21 = 1 / (p.kL * p.r56)    # Eq. 21 slope h1 = 1 / r56 from the fundamental alone; with A3 the Eq. 12
                                    # amplitude meets it too: (A1 + 3 A3) kL r56 = 0.98
    A1_eq12 = p.A1                  # Eq. 12 optimum: gives the printed b76
    sigma_delta_for_1kW = bisect(lambda s: 1000 - radiated_power_W(p, 0.095, s), 1.5e-3, 6e-3, 80)
    eps_one = emittance_eq20(p, A1_eq12)
    b76_one = bunching_gain(p, A1_eq12) * gaussian_b1(p, eps_one)
    sigma_delta_one = energy_spread_after_mod1(A1_eq12, math.sqrt(eps_one / p.beta_z0), p.A3_over_A1)
    return {
        "A1_eq21": A1_eq21,
        "A1_eq12": A1_eq12,
        "linear_condition_(A1+3A3)_kL_r56": (1 + 3 * p.A3_over_A1) * A1_eq12 * p.kL * p.r56,
        "rev_mod1_heating_over_eq21_budget": heating_ratio(p),
        "sigma_Dz_matched_nm": SIGMA_DZ_MATCHED * 1e9,
        "backaction_energy_conserving": ENERGY_CONSERVING,
        "eps_z_eq20_at_A1_eq21_nm": emittance_eq20(p, A1_eq21) * 1e9,
        "b76_eq12_printed_bunch": bunching_gain(p, A1_eq12) * gaussian_b1(p, 0.34e-9),
        "sigma_delta_eff_for_1kW": sigma_delta_for_1kW,
        "P_printed_b76_ring_spread_W": radiated_power_W(p, 0.095, 1.4e-3),
        "one_amplitude": {"eps_z_nm": eps_one * 1e9, "b76": b76_one,
                          "P_W": radiated_power_W(p, b76_one, sigma_delta_one)},
    }


def touschek() -> dict:
    """Lifetimes of the scan rows above 1 kW (0.9 GeV, 1 nm, 3 % acceptance, uniform line density) and the
    acceptance a 1 h lifetime needs at the Table 1 point."""
    scan_rows = {
        "6.2 kW (0.9 GeV, 1.2 kA peak, ring 71.7 m)": {"E": 0.9e9, "I_peak": 1200, "C_ring": 71.7},
        "17 to 19 kW (0.9 GeV, 2.4 kA peak, ring 56.0 m)": {"E": 0.9e9, "I_peak": 2400, "C_ring": 56.0},
    }
    out = {name: {"lifetime_h": lifetime_h(**ring), "acceptance_for_1h": acceptance_for(1.0, **ring)}
           for name, ring in scan_rows.items()}
    table1 = {"E": 600e6, "I_peak": 50, "C_ring": 311, "beta_low": 0.94, "L_low": 3.9,
              "enhancement": microbunch_enhancement(248e-9, 1.03e-6)}
    out["Table 1 point: acceptance for 1 h at ring beta 5, 10, 20 m"] = [
        acceptance_for(1.0, beta_ring=beta, **table1) for beta in (5.0, 10.0, 20.0)]
    return out


def print_tracking(tracking: dict) -> None:
    print(f"{'row':33s} {'P (Eq.23)':>10s} {'P via b1':>10s} {'b76':>7s} {'floor':>7s} {'eps_z nm':>8s} "
          f"{'sz nm':>6s}  note")
    for name, r in tracking.items():
        print(f"{name:33s} {r['power_W']:10.1f} {r['power_from_b1_W']:10.1f} {r['b_n']:7.4f} "
              f"{r['noise_floor_b_n']:7.4f} {r['eps_z'] * 1e9:8.3f} {r['sigma_z'] * 1e9:6.0f}  {r['note']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--particles", type=int, default=50_000)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--only", default="", help="comma-separated row names (default: all)")
    args = parser.parse_args()
    names = [name for name in ROWS if not args.only or name in args.only.split(",")]

    results = {"formulas": formulas(), "touschek": touschek()}
    with ProcessPoolExecutor(args.jobs) as pool:
        results["tracking"] = dict(zip(names, pool.map(run_row, names, repeat(args.particles))))
    with open("results.json", "w") as f:
        json.dump(results, f, indent=1)

    print(json.dumps(results["formulas"], indent=1))
    print(json.dumps(results["touschek"], indent=1))
    print_tracking(results["tracking"])


if __name__ == "__main__":
    main()
