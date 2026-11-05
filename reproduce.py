"""Reproduce every number in the note. Run: python reproduce.py [--particles 50000] [--jobs 8]
Writes results.json and prints each number next to the value quoted in the note."""
from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

from ws_repro.model import Params, P_eq23_W, S_bbhg, eps_eq20, sde_after_mod1
from ws_repro.touschek import acceptance_for, lifetime_h, microbunch_enhancement
from ws_repro.track import Options, track

TAU = 3093                       # turns per longitudinal damping time
REAL = dict(backaction=0.641, ibs=True, radiator_qe=4.05e-7)

# name: (parameter changes, option changes, value quoted in the note)
ROWS = {
    "paper_assumptions":        ({}, {}, "645 W"),
    "paper_sigmaDz_0":          ({"sigma_Dz": 0.0}, {}, "3.2 kW"),
    "A1_3.6e-3":                ({"A1": 3.6e-3}, {}, "9.5 W"),
    "A1_3.8e-3":                ({"A1": 3.8e-3}, {}, "5 W"),
    "A1_4.1e-3":                ({"A1": 4.1e-3}, {}, "37 W"),
    "A1_4.3e-3":                ({"A1": 4.3e-3}, {}, "55 W"),
    "paper_plus_backaction":    ({}, {"backaction": 0.641}, "34 W"),
    "realistic_no_csr_1tau":    ({"sigma_Dz": 0.0}, {**REAL, "turns": TAU}, "about 40 W"),
    "realistic_with_csr_1tau":  ({"sigma_Dz": 0.0}, {**REAL, "csr": True, "turns": TAU}, "about 3 W or less"),
    "paper_sigmaDz_0.091nm":    ({"sigma_Dz": 0.091e-9}, {}, "(124001 ideal line)"),
    "paper_sigmaDz_0.18nm":     ({"sigma_Dz": 0.18e-9}, {}, "(124001 at sigma_delta 1.4e-3)"),
}


def run(item):
    name, (pc, oc, quoted), particles = item
    r = track(replace(Params(), **pc), replace(Options(particles=particles), **oc))
    return name, {**r, "quoted": quoted}


def formulas() -> dict:
    """Paper formulas only (no tracking): how the printed 1 kW arises and what one amplitude gives."""
    p = Params()
    kL = p.kL()
    A1_eq21 = 1 / (kL * p.r56)                                  # Eq. 21 (h1 = 1 / r56): sets the printed eps_z
    A1_eq12 = p.A1                                              # Eq. 12 optimum: sets the printed b76
    b1g = lambda eps: math.exp(-kL ** 2 * eps * p.beta_z0 / 2)  # Gaussian b1 at the FF entrance
    lo, hi = 1.5e-3, 6e-3
    for _ in range(80):                                         # sigma_delta,eff giving exactly 1 kW at b76 0.095
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if P_eq23_W(p, 0.095, m) > 1000 else (lo, m)
    e1 = eps_eq20(p, A1_eq12)
    b_one = S_bbhg(p, A1_eq12) * b1g(e1)
    return {
        "A1_eq21": A1_eq21, "A1_eq12": A1_eq12,
        "eps_z_eq20_at_A1_eq21_nm": eps_eq20(p, A1_eq21) * 1e9,
        "b76_eq12_printed_bunch": S_bbhg(p, A1_eq12) * b1g(0.34e-9),
        "sigma_delta_eff_for_1kW": 0.5 * (lo + hi),
        "P_printed_b76_ring_spread_W": P_eq23_W(p, 0.095, 1.4e-3),
        "one_amplitude": {"eps_z_nm": e1 * 1e9, "b76": b_one,
                          "P_W": P_eq23_W(p, b_one, sde_after_mod1(A1_eq12, math.sqrt(e1 / p.beta_z0)))},
    }


def touschek() -> dict:
    """Lifetimes of the scan rows above 1 kW (0.9 GeV, 1 nm, 3 % acceptance, uniform line density = generous),
    and the acceptance a 1 h lifetime needs at the Weak-Strong Table 1 point."""
    rows = {"6.2 kW (0.9 GeV, 1.2 kA peak, ring 71.7 m)": dict(E=0.9e9, I_peak=1200, C_ring=71.7),
            "17-19 kW (0.9 GeV, 2.4 kA peak, ring 56.0 m)": dict(E=0.9e9, I_peak=2400, C_ring=56.0)}
    out = {k: {"lifetime_h": lifetime_h(**v), "acceptance_for_1h": acceptance_for(1.0, **v)} for k, v in rows.items()}
    enh = microbunch_enhancement(248e-9, 1.03e-6)
    out["Table 1 point: acceptance for 1 h, ring beta 5 / 10 / 20 m"] = [
        acceptance_for(1.0, E=600e6, I_peak=50, C_ring=311, beta_low=0.94, L_low=3.9, beta_ring=b, enhancement=enh)
        for b in (5.0, 10.0, 20.0)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--particles", type=int, default=50_000)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--only", default="", help="comma-separated row names")
    a = ap.parse_args()
    sel = [k for k in ROWS if not a.only or k in a.only.split(",")]
    res = {"formulas": formulas(), "touschek": touschek(), "tracking": {}}
    with ProcessPoolExecutor(a.jobs) as ex:
        for name, r in ex.map(run, [(k, ROWS[k], a.particles) for k in sel]):
            res["tracking"][name] = r
    with open("results.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res["formulas"], indent=1))
    print(json.dumps(res["touschek"], indent=1))
    print(f"{'row':26s} {'P (Eq.23)':>10s} {'P via b1':>10s} {'b76':>7s} {'floor':>7s} {'eps_z nm':>8s} {'sz nm':>6s}  note")
    for k in sel:
        r = res["tracking"][k]
        print(f"{k:26s} {r['P_W']:10.1f} {r['P_W_from_b1']:10.1f} {r['bR']:7.4f} {r['noise_floor_bR']:7.4f} "
              f"{r['eps'] * 1e9:8.3f} {r['sz'] * 1e9:6.0f}  {r['quoted']}")


if __name__ == "__main__":
    main()
