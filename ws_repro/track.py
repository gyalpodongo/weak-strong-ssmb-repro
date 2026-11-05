"""1D longitudinal multi-turn tracker of the Weak-Strong ring on one laser period (periodic window).

One turn, starting at the FF entrance (alpha_z = 0):
  [CSR kick] +R56 | Mod 1 | +r56 | radiator centre (readout; sigma_Dz noise; back-action; radiator QE) | -r56 |
  Rev Mod 1 | -R56 | ring -etaC0/2 | Mod 0 (sinusoidal) + damping + ring QE [+ IBS] | ideal RF | ring -etaC0/2.
Paper assumptions = lasers, damping, ring QE and the lumped sigma_Dz only. The "realistic" options add effects the
paper leaves out or lumps (radiator back-action, 1D steady-state CSR, IBS, radiator QE); see README for their status.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .model import C_LIGHT, QE, R_E, Params, P_eq23_W, S_bbhg, radiator


@dataclass
class Options:
    particles: int = 50_000
    turns: int = 10_000              # 3.2 damping times (1 tau = 1 / alpha_L = 3093 turns)
    seed: int = 1
    backaction: float = 0.0          # radiator coherent field, x Eq. 23 self-consistent loss (OURS 0.641: 3D calibration)
    csr: bool = False                # 1D steady-state CSR of the microbunch: wigglers + arc dipoles (OUTSIDE 1D validity)
    ibs: bool = False                # Bane high-energy IBS, round 1 nm beam, ring C0 311 m (OURS)
    radiator_qe: float = 0.0         # radiator quantum excitation, rms delta per pass (OURS 4.05e-7)
    every: int = 20


CSR_SETS = ((0.5, 0.79 * 52.0), (2.0, 2 * math.pi * 2.0))   # (rho [m], effective length [m]): wigglers, arc dipoles (OURS)
SIG_X, C0_RING, MMAX = math.sqrt(1e-9 * 2.0), 311.0, 80


def wrap(z, lam):
    return z - lam * np.floor(z / lam + 0.5)


def _csr_Z(k: float) -> complex:
    """Steady-state CSR impedance per set (Re = (sqrt3/2) 56.3 k^(1/3) rho^(-2/3) Ohm/m x L x transverse factor)."""
    tot = 0j
    for rho, L in CSR_SETS:
        re = math.sqrt(3) / 2 * 56.3 * k ** (1 / 3) * rho ** (-2 / 3) * L * min(1.0, ((1 / k) ** 2 * rho) ** (1 / 3) / SIG_X)
        tot += complex(re, -re / math.sqrt(3))
    return tot


def _csr_kick(z, kL: float, I_P: float, E0: float) -> np.ndarray:
    """delta change sum_m Re/Im of 2 I_P Z_m b_m* e^{i m kL z}, harmonics with |b_m| above 3/sqrt(N)."""
    noise = 3 / math.sqrt(z.size)
    e1 = np.exp(1j * kL * z)
    em, acc = np.ones_like(e1), np.zeros_like(e1)
    for m in range(1, MMAX + 1):
        em = em * e1
        b = complex(em.mean())
        if abs(b) > noise:
            S = _csr_Z(m * kL) * b.conjugate()
            acc += 2 * I_P * S * np.exp(-0.5j * math.pi) * em
    return acc.imag / E0


def _ibs_var(p: Params, sz: float) -> float:
    g = p.gamma()
    Ns = p.I_peak / (QE * C_LIGHT) * p.lambda_L / max(sz, 1e-9)
    return 2 * (C0_RING / C_LIGHT) * R_E ** 2 * C_LIGHT * Ns * 10.0 * 2 ** -0.5 / (16 * g ** 3 * p.eps_perp ** 1.5)


def track(p: Params, o: Options) -> dict:
    kL, lam, n = p.kL(), p.lambda_L, p.n
    R56, r56, A1, A3, A0 = p.R56(), p.r56, p.A1, p.A1 * p.A3_over_A1, p.A0()
    damp, sq = 1 - 2 * p.alpha_L(), math.sqrt(p.beta_avg_over_M0 * p.D_qe())
    eta2 = p.etaC0 / 2
    pk_per_b2 = P_eq23_W(p, 1.0, 0.0) * p.I_peak / p.I_avg            # peak W per |b|^2 (D_delta = 1)

    def mod(z, d, s):
        return d + s * (A1 * np.sin(kL * z) + A3 * np.sin(3 * kL * z))

    rng = np.random.default_rng(o.seed)
    z = rng.normal(0, math.sqrt(0.34e-9 * p.beta_z0), o.particles)
    d = rng.normal(0, math.sqrt(0.34e-9 / p.beta_z0), o.particles)
    N, hist, csr = z.size, [], np.zeros(o.particles)
    for t in range(o.turns):
        rec = t % o.every == 0 or t == o.turns - 1
        if rec:
            c = np.mean(np.exp(1j * kL * z))
            zc = wrap(z - np.angle(c) / kL, lam)
            dm = d - d.mean()
            h = {"turn": t, "b1": abs(c), "sz": math.sqrt(np.mean(zc * zc)),
                 "eps": math.sqrt(max(np.mean(zc * zc) * np.mean(dm * dm) - np.mean(zc * dm) ** 2, 0.0))}
        d = d + csr
        z = z + R56 * d
        d = mod(z, d, 1)
        z = z + r56 * d
        if rec or o.backaction > 0:
            B = np.mean(np.exp(-1j * n * kL * z))
        if rec:
            h.update(bR=abs(B), sdR=float(d.std()))
            hist.append(h)
        if o.backaction > 0:                                 # coherent loss locked to the bunching it radiates from
            amp = o.backaction * pk_per_b2 * abs(B) / p.I_peak / p.E0_eV
            d = d - amp * np.real(np.exp(-1j * n * kL * z) * np.conj(B)) / max(abs(B), 1e-30)
        if o.radiator_qe > 0:
            d = d + rng.normal(0, o.radiator_qe, N)
        if p.sigma_Dz > 0:
            z = z + rng.normal(0, p.sigma_Dz, N)
        z = z - r56 * d
        d = mod(z, d, -1)
        z = z - R56 * d
        z = z - eta2 * d
        d = (d + A0 * np.sin(kL * z)) * damp
        var = sq * sq
        if o.ibs:
            c = np.mean(np.exp(1j * kL * z))
            zc = wrap(z - np.angle(c) / kL, lam)
            var += _ibs_var(p, math.sqrt(np.mean(zc * zc)))
        d = d + rng.normal(0, math.sqrt(var), N)
        d = d - d.mean()                                     # ideal RF restores the mean energy loss
        z = wrap(z - eta2 * d, lam)
        csr = _csr_kick(z, kL, p.I_peak, p.E0_eV) if o.csr else csr
    return summarize(p, o, hist)


def summarize(p: Params, o: Options, hist: list[dict]) -> dict:
    """Average over the last damping time (or the last half of the run if shorter)."""
    tau = 1 / p.alpha_L()
    T = np.array([h["turn"] for h in hist])
    w = T >= T[-1] - min(tau, T[-1] / 2)
    m = {k: float(np.mean([h[k] for h, s in zip(hist, w) if s])) for k in ("b1", "sz", "eps", "bR", "sdR")}
    floor = 1 / math.sqrt(o.particles)
    b_est = S_bbhg(p, p.A1) * m["b1"] * math.exp(-0.5 * (p.n * p.kL() * p.sigma_Dz) ** 2)
    return {**m, "noise_floor_bR": floor, "P_W": P_eq23_W(p, m["bR"], m["sdR"]),
            "bR_from_b1": b_est, "P_W_from_b1": P_eq23_W(p, b_est, m["sdR"]),
            "turns": o.turns, "turns_over_tau": o.turns / tau, "particles": o.particles, "seed": o.seed,
            "N_u": radiator(p)["N_u"]}
