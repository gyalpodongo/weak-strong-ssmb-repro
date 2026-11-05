"""Weak-Strong SSMB (Deng et al., arXiv:2607.09243v1): inputs and the paper's own formulas.

Every input is tagged PAPER (Table 1 or text, page given) or OURS (not printed; why we chose it).
Sign convention (paper Methods): z' = z + R56 delta, delta' = delta + A sin(kL z); the ring outside the
final focus (FF) has momentum compaction -etaC0.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.special import erf, jv

ME_EV, R_E, LAMBAR_C, C_LIGHT, QE = 510998.95, 2.8179403262e-15, 3.8615926796e-13, 299792458.0, 1.602176634e-19


@dataclass
class Params:
    E0_eV: float = 600e6            # PAPER Table 1
    I_avg: float = 0.5              # PAPER Table 1
    I_peak: float = 50.0            # PAPER Sec. 2.3 (1 % fill)
    etaC0: float = 1.0e-4           # PAPER Table 1 |eta C0|; sign OURS (0 < h0 eta C0 << 1)
    r56: float = 50.1e-6            # PAPER Table 1
    n: int = 76                     # PAPER: radiated harmonic; R56 = n r56 (Methods), Table 1 prints 3.8 mm
    lambda_L: float = 1030e-9       # PAPER Table 1
    beta_z0: float = 181e-6         # PAPER Methods (beta_z at the FF entrance, alpha_z = 0 by symmetry)
    A3_over_A1: float = -0.064      # PAPER Methods (third-harmonic modulation)
    A1: float = 3.973e-3            # OURS (not printed): Eq. 12 optimum at the printed 248 nm bunch
    sigma_Dz: float = 0.14e-9       # PAPER Methods (effective FF path error per pass)
    beta_avg_over_M0: float = 1.2   # PAPER Methods (<beta_z> = 1.2 beta_zM0 for quantum excitation)
    Js: float = 2.0                 # PAPER Methods
    P_W: float = 97e3               # PAPER Table 1 (wiggler power); U0 = P_W / I_avg = 194 keV
    wig_B: float = 4.0              # PAPER Table 1
    wig_L: float = 52.0             # PAPER Table 1
    eps_perp: float = 1e-9          # PAPER Table 1
    rad_B: float = 0.8              # PAPER Table 1
    rad_L: float = 1.88             # PAPER Table 1
    rad_beta_over_L: float = 0.5    # PAPER (effective beta = L_uR / 2)

    def gamma(self) -> float:
        return self.E0_eV / ME_EV

    def kL(self) -> float:
        return 2 * math.pi / self.lambda_L

    def R56(self) -> float:
        return self.n * self.r56

    def alpha_L(self) -> float:
        """Longitudinal damping per turn (amplitude): U0 Js / (2 E0)."""
        return (self.P_W / self.I_avg) * self.Js / (2 * self.E0_eV)

    def x0(self) -> float:
        """h0 etaC0 such that the symmetric ring map D(-eta/2) K0 D(-eta/2) gives beta_z0 at the FF entrance."""
        q = (self.etaC0 / self.beta_z0) ** 2
        lo, hi = 1e-6, 2.0
        for _ in range(200):
            x = 0.5 * (lo + hi)
            lo, hi = (lo, x) if (x - x * x / 4) - q * (1 - x / 4) ** 2 > 0 else (x, hi)
        return 0.5 * (lo + hi)

    def A0(self) -> float:
        return self.x0() / (self.etaC0 * self.kL())

    def beta_zM0(self) -> float:
        return self.beta_z0 / (1 - self.x0() / 4)

    def D_qe(self) -> float:
        """Per-turn sigma_delta^2 increment from the damping wigglers (standard quantum-excitation coefficient)."""
        g, rho = self.gamma(), self.E0_eV / C_LIGHT / self.wig_B
        return 55 / (24 * math.sqrt(3)) * R_E * LAMBAR_C * g ** 5 * 4 / (3 * math.pi) * self.wig_L / rho ** 3


def radiator(p: Params) -> dict:
    """Planar radiator resonant at lambda_L / n with K = 0.9336 B[T] lambda_u[cm]; transverse factor of Eq. 23."""
    g, lam_R = p.gamma(), p.lambda_L / p.n
    lo, hi = 0.005, 0.05
    for _ in range(80):
        lu = 0.5 * (lo + hi)
        K = 0.9336 * p.rad_B * lu * 100
        lo, hi = (lu, hi) if lu / (2 * g * g) * (1 + K * K / 2) < lam_R else (lo, lu)
    Nu, chi = p.rad_L / lu, K * K / (4 + 2 * K * K)
    S = 2 * math.pi * p.eps_perp * p.rad_beta_over_L / lam_R
    FF = 2 / math.pi * (math.atan(1 / (2 * S)) + S * math.log((2 * S) ** 2 / ((2 * S) ** 2 + 1)))
    return {"lambda_u": lu, "K": K, "N_u": Nu, "chi": chi, "JJ": float(jv(0, chi) - jv(1, chi)), "FF_perp": FF}


def D_delta(Nu: float, sde: float) -> float:
    xi = Nu * sde
    return 1.0 if xi == 0 else float(erf(math.sqrt(2) * math.pi * xi) ** 2 / (8 * math.pi * xi * xi))


def P_eq23_W(p: Params, b: float, sde: float) -> float:
    """Paper Eq. 23 (peak power, kW) times the duty I_avg / I_peak, in W."""
    r = radiator(p)
    pk_kW = 1.183 * r["N_u"] * r["chi"] * r["JJ"] ** 2 * r["FF_perp"] * b * b * D_delta(r["N_u"], sde) * p.I_peak ** 2
    return pk_kW * 1e3 * p.I_avg / p.I_peak


def S_bbhg(p: Params, A1: float) -> float:
    """Eq. 12 with R56 = n r56: only the m = -1 term survives, so b_n = S(A1) x |b1 of the FF-entrance beam|."""
    th = -p.n * p.kL() * p.r56
    return abs(sum(jv(p.n + 1 - 3 * k, th * A1) * jv(k, th * p.A3_over_A1 * A1) for k in range(-30, 31)))


def sde_after_mod1(A1: float, sd_ring: float, a3: float = -0.064) -> float:
    return math.hypot(sd_ring, A1 * math.sqrt((1 + a3 ** 2) / 2))


def eps_eq20(p: Params, A1: float) -> float:
    """Paper Eq. 20: (ring quantum excitation + FF path error) / (2 alpha_L), in m."""
    two_a = 2 * p.alpha_L()
    ring = p.D_qe() * p.beta_avg_over_M0 * p.beta_zM0() / 2 / two_a
    beta_M1 = p.beta_z0 + p.R56() ** 2 / p.beta_z0
    ff = (A1 * p.kL() * p.sigma_Dz) ** 2 * (1 + 9 * p.A3_over_A1 ** 2) / 2 * beta_M1 / 2 / two_a
    return ring + ff
