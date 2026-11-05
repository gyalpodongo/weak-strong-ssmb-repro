"""Weak-Strong SSMB (X. Deng et al., arXiv:2607.09243v1): inputs and the paper's formulas.

Each input is marked PAPER (printed in Table 1 or the text) or OURS (not printed; reason given).
Sign convention (paper Methods): z' = z + R56 delta and delta' = delta + A sin(kL z).
The ring outside the final focus (FF) has momentum compaction -eta_C0.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.special import erf, jv

ME_EV = 510998.95                 # electron rest energy [eV]
R_E = 2.8179403262e-15            # classical electron radius [m]
LAMBDA_BAR_C = 3.8615926796e-13   # reduced Compton wavelength [m]
C_LIGHT = 299792458.0             # [m/s]
QE = 1.602176634e-19              # elementary charge [C]


def bisect(f, lo: float, hi: float, iterations: int = 100) -> float:
    """Root of an increasing function f on [lo, hi]."""
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class Params:
    E0_eV: float = 600e6              # PAPER Table 1: beam energy
    I_avg: float = 0.5                # PAPER Table 1: average current [A]
    I_peak: float = 50.0              # PAPER Sec. 2.3: peak current (1 % fill) [A]
    eta_C0: float = 1.0e-4            # PAPER Table 1 |eta C0|; sign OURS (0 < h0 eta_C0 << 1)
    r56: float = 50.1e-6              # PAPER Table 1 [m]
    n: int = 76                       # PAPER: radiated harmonic; R56 = n r56 (Methods; Table 1 prints 3.8 mm)
    lambda_L: float = 1030e-9         # PAPER Table 1: laser wavelength [m]
    beta_z0: float = 181e-6           # PAPER Methods: beta_z at the FF entrance, alpha_z = 0 by symmetry [m]
    A3_over_A1: float = -0.064        # PAPER Methods: third-harmonic share of the Mod 1 modulation
    A1: float = 3.973e-3              # OURS (not printed): Eq. 12 optimum at the printed 248 nm bunch
    sigma_Dz: float = 0.14e-9         # PAPER Methods: effective FF path error per pass [m]
    beta_avg_over_M0: float = 1.2     # PAPER Methods: <beta_z> = 1.2 beta_zM0 for quantum excitation
    J_s: float = 2.0                  # PAPER Methods: longitudinal damping partition
    wiggler_power_W: float = 97e3     # PAPER Table 1: damping wiggler power; U0 = P / I_avg = 194 keV
    wiggler_B: float = 4.0            # PAPER Table 1 [T]
    wiggler_L: float = 52.0           # PAPER Table 1 [m]
    eps_perp: float = 1e-9            # PAPER Table 1: transverse emittance [m]
    radiator_B: float = 0.8           # PAPER Table 1 [T]
    radiator_L: float = 1.88          # PAPER Table 1 [m]
    radiator_beta_over_L: float = 0.5 # PAPER: effective beta = L_uR / 2

    @property
    def gamma(self) -> float:
        return self.E0_eV / ME_EV

    @property
    def kL(self) -> float:
        return 2 * math.pi / self.lambda_L

    @property
    def R56(self) -> float:
        return self.n * self.r56

    @property
    def alpha_L(self) -> float:
        """Longitudinal amplitude damping per turn, U0 J_s / (2 E0)."""
        return (self.wiggler_power_W / self.I_avg) * self.J_s / (2 * self.E0_eV)

    @property
    def h0_eta_C0(self) -> float:
        """Mod 0 strength h0 eta_C0 for which the symmetric ring map D(-eta_C0/2) K0 D(-eta_C0/2) gives beta_z0."""
        q = (self.eta_C0 / self.beta_z0) ** 2
        return bisect(lambda x: (x - x * x / 4) - q * (1 - x / 4) ** 2, 1e-6, 2.0)

    @property
    def A0(self) -> float:
        """Mod 0 energy modulation amplitude."""
        return self.h0_eta_C0 / (self.eta_C0 * self.kL)

    @property
    def beta_zM0(self) -> float:
        """beta_z at Mod 0."""
        return self.beta_z0 / (1 - self.h0_eta_C0 / 4)

    @property
    def quantum_excitation(self) -> float:
        """Per-turn increase of sigma_delta^2 from the damping wigglers (standard coefficient)."""
        g, rho = self.gamma, self.E0_eV / C_LIGHT / self.wiggler_B
        return 55 / (24 * math.sqrt(3)) * R_E * LAMBDA_BAR_C * g ** 5 * 4 / (3 * math.pi) * self.wiggler_L / rho ** 3


@dataclass(frozen=True)
class Radiator:
    lambda_u: float   # period [m]
    K: float          # undulator parameter
    N_u: float        # number of periods
    chi: float        # K^2 / (4 + 2 K^2)
    JJ: float         # J0(chi) - J1(chi)
    F_perp: float     # transverse factor of Eq. 23


def radiator(p: Params) -> Radiator:
    """Planar radiator resonant at lambda_L / n, with K = 0.9336 B[T] lambda_u[cm]."""
    gamma, lambda_R = p.gamma, p.lambda_L / p.n

    def undulator_K(lambda_u: float) -> float:
        return 0.9336 * p.radiator_B * lambda_u * 100

    def detuning(lambda_u: float) -> float:
        K = undulator_K(lambda_u)
        return lambda_u / (2 * gamma * gamma) * (1 + K * K / 2) - lambda_R

    lambda_u = bisect(detuning, 0.005, 0.05)
    K = undulator_K(lambda_u)
    chi = K * K / (4 + 2 * K * K)
    S = 2 * math.pi * p.eps_perp * p.radiator_beta_over_L / lambda_R
    F_perp = 2 / math.pi * (math.atan(1 / (2 * S)) + S * math.log((2 * S) ** 2 / ((2 * S) ** 2 + 1)))
    return Radiator(lambda_u, K, p.radiator_L / lambda_u, chi, float(jv(0, chi) - jv(1, chi)), F_perp)


def D_delta(N_u: float, sigma_delta: float) -> float:
    """Eq. 23 energy-spread factor (1 for a cold beam)."""
    xi = N_u * sigma_delta
    if xi == 0:
        return 1.0
    return float(erf(math.sqrt(2) * math.pi * xi) ** 2 / (8 * math.pi * xi * xi))


def radiated_power_W(p: Params, b_n: float, sigma_delta: float) -> float:
    """Eq. 23 (peak power, kW) at bunching b_n and energy spread sigma_delta, times the duty I_avg / I_peak, in W."""
    r = radiator(p)
    peak_kW = (1.183 * r.N_u * r.chi * r.JJ ** 2 * r.F_perp * b_n * b_n
               * D_delta(r.N_u, sigma_delta) * p.I_peak ** 2)
    return peak_kW * 1e3 * p.I_avg / p.I_peak


def bunching_gain(p: Params, A1: float) -> float:
    """Eq. 12 at R56 = n r56: only the m = -1 term survives, so b_n = bunching_gain x |b1| at the FF entrance."""
    theta = -p.n * p.kL * p.r56
    return abs(sum(jv(p.n + 1 - 3 * k, theta * A1) * jv(k, theta * p.A3_over_A1 * A1) for k in range(-30, 31)))


def energy_spread_after_mod1(A1: float, sigma_delta_ring: float, A3_over_A1: float) -> float:
    """Ring energy spread plus the rms of the Mod 1 modulation."""
    return math.hypot(sigma_delta_ring, A1 * math.sqrt((1 + A3_over_A1 ** 2) / 2))


def emittance_eq20(p: Params, A1: float) -> float:
    """Eq. 20: equilibrium eps_z [m] = (ring quantum excitation + FF path error) / (2 alpha_L)."""
    two_alpha = 2 * p.alpha_L
    ring = p.quantum_excitation * p.beta_avg_over_M0 * p.beta_zM0 / 2 / two_alpha
    beta_M1 = p.beta_z0 + p.R56 ** 2 / p.beta_z0
    path_error = (A1 * p.kL * p.sigma_Dz) ** 2 * (1 + 9 * p.A3_over_A1 ** 2) / 2 * beta_M1 / 2 / two_alpha
    return ring + path_error
