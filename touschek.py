"""Touschek lifetime: Bruck formula, round beam, flat-top line density times a microbunch enhancement."""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad

from model import C_LIGHT, ME_EV, QE, R_E, bisect


def microbunch_enhancement(sigma_z: float, period: float, m_max: int = 20) -> float:
    """period <rho^2> / <rho>^2 for a train of Gaussian microbunches (rms sigma_z); 1 means uniform."""
    m = np.arange(-m_max, m_max + 1)
    return float(period / (2 * math.sqrt(math.pi) * sigma_z)
                 * np.sum(np.exp(-(m * period) ** 2 / (4 * sigma_z ** 2))))


def touschek_rate(E: float, eps: float, beta: float, line_density: float, acceptance: float) -> float:
    """Loss rate [1/s] at one location: round beam of size sqrt(eps beta), momentum acceptance as a fraction."""
    g = E / ME_EV
    sx, sxp = math.sqrt(eps * beta), math.sqrt(eps / beta)
    e = (acceptance / (g * sxp)) ** 2
    E1 = quad(lambda u: math.exp(-u) / u, e, np.inf)[0]
    E_log = quad(lambda u: math.log(u) * math.exp(-u) / u, e, np.inf)[0]
    D = math.sqrt(e) * (-1.5 * math.exp(-e) + e / 2 * E_log + 0.5 * (3 * e - e * math.log(e) + 2) * E1)
    return R_E ** 2 * C_LIGHT * math.sqrt(math.pi) * line_density / (4 * math.pi * sx * sx) / (g * g * acceptance ** 3) * D


def lifetime_h(E: float, I_peak: float, C_ring: float, acceptance: float = 0.03, eps: float = 1e-9,
               beta_low: float = 0.3, L_low: float = 2.3, beta_ring: float = 20.0, enhancement: float = 1.0) -> float:
    """Ring average: L_low metres at beta_low (radiator region), the rest at beta_ring.
    enhancement = 1 is a uniform (fully debunched) line density, the most generous case."""
    line_density = I_peak / (QE * C_LIGHT) * enhancement
    low = L_low * touschek_rate(E, eps, beta_low, line_density, acceptance)
    rest = (C_ring - L_low) * touschek_rate(E, eps, beta_ring, line_density, acceptance)
    return 1 / ((low + rest) / C_ring) / 3600


def acceptance_for(hours: float, **ring) -> float:
    """Momentum acceptance that gives a lifetime of `hours` (bisection in log acceptance)."""
    log_a = bisect(lambda u: lifetime_h(acceptance=math.exp(u), **ring) - hours, math.log(0.005), math.log(0.5), 60)
    return math.exp(log_a)
