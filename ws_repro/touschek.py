"""Touschek lifetime (Bruck formula, round beam, flat-top line density x microbunch enhancement)."""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad

from .model import C_LIGHT, ME_EV, QE, R_E


def microbunch_enhancement(sigma_z: float, lam: float, m_max: int = 20) -> float:
    """lambda <rho^2> / <rho>^2 over one period of a train of Gaussians (rms sigma_z, period lam); 1 = uniform."""
    m = np.arange(-m_max, m_max + 1)
    return float(lam / (2 * math.sqrt(math.pi) * sigma_z) * np.sum(np.exp(-(m * lam) ** 2 / (4 * sigma_z ** 2))))


def rate(E: float, eps: float, beta: float, line_density: float, acc: float) -> float:
    """Touschek loss rate [1/s] at one location: round beam sqrt(eps beta), momentum acceptance acc."""
    g = E / ME_EV
    sx, sxp = math.sqrt(eps * beta), math.sqrt(eps / beta)
    e = (acc / (g * sxp)) ** 2
    E1 = quad(lambda u: math.exp(-u) / u, e, np.inf)[0]
    El = quad(lambda u: math.log(u) * math.exp(-u) / u, e, np.inf)[0]
    D = math.sqrt(e) * (-1.5 * math.exp(-e) + e / 2 * El + 0.5 * (3 * e - e * math.log(e) + 2) * E1)
    return R_E ** 2 * C_LIGHT * math.sqrt(math.pi) * line_density / (4 * math.pi * sx * sx) / (g * g * acc ** 3) * D


def lifetime_h(E: float, I_peak: float, C_ring: float, acc: float = 0.03, eps: float = 1e-9, beta_low: float = 0.3,
               L_low: float = 2.3, beta_ring: float = 20.0, enhancement: float = 1.0) -> float:
    """Ring average: L_low metres at beta_low (radiator region), the rest at beta_ring. enhancement = 1 is a uniform
    (fully debunched) line density, the most generous case."""
    lam = I_peak / (QE * C_LIGHT) * enhancement
    r = (L_low * rate(E, eps, beta_low, lam, acc) + (C_ring - L_low) * rate(E, eps, beta_ring, lam, acc)) / C_ring
    return 1 / r / 3600


def acceptance_for(hours: float, **kw) -> float:
    lo, hi = 0.005, 0.5
    for _ in range(60):
        a = math.sqrt(lo * hi)
        lo, hi = (a, hi) if lifetime_h(acc=a, **kw) < hours else (lo, a)
    return math.sqrt(lo * hi)
