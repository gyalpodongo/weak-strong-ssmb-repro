"""1D longitudinal multi-turn tracking of the Weak-Strong ring over one laser period (periodic window).

One turn, starting at the FF entrance (alpha_z = 0):
    CSR kick, +R56, Mod 1, +r56,
    radiator centre: readout, back-action, radiator quantum excitation, sigma_Dz path error,
    -r56, Rev Mod 1, -R56, half ring (-eta_C0/2),
    Mod 0 with damping, ring quantum excitation and IBS, ideal RF, half ring (-eta_C0/2).
With every Options effect off this is the paper's model: lasers, damping, ring quantum excitation and the lumped
sigma_Dz. The effects (ours, not the paper's) add what the paper leaves out or lumps; the README gives their status.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from model import C_LIGHT, QE, R_E, Params, bunching_gain, radiated_power_W, radiator

INITIAL_EMITTANCE = 0.34e-9                    # printed eps_z; the start of tracking only [m]
CSR_BENDS = ((0.5, 0.79 * 52.0),               # OURS (rho [m], length [m]): damping wigglers
             (2.0, 2 * math.pi * 2.0))         # OURS: arc dipoles
CSR_MAX_HARMONIC = 80
BEAM_SIZE = math.sqrt(1e-9 * 2.0)              # OURS: rms transverse size for CSR shielding (1 nm, beta 2 m)
RING_CIRCUMFERENCE = 311.0                     # OURS: for IBS [m]
COULOMB_LOG = 10.0                             # OURS: for IBS


@dataclass(frozen=True)
class Options:
    particles: int = 50_000
    turns: int = 10_000          # 3.2 damping times (one damping time is 1 / alpha_L = 3093 turns)
    seed: int = 1
    backaction: float = 0.0      # radiator coherent field, as a fraction of the Eq. 23 loss (OURS: 0.641, 3D code)
    csr: bool = False            # 1D steady-state CSR of the microbunch (outside 1D validity at these harmonics)
    ibs: bool = False            # IBS, Bane high-energy form, round 1 nm beam
    radiator_qe: float = 0.0     # radiator quantum excitation, rms delta per pass (OURS: 4.05e-7)
    record_every: int = 20       # turns between recorded moments


def wrap(z: np.ndarray, lambda_L: float) -> np.ndarray:
    """Fold z into one laser period centred on 0."""
    return z - lambda_L * np.floor(z / lambda_L + 0.5)


def centred(z: np.ndarray, kL: float, lambda_L: float) -> tuple[complex, np.ndarray]:
    """Bunching b1 and z measured from the microbunch centre."""
    b1 = np.mean(np.exp(1j * kL * z))
    return b1, wrap(z - np.angle(b1) / kL, lambda_L)


def moments(z: np.ndarray, delta: np.ndarray, kL: float, lambda_L: float) -> dict:
    """|b1|, rms bunch length and rms emittance at the FF entrance."""
    b1, zc = centred(z, kL, lambda_L)
    dm = delta - delta.mean()
    eps = math.sqrt(max(np.mean(zc * zc) * np.mean(dm * dm) - np.mean(zc * dm) ** 2, 0.0))
    return {"b1": abs(b1), "sigma_z": math.sqrt(np.mean(zc * zc)), "eps_z": eps}


def modulate(z: np.ndarray, delta: np.ndarray, kL: float, A1: float, A3: float, sign: int) -> np.ndarray:
    """Mod 1 (sign +1) or Rev Mod 1 (sign -1): fundamental plus third harmonic."""
    return delta + sign * (A1 * np.sin(kL * z) + A3 * np.sin(3 * kL * z))


def csr_impedance(k: float) -> complex:
    """1D steady-state CSR impedance summed over CSR_BENDS, with a transverse shielding factor."""
    total = 0j
    for rho, length in CSR_BENDS:
        shielding = min(1.0, ((1 / k) ** 2 * rho) ** (1 / 3) / BEAM_SIZE)
        re = math.sqrt(3) / 2 * 56.3 * k ** (1 / 3) * rho ** (-2 / 3) * length * shielding
        total += complex(re, -re / math.sqrt(3))
    return total


def csr_kick(z: np.ndarray, kL: float, I_peak: float, E0_eV: float) -> np.ndarray:
    """Energy change per particle from the harmonics whose |b_m| is above the 3 / sqrt(N) noise level."""
    noise = 3 / math.sqrt(z.size)
    phasor = np.exp(1j * kL * z)
    phasor_m, total = np.ones_like(phasor), np.zeros_like(phasor)
    for m in range(1, CSR_MAX_HARMONIC + 1):
        phasor_m = phasor_m * phasor
        b_m = complex(phasor_m.mean())
        if abs(b_m) > noise:
            source = csr_impedance(m * kL) * b_m.conjugate()
            total += 2 * I_peak * source * np.exp(-0.5j * math.pi) * phasor_m
    return total.imag / E0_eV


def ibs_variance(p: Params, sigma_z: float) -> float:
    """Per-turn increase of sigma_delta^2 from IBS (Bane high-energy form)."""
    particles_per_period = p.I_peak / (QE * C_LIGHT) * p.lambda_L / max(sigma_z, 1e-9)
    return (2 * (RING_CIRCUMFERENCE / C_LIGHT) * R_E ** 2 * C_LIGHT * particles_per_period * COULOMB_LOG
            * 2 ** -0.5 / (16 * p.gamma ** 3 * p.eps_perp ** 1.5))


def track(p: Params, o: Options) -> dict:
    """Track o.turns turns from a Gaussian beam at the printed eps_z; return equilibrium averages."""
    kL, lambda_L, R56, r56 = p.kL, p.lambda_L, p.R56, p.r56
    A1, A3, A0 = p.A1, p.A1 * p.A3_over_A1, p.A0
    half_ring = p.eta_C0 / 2
    damping = 1 - 2 * p.alpha_L
    ring_qe_rms = math.sqrt(p.beta_avg_over_M0 * p.quantum_excitation)
    peak_W_per_b2 = radiated_power_W(p, 1.0, 0.0) * p.I_peak / p.I_avg

    rng = np.random.default_rng(o.seed)
    N = o.particles
    z = rng.normal(0, math.sqrt(INITIAL_EMITTANCE * p.beta_z0), N)
    delta = rng.normal(0, math.sqrt(INITIAL_EMITTANCE / p.beta_z0), N)
    csr = np.zeros(N)
    history = []
    for turn in range(o.turns):
        record = turn % o.record_every == 0 or turn == o.turns - 1
        if record:
            row = {"turn": turn, **moments(z, delta, kL, lambda_L)}

        delta = delta + csr
        z = z + R56 * delta
        delta = modulate(z, delta, kL, A1, A3, +1)
        z = z + r56 * delta

        if record or o.backaction > 0:
            phasor_n = np.exp(-1j * p.n * kL * z)
            b_n = np.mean(phasor_n)
        if record:
            row.update(b_n=abs(b_n), sigma_delta_R=float(delta.std()))
            history.append(row)
        if o.backaction > 0:
            # coherent loss, locked to the bunching it radiates from
            amplitude = o.backaction * peak_W_per_b2 * abs(b_n) / p.I_peak / p.E0_eV
            delta = delta - amplitude * np.real(phasor_n * np.conj(b_n)) / max(abs(b_n), 1e-30)
        if o.radiator_qe > 0:
            delta = delta + rng.normal(0, o.radiator_qe, N)
        if p.sigma_Dz > 0:
            z = z + rng.normal(0, p.sigma_Dz, N)

        z = z - r56 * delta
        delta = modulate(z, delta, kL, A1, A3, -1)
        z = z - R56 * delta
        z = z - half_ring * delta

        delta = (delta + A0 * np.sin(kL * z)) * damping
        variance = ring_qe_rms * ring_qe_rms
        if o.ibs:
            _, zc = centred(z, kL, lambda_L)
            variance += ibs_variance(p, math.sqrt(np.mean(zc * zc)))
        delta = delta + rng.normal(0, math.sqrt(variance), N)
        delta = delta - delta.mean()                            # ideal RF restores the mean energy loss
        z = wrap(z - half_ring * delta, lambda_L)
        if o.csr:
            csr = csr_kick(z, kL, p.I_peak, p.E0_eV)
    return summarize(p, o, history)


def summarize(p: Params, o: Options, history: list[dict]) -> dict:
    """Average the recorded moments over the last damping time (or the last half of a shorter run)."""
    tau = 1 / p.alpha_L
    turns = np.array([row["turn"] for row in history])
    window = turns >= turns[-1] - min(tau, turns[-1] / 2)
    keys = ("b1", "sigma_z", "eps_z", "b_n", "sigma_delta_R")
    mean = {k: float(np.mean([row[k] for row, keep in zip(history, window) if keep])) for k in keys}
    b_n_from_b1 = bunching_gain(p, p.A1) * mean["b1"] * math.exp(-0.5 * (p.n * p.kL * p.sigma_Dz) ** 2)
    return {
        **mean,
        "noise_floor_b_n": 1 / math.sqrt(o.particles),
        "power_W": radiated_power_W(p, mean["b_n"], mean["sigma_delta_R"]),
        "b_n_from_b1": b_n_from_b1,
        "power_from_b1_W": radiated_power_W(p, b_n_from_b1, mean["sigma_delta_R"]),
        "turns": o.turns,
        "damping_times": o.turns / tau,
        "particles": o.particles,
        "seed": o.seed,
        "N_u": radiator(p).N_u,
    }
