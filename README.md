# Weak-Strong SSMB: reproduction of the power estimate

The Weak-Strong SSMB paper (X. Deng et al., arXiv:2607.09243v1) prints 1 kW of radiated power. Running the
paper's own 1D model with its stated assumptions gives about 645 W. Tracking is statistical: one run lands
anywhere from 640 to 660 W depending on the random seed, and the committed `results.json` (seed 1, 5e4
particles) shows 647 W. So 645 W (the note) and 647 W (`results.json`) are the same result, not two results.
Adding effects the paper leaves out or lumps (radiator back-action, IBS, radiator quantum excitation, CSR)
brings it to about 3 to 40 W.

This is the code behind the note "Weak-Strong SSMB (arXiv:2607.09243): an independent check of the power
estimate" (G. Dongo, S. Alberdi). It tracks the paper's longitudinal model to equilibrium and evaluates the
paper's Eq. 23 on the tracked beam. Dependencies: numpy and scipy.

## Quick start

```
pip install -r requirements.txt
python3 check.py        # about 15 s: two rows and the Touschek lifetimes against the note
python3 reproduce.py    # about 1 min on 8 cores: every row, writes results.json
python3 reproduce.py --only paper_assumptions --particles 200000
```

## Headline numbers

| case | results.json (seed 1, 5e4 particles) | note |
|---|---|---|
| paper's own assumptions | 647 W | about 645 W (paper prints 1 kW) |
| realistic, no CSR | 36 W | 36 W |
| realistic, with 1D CSR | at the noise floor | 3 W or less |
| Touschek lifetime, 6.2 kW scan point | 0.094 h | 0.09 h |
| Touschek lifetime, 17 to 19 kW scan points | 0.037 h | 0.04 h |

## Rows

All rows are in `ROWS` in `reproduce.py`; `reproduce.py` prints each next to the value quoted in the note.

| row | model | note |
|---|---|---|
| `paper_assumptions` | the paper's model: thin Mod 1 and Rev Mod 1 with the third harmonic, sinusoidal Mod 0, damping, ring quantum excitation at <beta_z> = 1.2 beta_zM0, Gaussian FF path error sigma_Dz = 0.14 nm per pass; 3.2 damping times | 645 W |
| `paper_sigmaDz_0`, `_0.091nm`, `_0.18nm` | same, other sigma_Dz | 3.2, 1.25, 0.41 kW |
| `A1_*` | same, Mod 1 amplitude A1 moved by -9 to +8 % | 6 W at -4 %, 37 W at +3 % |
| `paper_plus_backaction` | paper model plus the radiator's coherent field acting back on the bunched beam | 35 W |
| `realistic_no_csr_1tau` | no lumped sigma_Dz; back-action, IBS and radiator quantum excitation explicit; 1 damping time | 36 W |
| `realistic_with_csr_1tau` | as above plus 1D steady-state CSR of the microbunch in the wigglers and arc dipoles | 3 W or less |

`formulas()` checks the printed numbers without tracking: Eq. 23 needs sigma_delta,eff = 2.78e-3 for exactly
1 kW at b76 = 0.095; Eq. 20 gives eps_z = 0.35 nm at the Eq. 21 amplitude A1 = 3.27e-3, while b76 = 0.096
needs A1 = 3.97e-3; with one amplitude the formulas give 0.47 nm, b76 = 0.061 and 333 W.
`touschek()` gives the Touschek lifetimes of the parameter-scan points above 1 kW (at a 3 % acceptance) and
the acceptance a 1 h lifetime needs at the Table 1 point (1.5 to 2.0 %).

## Code

| file | contents |
|---|---|
| `model.py` | `Params` (every input, marked PAPER or OURS) and the paper's formulas: Eq. 12, 20, 23 |
| `tracking.py` | the one-turn map and `track()`; `Options` switches the effects beyond the paper on |
| `touschek.py` | Touschek lifetime (Bruck formula) |
| `reproduce.py` | the rows, the formula and Touschek checks; writes `results.json` |
| `check.py` | quick self-check |

## Assumptions (ours where the paper is silent)

- Inputs from Table 1 and the text; R56 = n r56 exactly; A3/A1 = -0.064.
- A1 = 3.973e-3 (not printed): the Eq. 12 optimum at the printed 248 nm bunch, which reproduces the printed b76.
- Mod 0 sits opposite the final focus (FF) with symmetric ring halves of -eta_C0/2, so alpha_z = 0 at the FF
  entrance; h0 eta_C0 = 0.28 from beta_z0 = 181 um. U0 = P_W / I_A = 194 keV, J_s = 2.
- Eq. 23 is evaluated with the tracked b76 and the tracked rms energy spread at the radiator (after Mod 1);
  beta_perp = L_uR / 2, eps_perp = 1 nm. "P via b1" uses b76 = S(A1) b1 (Eq. 12 at R56 = n r56), which has
  a lower noise floor than the direct b76.
- One laser period, periodic; no transverse dynamics, no RF slope, ideal energy restoration each turn.
- Realistic rows (ours, not the paper's): back-action amplitude self-consistent with Eq. 23 times 0.641 (a 3D
  radiator code calibration); IBS in the Bane high-energy form with a 311 m ring; radiator quantum excitation
  4.05e-7 per pass; CSR with the 1D steady-state impedance, rho 0.5 m over 0.79 x 52 m of wiggler and rho 2 m
  over 2 pi of arc.
- Touschek: Bruck formula, round 1 nm beam, uniform line density (the most generous case).

## Limits

- At the laser harmonics the microbunch is outside the validity of 1D steady-state CSR, so the CSR row
  overstates CSR; its power is a soft lower end, not a prediction.
- No lattice: no Weak-Strong ring lattice is published (Pan et al., PRAB 29, 103401 (2026) is a different
  400 MeV single-laser ring), and this 1D model does not need one.
- Tracked numbers move with the seed by a few percent at 5e4 particles; `--particles 200000` tightens them.

License: MIT.
