# Weak-Strong SSMB power estimate

Code for our note checking the radiated power of the Weak-Strong SSMB ring (X. Deng et al., arXiv:2607.09243).

The paper quotes 1 kW at 13.5 nm. Tracking the paper's own longitudinal model with its stated assumptions, we
get about 645 W (runs vary by a few percent with the random seed; `results.json` has 647 W). When we model
radiator back-action, IBS and CSR explicitly instead of through the 0.14 nm path-length budget, we get about
3 to 70 W. That range is model-dependent: we have not been able to validate those models against measurement (our checks
against MLS measurements and KARA instability thresholds did not pass).

```
pip install -r requirements.txt
python3 check.py        # 15 s sanity check
python3 reproduce.py    # all cases, about 1 min, writes results.json
```

## Notation

- Mod 0: the weak laser modulator that holds the microbunches. Mod 1: the strong one before the final focus.
  Rev Mod 1: cancels Mod 1 after the radiator.
- A1, A3: Mod 1 energy-modulation amplitudes at the laser wavenumber kL and its third harmonic.
- sigma_Dz: rms path error per pass from Mod 1 to Rev Mod 1 (0.14 nm in the paper).
- eps_z: longitudinal emittance. b76: bunching at the 76th laser harmonic (13.5 nm) at the radiator.

## Results

| case | power |
|---|---|
| paper's assumptions | 647 W |
| same, sigma_Dz 0 / 0.091 / 0.18 nm | 3.2 / 1.25 / 0.41 kW |
| same, path-error heating matched to the Eq. 21 budget | 0.91 kW |
| Mod 1 amplitude -4 % / +3 % (-8 % / +6 % in laser power) | 6 / 37 W |
| plus radiator back-action | 36 W |
| back-action, IBS, radiator quantum excitation | 36 W |
| same, without IBS | 38 W |
| same, back-action normalised so the beam's coherent loss equals the Eq. 23 power | 70 W |
| all of the above plus 1D CSR | at most 3 W at 1 damping time |

With CSR the bunching sits at the noise floor. A longer run (3.2 damping times, not in this repo) gives at most 0.3 W.

## Why 645 W and not 1 kW

- A path error leaves an energy error after Rev Mod 1. The beam spans about five laser periods at Mod 1, so the
  error averages over the laser phase: kL^2 (A1^2 + 9 A3^2) sigma_Dz^2 / 2 per pass. Tracking matches this to 0.3 %.
- That is 1.53 times the (sigma_Dz / r56)^2 / 2 per pass that reproduces the printed U0 from Eq. 21. The text on
  p. 19 writes the growth as (h1 sigma_Dz)^2, twice that.
- With the heating matched to Eq. 21 we get 0.91 kW. The rest is the energy spread in Eq. 23.
- The tracked eps_z (0.25 nm) is below the printed 0.34 nm: the Mod 0 bucket edge caps the rms growth, and the
  excess becomes an unbunched background. So b76 drops through the bunch shape, not the rms bunch length.
- A1 = 3.97e-3 with its third harmonic meets the linear compression condition |h1 r56| = 1 (0.98).

The lattice-free scan points above 1 kW (0.9 GeV, 1.2 to 2.4 kA peak) were parameter crossings, not designs, and have
Touschek lifetimes of 0.04 to 0.09 h. For a 1 h lifetime at the paper's point, the Bruck form in `touschek.py`
asks for 1.5 to 2.0 % momentum acceptance; our fuller estimate is about 2 to 4 %.

## Files

- `model.py`: parameters (each marked as from the paper or our choice) and the paper's Eqs. 12, 20, 23
- `tracking.py`: the one-turn map and the optional effects
- `touschek.py`: Touschek lifetime (Bruck)
- `reproduce.py`: every case, plus checks of the printed numbers
- `check.py`: quick check against the note

## Assumptions

Where the paper is silent we chose:

- A1 = 3.973e-3, the Eq. 12 optimum at the printed 248 nm bunch (it reproduces the printed b76)
- Mod 0 opposite the final focus, with symmetric ring halves
- Eq. 23 evaluated with the tracked bunching and the tracked energy spread at the radiator
- one laser period, longitudinal only, ideal RF
- back-action scaled by 0.641 from a 3D radiator calculation (or normalised to the Eq. 23 power); IBS in Bane's high-energy form for a 311 m ring

The CSR case uses the 1D steady-state impedance, which is outside its validity at these harmonics, so it overstates
CSR. We could not find a public Weak-Strong lattice, and this 1D model does not need one.

MIT license.
