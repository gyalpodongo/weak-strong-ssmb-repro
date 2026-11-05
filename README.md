# Weak-Strong SSMB power estimate

Code for our note checking the radiated power of the Weak-Strong SSMB ring (X. Deng et al., arXiv:2607.09243).

The paper quotes 1 kW at 13.5 nm. Tracking the paper's own longitudinal model with its stated assumptions, we
get about 645 W (runs vary by a few percent with the random seed; `results.json` has 647 W). When we model
radiator back-action, IBS and CSR explicitly instead of through the 0.14 nm path-length budget, we get roughly
3 to 40 W.

```
pip install -r requirements.txt
python3 check.py        # 15 s sanity check
python3 reproduce.py    # all cases, about 1 min, writes results.json
```

## Results

| case | power |
|---|---|
| paper's assumptions | 647 W |
| same, path error 0 / 0.091 / 0.18 nm | 3.2 / 1.25 / 0.41 kW |
| Mod 1 amplitude -4 % / +3 % | 6 / 37 W |
| plus radiator back-action | 36 W |
| back-action, IBS, radiator quantum excitation | 36 W |
| all of the above plus 1D CSR | below 3 W |

The scan points above 1 kW (0.9 GeV, 1.2 to 2.4 kA peak) have Touschek lifetimes of 0.04 to 0.09 h.

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
- back-action scaled by 0.641 from a 3D radiator calculation; IBS in Bane's high-energy form for a 311 m ring

The CSR case uses the 1D steady-state impedance, which is outside its validity at these harmonics, so treat
the low end of the range as soft. No Weak-Strong lattice is published, and this model doesn't need one.

MIT license.
