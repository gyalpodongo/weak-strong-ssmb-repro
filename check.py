"""Quick self-check (about 10 s): reruns two rows at 2e4 particles and checks them against the note."""
from reproduce import ROWS, run, touschek

checks = {"paper_assumptions": (600, 700), "realistic_no_csr_1tau": (25, 45)}
for name, (lo, hi) in checks.items():
    _, r = run((name, ROWS[name], 20_000))
    assert lo <= r["P_W"] <= hi, f"{name}: {r['P_W']:.1f} W outside {lo} to {hi} W"
    print(f"ok  {name:24s} {r['P_W']:7.1f} W  (note: {r['quoted']})")
t = touschek()
for k, v in t.items():
    if isinstance(v, dict):
        assert v["lifetime_h"] < 0.2, k
        print(f"ok  Touschek {k}: {v['lifetime_h']:.3f} h")
print("all checks passed")
