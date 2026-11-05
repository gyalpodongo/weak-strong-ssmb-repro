"""Quick self-check (about 15 s): two rows at 2e4 particles and the Touschek lifetimes, against the note."""
from reproduce import run_row, touschek

EXPECTED_W = {"paper_assumptions": (600, 700), "realistic_no_csr_1tau": (25, 45)}

for name, (low, high) in EXPECTED_W.items():
    r = run_row(name, particles=20_000)
    assert low <= r["power_W"] <= high, f"{name}: {r['power_W']:.1f} W outside {low} to {high} W"
    print(f"ok  {name:24s} {r['power_W']:7.1f} W  (note: {r['note']})")

for name, row in touschek().items():
    if isinstance(row, dict):
        assert row["lifetime_h"] < 0.2, name
        print(f"ok  Touschek {name}: {row['lifetime_h']:.3f} h")
print("all checks passed")
