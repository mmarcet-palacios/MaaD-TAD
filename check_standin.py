"""check_standin.py -- stand-in DNA clash check against the CURRENT ring waypoints.
Exit code 1 if any far DNA/histone overlap remains (used by REPAIR_TETHERED.bat)."""
import argparse, sys, numpy as np, dna_proxy as DP
ap = argparse.ArgumentParser()
ap.add_argument("--frames", default="frames_tethered.npz")
ap.add_argument("--tethers", default="data/tethers.json")
a = ap.parse_args()
F = np.load(a.frames)
tw, ex = DP.load_tw(a.tethers)
c = DP.clashes(F["origins"], F["rotations"], DP.local_templates(), ex, tw)
for e1, e2, k, d in c:
    print(f"  {DP.label(e1)} | {DP.label(e2)}  {k}  {d:.1f} A")
print(f"stand-in check with current waypoints: {len(c)} clash(es)")
sys.exit(1 if c else 0)
