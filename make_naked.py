"""
Layer 3, step 1: make the DNA naked at each cohesin site and its mirror.

    python make_naked.py                         -> data/nucleosomes_tethered.bed
    python make_naked.py --clear 4 --frames f.npz   (also remove nucleosomes that
                                                     still collide in tether zones)

A nucleosome is removed if its 147 bp core overlaps the cohesin ChIP peak, or
the same-width window centred on the mirror copy (both duplexes pass through
the ring, so both must be free of histones). Removed nucleosomes simply
lengthen the linker between their neighbours; genomic numbering is unchanged.

--clear N: after a stage-1 run, nucleosomes inside a tether zone (within
--zone-kb of a site) that still collide (<110 A) are removed too, at most N per
site side, closest to the site first. Every removal is listed with its reason.
"""
import argparse, json
import numpy as np
from scipy.spatial import cKDTree


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default="data/nucleosomes_GRCh38_resolved.bed")
    ap.add_argument("--restraints", default="data/homology_restraints.json")
    ap.add_argument("--out", default="data/nucleosomes_tethered.bed")
    ap.add_argument("--frames", help="stage-1 frames built from --out, for --clear")
    ap.add_argument("--clear", type=int, default=0)
    ap.add_argument("--zone-kb", type=float, default=3.0)
    ap.add_argument("--min-window", type=int, default=400,
                    help="naked stretch at least this long (bp), centred on each site, so both "
                         "duplexes have room to thread through the ring")
    a = ap.parse_args()

    R = json.load(open(a.restraints))
    header, rows = [], []
    for l in open(a.bed):
        (header if l.startswith("#") else rows).append(l)
    dy = np.array([int(r.split("\t")[6]) for r in rows])
    sites = []
    for t in R["tethers"]:
        lo, hi = t["cohesin"]; w = hi - lo
        W = max(w, a.min_window); c = (lo + hi) // 2
        sites.append((t["id"], "cohesin", c - W // 2, c + W - W // 2))
        sites.append((t["id"], "mirror", t["mirror"] - W // 2, t["mirror"] + W - W // 2))
    drop = {}
    for tid, kind, lo, hi in sites:
        hit = np.where((dy + 73 >= lo) & (dy - 73 <= hi))[0]
        for i in hit:
            drop[int(i)] = f"{tid} {kind} {lo:,}-{hi:,} ({hi-lo} bp)"

    if a.clear and a.frames:
        F = np.load(a.frames); O = F["origins"]; fdy = F["dyad_pos"]
        p = cKDTree(O).query_pairs(110.0, output_type="ndarray")
        p = p[np.abs(p[:, 0] - p[:, 1]) > 2]
        bad = set(p.ravel().tolist())
        for tid, kind, lo, hi in sites:
            c = 0.5 * (lo + hi)
            cand = [i for i in bad if abs(fdy[i] - c) < a.zone_kb * 1000]
            for side in (-1, 1):
                sc = sorted([i for i in cand if np.sign(fdy[i] - c) == side], key=lambda i: abs(fdy[i] - c))
                for i in sc[:a.clear]:
                    j = int(np.where(dy == fdy[i])[0][0])
                    drop.setdefault(j, f"{tid} {kind}: still colliding at the pinch, "
                                       f"{abs(fdy[i]-c)/1e3:.1f} kb from site")
    keep = [r for i, r in enumerate(rows) if i not in drop]
    with open(a.out, "w") as f:
        f.writelines(header)
        f.write(f"# layer 3: {len(drop)} nucleosomes removed at cohesin tethers (make_naked.py)\n")
        f.writelines(keep)
    print(f"removed {len(drop)} of {len(rows)} nucleosomes -> {a.out} ({len(keep)} remain)")
    for i in sorted(drop):
        print(f"  dyad {dy[i]:,}  <- {drop[i]}")
    kd = np.array([int(r.split('\t')[6]) for r in keep])
    nrl = np.diff(kd)
    print(f"  longest linker now {nrl.max()-147} bp (NRL {nrl.max()})")


if __name__ == "__main__":
    main()
