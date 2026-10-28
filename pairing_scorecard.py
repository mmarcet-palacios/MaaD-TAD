"""
Nucleosome-level scorecard for homology pairing and tethers.

    python pairing_scorecard.py --frames frames.npz

For every block in data/homology_restraints.json, sampled every 1 kb along the
arms (positions interpolated between nucleosome centres by genomic position):
  spacing     median arm-to-arm distance and its variation (CV) -> parallel
  direction   angle between the arms' local directions (0 = perfectly parallel
              for direct, perfectly antiparallel for inverted)
  register    for each point on arm A, the genomic offset between the point of
              arm C that is actually closest and its homologous partner
              (0 bp = homologous positions face each other)
Tethers: distance between the cohesin site and its mirror.
Collisions: nucleosome pairs < 110 A (>2 apart), split into inside / outside
the tether zones (tether zones are expected to clash until layer 3 removes the
nucleosomes over the cohesin peaks).
"""
import argparse, json
import numpy as np
from scipy.spatial import cKDTree


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--restraints", default="data/homology_restraints.json")
    ap.add_argument("--pinch-kb", type=float, default=6.0)
    a = ap.parse_args()
    F = np.load(a.frames); O = F["origins"]; dy = F["dyad_pos"].astype(float)
    R = json.load(open(a.restraints))
    pos = lambda g: np.column_stack([np.interp(g, dy, O[:, k]) for k in range(3)])
    T = [(t["id"], 0.5 * (t["cohesin"][0] + t["cohesin"][1]), t["mirror"]) for t in R["tethers"]]

    print(f"{'block':<5}{'type':<9}{'ident':>6}{'spacing':>9}{'CV':>6}{'dir':>7}{'register':>12}")
    for b in R["blocks"]:
        A0, A1 = b["armA"]; C0, C1 = b["armC"]
        n = int((A1 - A0) / 1000) + 1
        t = np.linspace(0, 1, n)
        gA = A0 + t * (A1 - A0)
        gC = C0 + t * (C1 - C0) if b["type"] == "direct" else C1 - t * (C1 - C0)
        # drop samples inside any tether pinch zone
        keep = np.ones(n, bool)
        for _, gc, gm in T:
            for g in (gA, gC):
                keep &= np.minimum(np.abs(g - gc), np.abs(g - gm)) > a.pinch_kb * 1000
        if keep.sum() < 3:
            print(f"{b['id']:<5}{b['type']:<9}{b['ident']:6.1f}   (arm lies within a tether pinch zone)")
            continue
        PA, PC = pos(gA), pos(gC)
        d = np.linalg.norm(PA - PC, axis=1)[keep]
        h = 1000.0
        tA = pos(gA + h) - pos(gA - h); tC = pos(gC + h) - pos(gC - h)
        tA /= np.linalg.norm(tA, axis=1)[:, None]; tC /= np.linalg.norm(tC, axis=1)[:, None]
        sgn = 1 if b["type"] == "direct" else -1
        ang = np.degrees(np.arccos(np.clip(sgn * np.sum(tA * tC, 1), -1, 1)))[keep]
        # register: closest point on arm C (dense, 50 bp) to each arm-A point
        gCd = np.arange(C0, C1 + 1, 50.0); PCd = pos(gCd)
        near = gCd[cKDTree(PCd).query(PA[keep])[1]]
        reg = near - gC[keep]
        print(f"{b['id']:<5}{b['type']:<9}{b['ident']:6.1f}{np.median(d):8.0f}A{d.std()/d.mean():6.1%}"
              f"{np.median(ang):6.1f}d{np.median(np.abs(reg)):8.0f} bp")
    print()
    for tid, gc, gm in T:
        print(f"tether {tid}: cohesin-site to mirror {np.linalg.norm(pos([gc]) - pos([gm])):.0f} A "
              f"(nucleosome centres; the naked DNA + ring close the rest in layer 3)")
    p = cKDTree(O).query_pairs(110.0, output_type="ndarray")
    p = p[np.abs(p[:, 0] - p[:, 1]) > 2]
    zone = np.zeros(len(O), bool)
    for _, gc, gm in T:
        zone |= (np.abs(dy - gc) < a.pinch_kb * 1000) | (np.abs(dy - gm) < a.pinch_kb * 1000)
    inz = zone[p[:, 0]] & zone[p[:, 1]]
    dist = np.linalg.norm(O[p[:, 0]] - O[p[:, 1]], axis=1)
    print(f"\ncollisions (<110 A, >2 apart): {len(p)} total -> {inz.sum()} inside tether zones, "
          f"{(~inz).sum()} outside" + (f" (closest outside {dist[~inz].min():.0f} A)" if (~inz).any() else ""))


if __name__ == "__main__":
    main()
