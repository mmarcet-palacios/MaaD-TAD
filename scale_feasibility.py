"""
Is there a domain diameter at which the model is both reachable and sterically
valid? Run this BEFORE stage 1 -- it decides whether a build is possible at all.

Two constraints pull in opposite directions as the domain is scaled.

REACHABILITY. Each guide segment must be given at least as much fibre axis as
its bead gap needs, i.e. meander >= 1. Meander scales as 1/D, so growing the
domain eventually leaves sections the fibre cannot span.

STERICS. The fibre is ~220 A thick and any coil wraps further out, so the guide
curve's own self-approach must exceed 2*(coil radius + fibre radius). Self-
approach scales as D, so a SMALL domain leaves no room. And the coil radius
cannot be made arbitrarily small, because axis curvature goes as 1/r: a tight
coil is a sharply bent fibre, and that bend adds to every linker.

If the two windows do not overlap there is no valid scale, and that is a
statement about the guide curve and the nucleosome content together -- not a
solver failure.
"""
import argparse, sys
import numpy as np
sys.path.insert(0, "/tmp/dna_poc/bpstep")
import fibre_adaptive as fa
from build_zigzag import canonical_1kx5, place

ROI, BIN = 12_771_273, 2000
ADV = 0.3023
FIBRE_RADIUS_A = 110.0
ADV_NM_PER_NUC = 66.9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--guide", required=True)
    ap.add_argument("--gaps", required=True)
    ap.add_argument("--max-steer", type=float, default=8.0,
                    help="deg per nucleosome the linkers can absorb")
    args = ap.parse_args()

    le, lx = canonical_1kx5()
    pos = []
    for l in open(args.bed):
        if not l.startswith("#"):
            pos.append(int(l.split("\t")[6]))
    pos = np.array(pos); nrl = np.diff(pos)
    cache = {}

    def sl(b):
        b = int(np.clip(b, 10, 200))
        if b not in cache:
            cache[b] = fa.solve_length(b, place, le, lx)[0]
        return cache[b]

    adv = np.array([sl(int(x) - 147) for x in nrl]) * ADV

    xyz = []
    for l in open(args.guide):
        p = l.split()
        if len(p) >= 4:
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
    xyz = np.array(xyz)
    gaps = {int(l.strip()) for l in open(args.gaps) if l.strip()}
    keep = ~np.isin(np.arange(len(xyz)) + 163, list(gaps))
    dmax = np.linalg.norm(xyz[keep][:, None] - xyz[keep][None], axis=-1).max()
    g = xyz * (5000.0 / dmax)
    chord = np.linalg.norm(np.diff(g, axis=0), axis=1)
    sec_of = np.clip((pos[:-1] - ROI) // BIN, 0, len(xyz) - 2)
    need = np.array([adv[sec_of == k].sum() for k in range(len(xyz) - 1)])
    ok = need > 0
    ratio = need[ok] / chord[ok]
    D = np.linalg.norm(g[:, None] - g[None], axis=-1)
    i = np.arange(len(g))
    clr = D[np.abs(i[:, None] - i[None, :]) > 3].min()

    print(f"nucleosomes {len(pos)}, fibre axis required {adv.sum()/1e4:.2f} um")
    print(f"need/chord at 500 nm: min {ratio.min():.2f}, median {np.median(ratio):.2f}")
    print(f"guide self-approach at 500 nm: {clr:.0f} A\n")
    print(f"{'domain':>8} {'unreachable':>12} {'coil r fits':>12} "
          f"{'coil r needed':>14} {'verdict':>12}")
    feas = []
    for Dn in range(400, 1601, 50):
        m = ratio * 500.0 / Dn
        unreach = int((m < 1.0).sum())
        rmax = clr * (Dn / 500.0) / 2 - FIBRE_RADIUS_A
        mg = float(np.median(np.maximum(m, 1.0001)))
        rmin = (ADV_NM_PER_NUC / np.radians(args.max_steer) * (1 - 1 / mg ** 2)
                if mg > 1.0001 else 0.0)
        good = unreach == 0 and rmax >= rmin
        if good:
            feas.append(Dn)
        if Dn % 100 == 0:
            v = "OK" if good else ("unreachable" if unreach else "no room")
            print(f"{Dn:6d}nm {unreach:9d}/{len(ratio)} {rmax:10.0f} A "
                  f"{rmin:12.0f} A {v:>12}")
    print()
    if feas:
        print(f"FEASIBLE domain diameters: {min(feas)}-{max(feas)} nm")
    else:
        print("NO FEASIBLE DOMAIN DIAMETER.")
        print("  Reachability caps the domain; sterics require a larger one.")
        print("  The guide curve and the nucleosome content are inconsistent")
        print("  under a coiled-fibre model. See RUN.md for the three options.")


if __name__ == "__main__":
    main()
