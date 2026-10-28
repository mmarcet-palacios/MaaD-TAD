"""
Measure the geometry of the b2_onelayer array against three stated claims:

  1. it is NOT a flat zigzag
  2. each core sits about as far from the others as it can
  3. nucleosomes always face the centre of the fibre trajectory

and assess the fourth, mechanical claim: that this arrangement can be
stretched without forcing tight linker curvature.
"""
import numpy as np
import pymol2
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.interpolate import splprep, splev

PSE = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "data", "b2_onelayer.pse")


def load():
    p = pymol2.PyMOL(); p.start(); cmd = p.cmd
    cmd.load(PSE)
    obj = cmd.get_names("objects")[0]

    prot = {}
    for ch in cmd.get_chains(obj):
        sel = f"{obj} and chain {ch} and polymer.protein and name CA"
        if cmd.count_atoms(sel) > 40:
            prot[ch] = cmd.get_coords(sel)

    dna = {}
    for ch in ("A", "B"):
        sel = f"{obj} and chain {ch} and polymer.nucleic and name C1'"
        if cmd.count_atoms(sel) == 0:
            continue
        idx = []
        cmd.iterate(sel, "idx.append(resi)", space={"idx": idx})
        xyz = cmd.get_coords(sel)
        order = np.argsort([int(i) for i in idx])
        dna[ch] = (np.array([int(idx[i]) for i in order]), xyz[order])
    p.stop()
    return prot, dna


def nucleosome_frames(prot, dna):
    cents = np.array([v.mean(axis=0) for v in prot.values()])
    n_nuc = len(cents) // 8
    lab = fcluster(linkage(cents, "ward"), n_nuc, criterion="maxclust")
    oct_c = np.array([cents[lab == k].mean(axis=0) for k in range(1, n_nuc + 1)])

    resi, xyz = dna["A"]
    d = np.linalg.norm(xyz[:, None, :] - oct_c[None, :, :], axis=-1)
    near = np.argmin(d, axis=1)
    mind = d[np.arange(len(xyz)), near]
    assign = np.where(mind < 60.0, near, -1)

    runs, cur, start = [], assign[0], 0
    for i in range(1, len(assign) + 1):
        if i == len(assign) or assign[i] != cur:
            if cur >= 0 and (i - start) >= 90:
                runs.append((cur, start, i))
            if i < len(assign):
                cur, start = assign[i], i

    frames, meta = [], []
    for nuc, s, e in runs:
        pts = xyz[s:e]
        c = pts.mean(axis=0)
        _, _, vt = np.linalg.svd(pts - c)
        z = vt[2]                                  # disc normal
        dyad = xyz[(s + e) // 2]
        y = dyad - c; y -= np.dot(y, z) * z; y /= np.linalg.norm(y)
        frames.append((c, np.column_stack([np.cross(y, z), y, z])))
        meta.append({"wrapped_bp": e - s, "first": int(resi[s]), "last": int(resi[e - 1]),
                     "entry": pts[0], "exit": pts[-1]})
    return frames, meta, oct_c


def main():
    prot, dna = load()
    frames, meta, oct_c = nucleosome_frames(prot, dna)
    C = np.array([f[0] for f in frames])
    Z = np.array([f[1][:, 2] for f in frames])
    n = len(C)

    print(f"{len(prot)} histone chains, {len(dna['A'][0])} bp per strand "
          f"-> {n} nucleosomes ordered along the DNA\n")
    linkers = [meta[i + 1]["first"] - meta[i]["last"] - 1 for i in range(n - 1)]
    print(f"wrapped bp   : {[m['wrapped_bp'] for m in meta]}")
    print(f"linker bp    : {linkers}")
    nrl = np.mean([meta[i]['wrapped_bp'] for i in range(n)]) + np.mean(linkers)
    print(f"implied NRL  : {nrl:.0f} bp\n")

    # ---- claim 1: not a flat zigzag ------------------------------------
    cen = C.mean(axis=0)
    u, s, vt = np.linalg.svd(C - cen)
    out_of_plane = np.abs((C - cen) @ vt[2])
    print("CLAIM 1 -- is it flat?")
    print(f"  principal extents (A)      : {s[0]/np.sqrt(n):.1f}  {s[1]/np.sqrt(n):.1f}  {s[2]/np.sqrt(n):.1f}")
    print(f"  out-of-plane deviation     : {out_of_plane.mean():.1f} +/- {out_of_plane.std():.1f} A"
          f"  (max {out_of_plane.max():.1f})")
    flat = s[2] / s[0]
    print(f"  3rd/1st extent ratio       : {flat:.3f}"
          f"   -> {'FLAT (planar)' if flat < 0.08 else 'NOT FLAT -- genuinely 3D'}\n")

    # ---- claim 2: maximally separated ----------------------------------
    print("CLAIM 2 -- are cores pushed as far apart as possible?")
    d1 = np.linalg.norm(np.diff(C, axis=0), axis=1)
    print(f"  i -> i+1 centre distance   : {d1.mean():.1f} +/- {d1.std():.1f} A")
    for k in (2, 3, 4):
        if n > k:
            dk = np.linalg.norm(C[k:] - C[:-k], axis=1)
            print(f"  i -> i+{k} centre distance   : {dk.mean():.1f} +/- {dk.std():.1f} A")
    allp = np.linalg.norm(C[:, None] - C[None], axis=-1)
    np.fill_diagonal(allp, np.inf)
    nn = allp.min(axis=1)
    print(f"  nearest neighbour (any i)  : {nn.mean():.1f} +/- {nn.std():.1f} A")
    theo = np.array([lk * 3.34 + 110.0 for lk in linkers])   # extended linker + disc
    print(f"  max possible for the linker: {theo.mean():.1f} A  "
          f"-> achieving {100*d1.mean()/theo.mean():.0f}% of maximum separation\n")

    # ---- claim 3: facing the fibre centre ------------------------------
    print("CLAIM 3 -- do nucleosomes face the trajectory centre?")
    tck, _ = splprep([C[:, 0], C[:, 1], C[:, 2]], s=len(C) * 50.0, k=min(3, n - 1))
    tt = np.linspace(0, 1, 400)
    curve = np.array(splev(tt, tck)).T
    radial, ang_norm, ang_dyad = [], [], []
    for i in range(n):
        j = np.argmin(np.linalg.norm(curve - C[i], axis=1))
        r = curve[j] - C[i]
        nr = np.linalg.norm(r)
        if nr < 1e-6:
            continue
        r /= nr
        radial.append(nr)
        a = np.degrees(np.arccos(abs(np.clip(np.dot(Z[i], r), -1, 1))))
        ang_norm.append(a)
        y_i = frames[i][1][:, 1]
        ang_dyad.append(np.degrees(np.arccos(abs(np.clip(np.dot(y_i, r), -1, 1)))))
    radial = np.array(radial); ang_norm = np.array(ang_norm); ang_dyad = np.array(ang_dyad)
    print(f"  distance from trajectory   : {radial.mean():.1f} +/- {radial.std():.1f} A")
    print(f"  disc NORMAL vs radial dir  : {ang_norm.mean():.1f} +/- {ang_norm.std():.1f} deg")
    print(f"  dyad AXIS  vs radial dir   : {ang_dyad.mean():.1f} +/- {ang_dyad.std():.1f} deg")
    print("   (0 deg = that axis points at the trajectory; 90 deg = perpendicular to it)")
    if ang_norm.mean() < 35:
        print("   -> disc FACE points at the trajectory centre")
    elif ang_norm.mean() > 55:
        print("   -> disc EDGE faces the trajectory; the face looks outward/along")
    else:
        print("   -> intermediate")
    print()

    # ---- claim 4: stretch headroom -------------------------------------
    print("CLAIM 4 -- stretch headroom in the linkers")
    spans, contours = [], []
    for i in range(n - 1):
        span = np.linalg.norm(meta[i + 1]["entry"] - meta[i]["exit"])
        spans.append(span); contours.append(linkers[i] * 3.34)
    spans = np.array(spans); contours = np.array(contours)
    frac = spans / contours
    print(f"  linker exit->entry span    : {spans.mean():.1f} +/- {spans.std():.1f} A")
    print(f"  linker contour available   : {contours.mean():.1f} A")
    print(f"  extension  (span/contour)  : {frac.mean():.2f} +/- {frac.std():.2f}")
    print(f"  headroom before taut       : {100*(1-frac.mean()):.0f}% of contour still in reserve")
    return frames, meta, C, Z, linkers


if __name__ == "__main__":
    main()
