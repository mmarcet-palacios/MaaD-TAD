"""
Cheap DNA/histone stand-in for pre-build clash checks (used by relax_frames.py).

Each nucleosome contributes
  core DNA axis points : midpoints of paired C1' atoms of 1KX5's wrapped DNA,
                         every 5 bp (~30 points), in the nucleosome frame
  histone points       : every 4th histone C-alpha (~240 points)
Each linker contributes
  axis points          : a cubic Bezier from the exit bp to the entry bp, with
                         handles along the exit/entry helix axes (|span|/3),
                         sampled every ~10 A -- a smooth stand-in for the
                         elastic path the stage-2 solver will find.
Element ids: core k -> 2k, linker k -> 2k+1; "far" means |id difference| > 4,
i.e. more than two nucleosomes apart.
"""
import numpy as np
from build_zigzag import canonical_1kx5
from build_array9b import kx5_parts

AXIS_MIN = 22.0      # duplex axis to duplex axis: 2 x 8.9 A (P radius) + ~4 A gap
AXIS_CA_MIN = 11.0   # duplex axis to histone C-alpha


def local_templates():
    parts = kx5_parts()
    I = {sq: at for ch, rn, sq, at in parts["dnaI"]}
    J = {sq: at for ch, rn, sq, at in parts["dnaJ"]}
    ax = [0.5 * (I[k]["C1'"] + J[-k]["C1'"]) for k in range(-70, 71, 5) if k in I and -k in J
          and "C1'" in I[k] and "C1'" in J[-k]]
    ca = [at["CA"] for ch, rn, sq, at in parts["protein"] if "CA" in at][::4]
    le, lx = canonical_1kx5()
    return np.array(ax), np.array(ca), le, lx


def linker_curve(p0, t0, p3, t3, step=10.0):
    span = np.linalg.norm(p3 - p0)
    h = span / 3.0
    p1, p2 = p0 + t0 * h, p3 - t3 * h
    n = max(3, int(span / step) + 1)
    s = np.linspace(0, 1, n)[1:-1, None]
    return ((1 - s) ** 3) * p0 + 3 * ((1 - s) ** 2) * s * p1 + 3 * (1 - s) * s ** 2 * p2 + s ** 3 * p3


def bez(p0, t0, p3, t3, step=10.0):
    return linker_curve(p0, t0, p3, t3, step)


def build(O, R, tmpl, tw=None):
    """tw: {linker: (waypoint, dir)} -- tether linkers are drawn as two curves
    meeting at their ring waypoint, as stage 2 builds them."""
    ax, ca, le, lx = tmpl
    tw = tw or {}
    n = len(O)
    P, E, K = [], [], []          # points, element id, kind (0 axis, 1 CA)
    for k in range(n):
        P.append(np.einsum("ij,nj->ni", R[k], ax) + O[k]); E.append(np.full(len(ax), 2 * k)); K.append(np.zeros(len(ax), int))
        P.append(np.einsum("ij,nj->ni", R[k], ca) + O[k]); E.append(np.full(len(ca), 2 * k)); K.append(np.ones(len(ca), int))
        if k < n - 1:
            p0 = R[k] @ lx[0] + O[k]; t0 = R[k] @ lx[1][:, 2]
            p3 = R[k + 1] @ le[0] + O[k + 1]; t3 = R[k + 1] @ le[1][:, 2]
            if k in tw:
                W, dW = tw[k]
                c = np.vstack([bez(p0, t0, W, dW), W[None], bez(W, dW, p3, t3)])
            else:
                c = linker_curve(p0, t0, p3, t3)
            P.append(c); E.append(np.full(len(c), 2 * k + 1)); K.append(np.zeros(len(c), int))
    return np.vstack(P), np.concatenate(E), np.concatenate(K)


def load_tw(path):
    import json, os
    if not path or not os.path.exists(path):
        return {}, []
    T = json.load(open(path))
    tw = {int(L["linker"]): (np.array(L["waypoint"]), np.array(L["dir"])) for t in T for L in t["linkers"]}
    return tw, [(t["linkers"][0]["linker"], t["linkers"][1]["linker"]) for t in T]


def clashes(O, R, tmpl, exempt_links=(), tw=None):
    """Far element pairs whose stand-ins overlap. Returns list of (elem_a, elem_b, kind, dist)."""
    from scipy.spatial import cKDTree
    P, E, K = build(O, R, tmpl, tw)
    t = cKDTree(P)
    pr = t.query_pairs(AXIS_MIN, output_type="ndarray")
    ex = set()
    for a, b in exempt_links:
        ex.add((min(2 * a + 1, 2 * b + 1), max(2 * a + 1, 2 * b + 1)))
    out = {}
    for i, j in pr:
        ei, ej = E[i], E[j]
        if abs(ei - ej) <= 4 or (K[i] == 1 and K[j] == 1):
            continue
        d = np.linalg.norm(P[i] - P[j])
        kind = "axis-axis" if K[i] == 0 and K[j] == 0 else "axis-CA"
        if kind == "axis-CA" and d >= AXIS_CA_MIN:
            continue
        key = (min(ei, ej), max(ei, ej))
        if key in ex:
            continue
        if key not in out or d < out[key][1]:
            out[key] = (kind, d)
    return [(a, b, k, d) for (a, b), (k, d) in sorted(out.items())]


def label(e):
    return f"core {e // 2}" if e % 2 == 0 else f"linker {e // 2}"
