"""
Option 2: absorb the excess fibre contour with a PLANAR serpentine (local
backfolding) instead of a helical coil.

Why this might help: a coil needs clearance in EVERY radial direction; a flat
serpentine needs width in one direction only and just the fibre radius in the
other. Each section's fold plane is oriented to put its thin direction toward
the nearest crowding guide strands.

What cannot change: the fibre's minimum bend radius. At adv ~66.9 A/nucleosome
and a steering budget of theta deg/nuc, R_min = adv / theta (8 deg -> ~480 A).
The serpentine is built from arcs of exactly that radius plus straight legs.

    python backfold_axis.py --bed ... --guide ... --gaps ... [--max-steer 8]
scans domain diameters and reports, per scale: max steering, how well the
per-section contour demand is met, and the minimum fibre-fibre clearance
between non-local parts of the axis (must exceed 2 x fibre radius = 220 A).
"""
import argparse
import numpy as np
from scipy.spatial import cKDTree

import fibre_adaptive as fa
import fibre_builder as fb
from build_zigzag import canonical_1kx5, place
from stage1_frames import load_dyads, load_guide, ROI_START, BIN, ADV_PER_STEP

FIBRE_R = 110.0
SPP = 40                 # guide samples per section


def demand(bed, guide):
    pos, _ = load_dyads(bed)
    nrl = np.diff(pos)
    le, lx = canonical_1kx5()
    cache = {}

    def sl(b):
        b = int(np.clip(b, 10, 200))
        if b not in cache:
            cache[b] = fa.solve_length(b, place, le, lx)[0]
        return cache[b]
    Ls = np.array([sl(int(x) - 147) for x in nrl])
    adv = Ls * ADV_PER_STEP
    sec_of = np.clip((pos[:-1] - ROI_START) // BIN, 0, len(guide) - 2)
    need = np.array([adv[sec_of == k].sum() for k in range(len(guide) - 1)])
    return need, adv


def dense_guide(guide):
    pts, sid = [], []
    for k in range(len(guide) - 1):
        t = np.linspace(0, 1, SPP, endpoint=False)
        pts.append(guide[k] + t[:, None] * (guide[k + 1] - guide[k]))
        sid.append(np.full(SPP, k))
    pts.append(guide[-1][None]); sid.append([len(guide) - 2])
    return np.vstack(pts), np.concatenate(sid)


def fold_direction(pts, T, N1, N2, reach, excl, smooth=121):
    """Lateral unit vector u(s) in the normal plane, chosen PERPENDICULAR to
    the local crowding so the serpentine's thin side faces the neighbours."""
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    tree = cKDTree(pts)
    psi = np.zeros(len(pts))
    for i, nb in enumerate(tree.query_ball_point(pts, reach)):
        nb = np.array([j for j in nb if abs(s[j] - s[i]) > excl], int)
        if not len(nb):
            psi[i] = np.nan; continue
        d = pts[nb] - pts[i]
        a, b = d @ N1[i], d @ N2[i]
        w = 1.0 / np.maximum(np.hypot(a, b), 50.0) ** 2
        C = np.array([[np.sum(w*a*a), np.sum(w*a*b)], [np.sum(w*a*b), np.sum(w*b*b)]])
        ev, V = np.linalg.eigh(C)
        v = V[:, 0]                       # least-crowded in-plane direction
        psi[i] = np.arctan2(v[1], v[0])
    ok = ~np.isnan(psi)
    psi[~ok] = np.interp(np.where(~ok)[0], np.where(ok)[0], psi[ok]) if ok.any() else 0.0
    # u and -u are the same fold plane: work with 2*psi, unwrap, smooth
    z = np.unwrap(2 * psi)
    k = np.ones(smooth) / smooth
    z = np.convolve(np.pad(z, smooth // 2, mode="edge"), k, mode="valid")
    psi = z / 2
    return np.cos(psi)[:, None] * N1 + np.sin(psi)[:, None] * N2


def serpentine_2d(xlen, m_of_x, R, ds=10.0):
    """Planar curve in (x along guide, y lateral), built from arcs of radius R
    and straight legs, whose local contour/advance ratio tracks m(x).
    For m <= pi/2: arcs only, swinging between headings +-Theta, Theta/sin = m.
    For m >  pi/2: Theta = 90 deg, plus lateral legs of length R*(2m - pi)."""
    def theta_leg(m):
        if m <= 1.0005:
            return 0.0, 0.0
        if m < np.pi / 2:
            lo, hi = 1e-4, np.pi / 2
            for _ in range(50):
                mid = (lo + hi) / 2
                (lo, hi) = (mid, hi) if mid / np.sin(mid) < m else (lo, mid)
            return lo, 0.0
        return np.pi / 2, R * (2 * m - np.pi)

    x, y, h = 0.0, 0.0, 0.0
    out = [(x, y)]
    sign = 1.0
    while x < xlen:
        Th, leg = theta_leg(m_of_x(x))
        if Th == 0.0:                      # no excess: go straight, re-centre
            h = -np.clip(y, -R, R) / (4 * R)
            x += ds * np.cos(h); y += ds * np.sin(h); out.append((x, y)); continue
        target = sign * Th
        # arc toward target heading
        while (target - h) * sign > 1e-9 and x < xlen:
            dh = min(ds / R, (target - h) * sign) * sign
            h += dh
            x += ds * np.cos(h); y += ds * np.sin(h); out.append((x, y))
        # straight leg; stop at +-leg/2 lateral so the wave stays centred
        ytarget = sign * leg / 2
        while leg > 0 and (ytarget - y) * sign > 0 and x < xlen:
            x += ds * np.cos(h); y += ds * np.sin(h); out.append((x, y))
        sign = -sign
    return np.array(out)


def build_axis(guide, need, R, reach=2500.0):
    pts, sid = dense_guide(guide)
    T, N1, N2 = fb.parallel_frames(pts)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    chord = np.linalg.norm(np.diff(guide, axis=0), axis=1)
    u = fold_direction(pts, T, N1, N2, reach, excl=3 * R)
    sec_start = np.concatenate([[0], np.cumsum(chord)])
    m = np.maximum(need / np.maximum(chord, 1e-9), 1.0)
    # smooth m a little (3 sections), preserving total demand
    mk = np.convolve(np.pad(m, 1, mode="edge"), np.ones(3) / 3, mode="valid")
    scale = 1.0
    for _ in range(6):                     # match total contour
        mm = mk * scale
        mf = lambda x: mm[int(np.clip(np.searchsorted(sec_start, x) - 1, 0, len(mm) - 1))]
        c2 = serpentine_2d(s[-1], mf, R)
        L = np.linalg.norm(np.diff(c2, axis=0), axis=1).sum()
        scale *= need.sum() / L
    # map (x, y) onto the guide: g(x) + y * u(x)
    X = np.clip(c2[:, 0], 0, s[-1])
    idx = np.clip(np.searchsorted(s, X) - 1, 0, len(pts) - 2)
    f = ((X - s[idx]) / np.maximum(s[idx + 1] - s[idx], 1e-9))[:, None]
    g = pts[idx] * (1 - f) + pts[idx + 1] * f
    uu = u[idx] * (1 - f) + u[idx + 1] * f
    uu /= np.linalg.norm(uu, axis=1, keepdims=True)
    axis = g + c2[:, 1:2] * uu
    return axis, X, sec_start


def evaluate(axis, X, sec_start, need, adv_mean):
    seg = np.linalg.norm(np.diff(axis, axis=0), axis=1)
    sc = np.concatenate([[0], np.cumsum(seg)])
    # curvature -> steering per nucleosome
    t = np.diff(axis, axis=0) / seg[:, None]
    ang = np.arccos(np.clip(np.sum(t[1:] * t[:-1], axis=1), -1, 1))
    kappa = ang / (0.5 * (seg[1:] + seg[:-1]))
    # average over one nucleosome's advance (steering is per nucleosome)
    w = max(1, int(round(adv_mean / np.median(seg))))
    kap_n = np.convolve(kappa, np.ones(w) / w, mode="same")
    steer = np.degrees(kap_n * adv_mean)
    # per-section contour actually delivered
    sec = np.clip(np.searchsorted(sec_start, X[:-1]) - 1, 0, len(need) - 1)
    got = np.bincount(sec, weights=seg, minlength=len(need))
    err = np.abs(got - need) / np.maximum(need, 1)
    # non-local clearance (contour separation > 3*fibre diameter)
    tree = cKDTree(axis)
    pairs = tree.query_pairs(2 * FIBRE_R, output_type="ndarray")
    bad = pairs[np.abs(sc[pairs[:, 0]] - sc[pairs[:, 1]]) > 660]
    d, j = tree.query(axis, k=40)
    far = np.abs(sc[j] - sc[:, None]) > 660
    dmin = np.where(far, d, np.inf).min()
    clash_secs = len(np.unique(sec[np.clip(bad[:, 0], 0, len(sec) - 1)])) if len(bad) else 0
    return dict(len_um=sc[-1] / 1e4, steer_med=np.median(steer),
                steer95=np.percentile(steer, 95), steer_max=steer.max(),
                err_med=np.median(err), dmin=dmin, clash_secs=clash_secs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True); ap.add_argument("--guide", required=True)
    ap.add_argument("--gaps", required=True)
    ap.add_argument("--max-steer", type=float, default=8.0)
    ap.add_argument("--domains", default="500,600,700,800")
    a = ap.parse_args()
    for D in [float(x) for x in a.domains.split(",")]:
        guide = load_guide(a.guide, a.gaps, D)
        need, adv = demand(a.bed, guide)
        R = adv.mean() / np.radians(a.max_steer)
        axis, X, ss = build_axis(guide, need, R)
        r = evaluate(axis, X, ss, need, adv.mean())
        print(f"{D:5.0f} nm  R={R:4.0f} A  axis {r['len_um']:.2f} um  "
              f"steer med/95/max {r['steer_med']:.1f}/{r['steer95']:.1f}/{r['steer_max']:.1f} deg  "
              f"arc err {r['err_med']:.0%}  min clearance {r['dmin']:.0f} A  "
              f"clashing sections {r['clash_secs']}/325", flush=True)


if __name__ == "__main__":
    main()
