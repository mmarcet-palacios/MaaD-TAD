"""
Stage 1: nucleosome frames along the meandering fibre axis.

This is the only genuinely new code in the full build, and the only place where
an error is structural rather than cosmetic. Validate it before writing atoms.

  guide curve (326 beads, arbitrary units)
      -> scaled to the 500 nm domain
      -> per-section meander from the nucleosome counts
      -> coiled axis (continuous phase across section joins)
      -> nucleosome frames, steered to follow that axis

Usage:
    python3 stage1_frames.py --out frames.npz
"""
import argparse, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np

import fibre_adaptive as fa
import fibre_builder as fb
import closure as cl
from build_zigzag import canonical_1kx5, place

ROI_START, BIN, FIRST_SECTION = 12_771_273, 2000, 163
DOMAIN_NM = 500.0
COIL_RADIUS_A = 800.0
ADV_PER_STEP = 0.3023      # axis advance / step length, fixed by the b2 rotation
SAMPLES_PER_SECTION = 40


def load_dyads(path):
    pos, src = [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pos.append(int(f[6])); src.append(f[7])
    return np.array(pos), np.array(src)


def load_guide(path, gaps_path, domain_nm=DOMAIN_NM):
    xyz = []
    for line in open(path):
        p = line.split()
        if len(p) >= 4:
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
    xyz = np.array(xyz)
    gaps = {int(l.strip()) for l in open(gaps_path) if l.strip()}
    sec = np.arange(len(xyz)) + FIRST_SECTION
    real = xyz[~np.isin(sec, list(gaps))]
    dmax = np.linalg.norm(real[:, None, :] - real[None, :, :], axis=-1).max()
    return xyz * (domain_nm * 10.0 / dmax)


def coiled_axis(guide, need, coil_r=COIL_RADIUS_A, spp=SAMPLES_PER_SECTION):
    """Helical meander around the guide, phase continuous across joins.

    The phase rate is solved NUMERICALLY per section so the coil's measured arc
    length equals what that section needs. The analytic pitch formula
    (2*pi*r/p = sqrt(m^2-1)) assumes a helix around a STRAIGHT axis; around a
    curving one the rotation-minimising normal frame contributes extra length,
    and using the formula directly overshot the target by 15%.

    Phase is ADVANCED across section boundaries rather than recomputed, because
    recomputing restarts it at every bead and puts a kink there.
    """
    # Interpolate the guide with a C2 cubic spline through the beads (chord-
    # length parameterised). Straight segments put a tangent kink at every
    # bead, which the steering then has to absorb in a single nucleosome.
    from scipy.interpolate import CubicSpline
    u = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(guide, axis=0), axis=1))])
    spl = CubicSpline(u, guide, bc_type="natural")
    pts, sec_id = [], []
    for k in range(len(guide) - 1):
        t = np.linspace(u[k], u[k + 1], spp, endpoint=False)
        pts.append(spl(t))
        sec_id.append(np.full(spp, k))
    pts.append(guide[-1][None, :]); sec_id.append([len(guide) - 2])
    pts = np.vstack(pts); sec_id = np.concatenate(sec_id)

    T, N1, N2 = fb.parallel_frames(pts)
    ds = np.concatenate([[0.0], np.linalg.norm(np.diff(pts, axis=0), axis=1)])

    def arc_of(k, inv_b, theta0):
        m = sec_id[:-1] == k
        idx = np.where(m)[0]
        if not len(idx):
            return 0.0, theta0, None
        lo, hi = idx[0], idx[-1] + 2
        th = theta0 + np.cumsum(np.concatenate([[0.0], ds[lo + 1:hi]])) * inv_b
        q = pts[lo:hi] + coil_r * (np.cos(th)[:, None] * N1[lo:hi]
                                   + np.sin(th)[:, None] * N2[lo:hi])
        return float(np.linalg.norm(np.diff(q, axis=0), axis=1).sum()), th[-1], (lo, hi, th)

    coil = np.zeros_like(pts)
    theta = 0.0
    filled = np.zeros(len(pts), bool)
    for k in range(len(guide) - 1):
        want = need[k]
        lo_b, hi_b = 0.0, 0.5           # inverse pitch parameter 1/b
        for _ in range(60):
            mid = 0.5 * (lo_b + hi_b)
            a, _, _ = arc_of(k, mid, theta)
            if a < want:
                lo_b = mid
            else:
                hi_b = mid
        inv_b = 0.5 * (lo_b + hi_b)
        a, theta_end, got = arc_of(k, inv_b, theta)
        if got is None:
            continue
        lo, hi, th = got
        coil[lo:hi] = pts[lo:hi] + coil_r * (np.cos(th)[:, None] * N1[lo:hi]
                                             + np.sin(th)[:, None] * N2[lo:hi])
        filled[lo:hi] = True
        theta = theta_end
    if not filled.all():
        coil[~filled] = pts[~filled]
    return coil, sec_id


def screw_of_step(st):
    """Decompose the b2 step (R, t_hat*L) into a screw about the fibre axis.

    Returns, in the nucleosome's own (body) frame: axis direction a, twist phi
    per step, the axis point for unit step length c1 (the axis offset scales
    linearly with L because R is fixed), and B, whose columns are the axis
    frame (x from axis toward the nucleosome, z along the axis)."""
    R, th = st["R"], st["t_hat"]
    w, V = np.linalg.eig(R)
    a = np.real(V[:, np.argmin(np.abs(w - 1.0))]); a /= np.linalg.norm(a)
    if np.dot(th, a) < 0:
        a = -a                                    # axis points the way the fibre advances
    t_perp = th - np.dot(th, a) * a
    c1 = np.linalg.lstsq(np.eye(3) - R, t_perp, rcond=None)[0]
    c1 -= np.dot(c1, a) * a                       # axis point nearest the origin, L = 1
    x = -c1 / np.linalg.norm(c1)                  # axis -> nucleosome
    y = np.cross(a, x)
    B = np.column_stack([x, y, a])
    # twist about a, measured in the axis frame
    Rax = B.T @ R @ B
    phi = np.arctan2(Rax[1, 0], Rax[0, 0])
    return a, phi, c1, B


def axial_advance(Ls, st):
    """Advance along the fibre axis for each step, accounting for the fibre
    radius changing with repeat length (see direct_frames). Replaces the
    uniform-screw estimate L * 0.302, which undercounts by ~5% here."""
    a, phi, c1, B = screw_of_step(st)
    rho1 = np.linalg.norm(c1)
    Lloc = np.concatenate([[Ls[0]], 0.5 * (Ls[:-1] + Ls[1:]), [Ls[-1]]])
    r = rho1 * Lloc
    lat2 = r[:-1] ** 2 + r[1:] ** 2 - 2 * r[:-1] * r[1:] * np.cos(phi)
    return np.sqrt(np.maximum(Ls ** 2 - lat2, (0.05 * Ls) ** 2))


def direct_frames(path, adv, Ls, st):
    """Place every nucleosome directly on the path: no feedback, no drift.

    Nucleosome i sits at arclength s_i = sum(adv[:i]) along the path, rotated
    by i*phi about the path's local tangent (rotation-minimising frame) and
    offset from it by the fibre radius for its local step length. On a straight
    path this reproduces the b2 step exactly; on a curved one, the difference
    from the b2 step is precisely the bend each linker has to absorb."""
    a, phi, c1, B = screw_of_step(st)
    rho1 = np.linalg.norm(c1)
    T, N1, N2 = fb.parallel_frames(path)
    seg = np.linalg.norm(np.diff(path, axis=0), axis=1)
    sp = np.concatenate([[0.0], np.cumsum(seg)])
    n = len(Ls) + 1
    Lloc = np.concatenate([[Ls[0]], 0.5 * (Ls[:-1] + Ls[1:]), [Ls[-1]]])
    # Repeat length varies, so the fibre radius (rho ~ L) varies nucleosome to
    # nucleosome. Choose each step's advance along the axis so the centre-to-
    # centre distance is exactly L_i: the lateral part is fixed by the two radii
    # and the twist, the axial part takes up the rest.
    r = rho1 * Lloc
    d = axial_advance(Ls, st)
    s = np.concatenate([[0.0], np.cumsum(d)])
    print(f"  direct: axis the chain needs {s[-1]/1e4:.2f} um, path provides "
          f"{sp[-1]/1e4:.2f} um (scale {sp[-1]/s[-1]:.3f})")
    s *= sp[-1] / s[-1]                            # map onto the path exactly
    O = np.zeros((n, 3)); Rf = np.zeros((n, 3, 3))
    for i in range(n):
        j = int(np.clip(np.searchsorted(sp, s[i]) - 1, 0, len(path) - 2))
        f = (s[i] - sp[j]) / max(seg[j], 1e-9)
        p = path[j] * (1 - f) + path[j + 1] * f
        t = T[j] * (1 - f) + T[j + 1] * f; t /= np.linalg.norm(t)
        u = N1[j] * (1 - f) + N1[j + 1] * f; u -= np.dot(u, t) * t; u /= np.linalg.norm(u)
        Q = np.column_stack([u, np.cross(t, u), t])
        th = i * phi
        Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        Rf[i] = Q @ Rz @ B.T
        O[i] = p + Q @ Rz @ np.array([r[i], 0.0, 0.0])
    return O, Rf, rho1 * np.median(Ls)


def rot_between(a, b, damp=1.0):
    """Rotation taking unit vector a onto unit vector b, optionally damped."""
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = float(np.dot(a, b))
    if s < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    ang = np.arctan2(s, c) * damp
    return cl._rodrigues(v / s * ang)


def screw_walk(n, adv, Ls, st, point_at, tangent_at, s_end, look_n=4.0, damp=0.7, lref_win=1,
               ext_of=None, max_ext=0.95, s_target=None, min_clear=120.0):
    """Exact b2 steps, steered from the EXACT fibre axis.

    Every step is the validated b2 transform, so local geometry is identical to
    the 15-nucleosome segment. Before each step the chain's own axis point and
    direction are computed from the screw decomposition (no averaging, no lag),
    and a small rotation about the current nucleosome aims that axis at a point
    look_n advances ahead on the path. The rotation angle is what the linker
    must absorb."""
    a, phi, c1, B = screw_of_step(st)
    if s_target is not None:
        s = np.asarray(s_target, float)   # genome-anchored: nucleosome i aims at ITS genomic spot
    else:
        s = np.concatenate([[0.0], np.cumsum(adv)])
        s *= s_end / s[-1]          # progress targets along the path, not past its end
    look = look_n * float(np.mean(adv))
    # Axis reference uses a running mean step length. The instantaneous axis
    # jumps sideways by rho*dL/L whenever the repeat length changes; that is
    # real, but chasing it turns every NRL change into a steering spike.
    w = int(lref_win) | 1
    Lref = np.convolve(np.pad(Ls, w // 2, mode="edge"), np.ones(w) / w, mode="valid")
    R = rot_between(a, tangent_at(0.0))
    o = point_at(0.0) - R @ (c1 * Lref[0])       # put the axis ON the path
    O = np.zeros((n, 3)); Rf = np.zeros((n, 3, 3)); corr = np.zeros(n - 1)
    guarded = np.zeros(n - 1, bool); tight = np.zeros(n - 1, bool)
    O[0], Rf[0] = o, R
    for i in range(n - 1):
        p_ax = o + R @ (c1 * Lref[i])            # axis point (smoothed step length)
        sl = s[i] + look
        tgt = point_at(min(sl, s_end))
        if sl > s_end:                            # past the end: continue straight
            tgt = tgt + tangent_at(s_end) * (sl - s_end)
        want = tgt - p_ax
        nw = np.linalg.norm(want)
        want = want / nw if nw > 1e-9 else tangent_at(s[i])
        # Guard: the correction changes linker i's geometry. If it would
        # stretch the linker past max_ext of its B-DNA contour, scale it back;
        # the uncorrected b2 step is always reachable (ext ~0.90 by design).
        # Guards, in order: the linker must stay reachable (ext <= max_ext) and the
        # new nucleosome must not collide with any already placed (>2 back).
        # Try the planned correction, then softer, none, and stronger; take the
        # first that passes both, else the reachable one with most clearance.
        best = None
        for scale in (1.0, 0.6, 0.3, 0.0, 1.6, 2.2):
            Rc_ = rot_between(R @ a, want, min(1.0, damp * scale))
            R1_ = Rc_ @ R
            o1_ = o + R1_ @ (st["t_hat"] * Ls[i])
            R1n_ = R1_ @ st["R"]
            reach = ext_of is None or ext_of(i, (O[i], Rf[i]), (o1_, R1n_)) <= max_ext
            clr = np.inf if i < 2 else float(np.min(np.linalg.norm(O[:i - 1] - o1_, axis=1)))
            cand = (reach, clr, scale, Rc_, R1_, o1_, R1n_)
            if reach and clr >= min_clear:
                best = cand; break
            if reach and (best is None or clr > best[1]):
                best = cand
        if best is None:
            best = cand
        _, clr, scale, Rc, R1, o1, R1n = best
        guarded[i] = scale != 1.0
        tight[i] = clr < min_clear
        corr[i] = np.degrees(np.arccos(np.clip((np.trace(Rc) - 1) / 2, -1, 1)))
        R = R1; o = o1
        R = R1n
        O[i + 1], Rf[i + 1] = o, R
    print(f"  screw walker: {guarded.sum()} steps adjusted by the reach/collision guards; "
          f"{tight.sum()} could not reach {min_clear:.0f} A clearance")
    return O, Rf, corr


def along_path_offset(frames_o, s_anchor, point_at, s_end, window=3000.0, step=10.0):
    """Arclength of each nucleosome's closest path point, searched only within
    +-window of its own genomic spot (a global nearest point can belong to a
    different part of a folded path, e.g. at a tether), minus that spot."""
    dense_s = np.arange(0.0, s_end + step, step)
    dense_p = np.array([point_at(v) for v in dense_s])
    lag = np.zeros(len(frames_o))
    for i, (o, sa) in enumerate(zip(frames_o, s_anchor)):
        lo = int(max(0, (sa - window) / step)); hi = int(min(len(dense_s), (sa + window) / step + 1))
        j = lo + int(np.argmin(np.sum((dense_p[lo:hi] - o) ** 2, axis=1)))
        lag[i] = dense_s[j] - sa
    return lag


S_REF=[None]; PA=[None]; SE=[None]
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default="nucleosomes_GRCh38_resolved.bed")
    ap.add_argument("--guide", default="gap-filled.xyz")
    ap.add_argument("--gaps", default="gap_sections.txt")
    ap.add_argument("--out", default="frames.npz")
    ap.add_argument("--damp", type=float, default=0.5)
    ap.add_argument("--coil-radius", type=float, default=COIL_RADIUS_A)
    ap.add_argument("--domain", type=float, default=DOMAIN_NM,
                    help="domain diameter in nm; the guide is scaled to it")
    ap.add_argument("--smooth", type=int, default=9,
                    help="sections to smooth the meander demand over")
    ap.add_argument("--walker", choices=["screw", "direct", "pursuit"], default="screw",
                    help="screw: exact b2 steps steered from the exact axis (default); "
                         "direct: place on path; pursuit: original averaged-hub steering")
    ap.add_argument("--calib", type=int, default=3,
                    help="screw walker + anchoring: rounds of re-coiling to remove along-path offset")
    ap.add_argument("--no-anchor", dest="anchor", action="store_false",
                    help="disable genomic anchoring (old behaviour: nucleosomes may slip along the path)")
    ap.add_argument("--max-ext", type=float, default=0.95,
                    help="screw walker: max linker extension allowed when steering")
    ap.add_argument("--look", type=float, default=16.0,
                    help="screw walker: lookahead, in nucleosome advances")
    ap.add_argument("--lref-win", type=int, default=1,
                    help="screw walker: steps averaged for the axis reference")
    ap.add_argument("--win", type=int, default=13,
                    help="nucleosomes averaged to estimate the fibre axis")
    args = ap.parse_args()

    pos, src = load_dyads(args.bed)
    guide = load_guide(args.guide, args.gaps, args.domain)
    n = len(pos)
    nrl = np.diff(pos)
    loc_entry, loc_exit = canonical_1kx5()
    print(f"dyads {n}, guide {len(guide)} beads, "
          f"{np.linalg.norm(np.diff(guide,axis=0),axis=1).sum()/1e4:.2f} um")

    # ---- step lengths and advances --------------------------------------
    Lcache = {}

    def step_len(link_bp):
        b = int(np.clip(link_bp, 10, 200))
        if b not in Lcache:
            Lcache[b] = fa.solve_length(b, place, loc_entry, loc_exit)[0]
        return Lcache[b]

    Ls = np.array([step_len(int(x) - 147) for x in nrl])
    adv = Ls * ADV_PER_STEP
    if args.walker == "direct":
        adv = axial_advance(Ls, fa.b2_step())
        print(f"axis demand (variable-radius fibre): {adv.sum()/1e4:.2f} um")

    # ---- per-section meander --------------------------------------------
    sec_of = np.clip((pos[:-1] - ROI_START) // BIN, 0, len(guide) - 2)
    chord = np.linalg.norm(np.diff(guide, axis=0), axis=1)
    need = np.zeros(len(guide) - 1)
    for k in range(len(guide) - 1):
        need[k] = adv[sec_of == k].sum()
    meander = np.where(need > 0, need / chord, 1.0)
    meander = np.maximum(meander, 1.0001)
    # Smooth the demand before it becomes coil pitch. A fibre cannot change its
    # coil pitch abruptly, and an unsmoothed profile makes the phase RATE jump
    # at every bead: phase stays continuous but its derivative does not, so the
    # axis tangent kinks 325 times and the steering spikes (measured up to
    # 128 deg at one step). Smoothing trades a small per-section arc error,
    # reported below, for a differentiable axis.
    if args.smooth > 1:
        w = int(args.smooth) | 1
        k = np.ones(w) / w
        need_s = np.convolve(np.pad(need, w//2, mode="edge"), k, mode="valid")
        need_s *= need.sum() / need_s.sum()        # preserve total length
        print(f"  smoothing over {w} sections: per-section arc error "
              f"median {np.median(np.abs(need_s-need)/np.maximum(need,1)):.1%}, "
              f"max {np.max(np.abs(need_s-need)/np.maximum(need,1)):.1%}")
        need = need_s
    print(f"meander: global {need.sum()/chord.sum():.2f}x, "
          f"median {np.median(meander):.2f}x, "
          f"5-95% {np.percentile(meander,5):.2f}-{np.percentile(meander,95):.2f}x, "
          f"max {meander.max():.2f}x")
    tol = 0.03 if args.coil_radius == 0 else 1e-6   # no coil: path IS the axis; walker rescales
    assert (need[need > 0] / chord[need > 0] >= 1.0 - tol).all(), \
        "a section needs LESS axis than its bead gap -- infeasible"

    need_arc = need.copy()
    for cal in range(args.calib + 1 if (args.walker == 'screw' and args.anchor) else 1):
        # ---- coiled axis -----------------------------------------------------
        coil, coil_sec = coiled_axis(guide, need_arc, args.coil_radius)
        seg = np.linalg.norm(np.diff(coil, axis=0), axis=1)
        s_coil = np.concatenate([[0.0], np.cumsum(seg)])
        print(f"coiled axis: {s_coil[-1]/1e4:.2f} um "
              f"(needed {adv.sum()/1e4:.2f} um, ratio {s_coil[-1]/adv.sum():.3f})")

        def tangent_at(s):
            i = int(np.clip(np.searchsorted(s_coil, s) - 1, 0, len(coil) - 2))
            v = coil[i + 1] - coil[i]
            return v / np.linalg.norm(v)

        def point_at(s):
            i = int(np.clip(np.searchsorted(s_coil, s) - 1, 0, len(coil) - 2))
            f = (s - s_coil[i]) / max(seg[i], 1e-9)
            return coil[i] + np.clip(f, 0, 1) * (coil[i + 1] - coil[i])

        # ---- walk the chain --------------------------------------------------
        st = fa.b2_step()
        w, V = np.linalg.eig(st["R"])
        axis_local = np.real(V[:, np.argmin(np.abs(w - 1.0))])
        axis_local /= np.linalg.norm(axis_local)

        # Pursuit steering. Aligning the local fibre axis with the axis TANGENT
        # controls direction but never position, so tiny errors integrate: over
        # 3,376 steps that produced 3,300 A of drift and a chain folded back on
        # itself at 14 A centre-to-centre. Steering instead toward a point AHEAD on
        # the axis closes the loop on position.
        #
        # The comparison is made against the fibre's own axis, estimated as the
        # running mean of the last WIN nucleosome centres, not against a single
        # nucleosome -- individual centres orbit the axis at the fibre radius, so
        # chasing one of them would inject that orbit into the steering.
        if args.walker == "screw":
            def ext_of(i, fa_, fb_):
                ex = place(fa_, loc_exit); en = place(fb_, loc_entry)
                return float(np.linalg.norm(en[0] - ex[0])) / ((int(nrl[i]) - 147 + 1) * cl.B_RISE)
            # genomic anchoring: arclength at each section start, then each dyad at
            # its fractional position within its 2 kb section
            nsec = len(guide) - 1
            sec_s = np.array([s_coil[np.argmax(coil_sec == k)] for k in range(nsec)] + [s_coil[-1]])
            u = (pos.astype(float) - ROI_START) / BIN
            kk = np.clip(np.floor(u).astype(int), 0, nsec - 1)
            s_anchor = sec_s[kk] + (u - kk) * (sec_s[kk + 1] - sec_s[kk])
            S_REF[0] = s_anchor.copy(); PA[0] = point_at; SE[0] = s_coil[-1]
            if not args.anchor:
                s_anchor = None
            frames_o, frames_R, corr = screw_walk(n, adv, Ls, st, point_at, tangent_at,
                                                  s_coil[-1], look_n=args.look, damp=args.damp,
                                                  lref_win=args.lref_win, ext_of=ext_of,
                                                  max_ext=args.max_ext, s_target=s_anchor)
        elif args.walker == "direct":
            frames_o, frames_R, rho = direct_frames(coil, adv, Ls, st)
            print(f"direct placement: fibre radius (nucleosome centre to axis) {rho:.0f} A")
            # rotation each linker must absorb beyond the ideal b2 step
            corr = np.array([np.degrees(np.arccos(np.clip(
                (np.trace(st["R"].T @ frames_R[i].T @ frames_R[i + 1]) - 1) / 2, -1, 1)))
                for i in range(n - 1)])
        else:
            WIN = args.win
            LOOKAHEAD = 4.0          # in units of mean nucleosome advance

            o = point_at(0.0)
            R = rot_between(axis_local, tangent_at(0.0))
            frames_o = np.zeros((n, 3)); frames_R = np.zeros((n, 3, 3))
            frames_o[0], frames_R[0] = o, R
            s_pos = 0.0
            corr, recent = [], [o.copy()]
            look = LOOKAHEAD * float(np.mean(adv))
            for i in range(n - 1):
                hub = np.mean(recent[-WIN:], axis=0)
                target = point_at(min(s_pos + look, s_coil[-1]))
                desired = target - hub
                nd = np.linalg.norm(desired)
                desired = desired / nd if nd > 1e-9 else tangent_at(s_pos)
                a_world = R @ axis_local
                corr.append(np.degrees(np.arccos(np.clip(np.dot(a_world, desired), -1, 1))))
                R = rot_between(a_world, desired, args.damp) @ R
                o = o + R @ (st["t_hat"] * Ls[i])
                R = R @ st["R"]
                s_pos += adv[i]
                recent.append(o.copy())
                frames_o[i + 1], frames_R[i + 1] = o, R
            corr = np.array(corr)


        # ---- calibration: lengthen the path where the chain ran ahead -------
        if cal < args.calib and args.walker == "screw" and args.anchor:
            lag = along_path_offset(frames_o, s_anchor, point_at, s_coil[-1])
            Lk = np.array([lag[kk == k].mean() if (kk == k).any() else np.nan for k in range(nsec)])
            okk = ~np.isnan(Lk)
            Lk = np.interp(np.arange(nsec), np.where(okk)[0], Lk[okk])
            inc = np.diff(np.concatenate([[0.0], Lk]))     # how far the chain over-ran in each section
            need_arc = np.maximum(need_arc + inc, 0.3 * need_arc)
            print(f"  calibration {cal+1}/{args.calib}: along-path offset median {np.median(lag):+.0f} A, "
                  f"max |{np.abs(lag).max():.0f}| A -> path re-coiled", flush=True)

    # ---- validation ------------------------------------------------------
    print("\nVALIDATION")
    print(f"  steering correction : median {np.median(corr):.2f} deg/nuc, "
          f"95th {np.percentile(corr,95):.2f}, max {corr.max():.2f}"
          f"   [expect ~4 deg; tail above 8 = linker load]")
    # distance from each nucleosome centre to the NEAREST point of the path:
    # should sit near fibre radius (~113 A) + coil radius, and never grow
    from scipy.spatial import cKDTree as _KD
    dense = np.vstack([coil[:-1] + f * (coil[1:] - coil[:-1]) for f in np.linspace(0, 1, 6)])
    drift = _KD(dense).query(frames_o)[0]
    if args.walker == "screw" and args.anchor:
        reg = np.linalg.norm(frames_o - np.array([point_at(v) for v in s_anchor]), axis=1)
        print(f"  genomic registration: nucleosome to its own genomic spot on the path "
              f"median {np.median(reg):.0f} A, 95th {np.percentile(reg,95):.0f}, max {reg.max():.0f} A"
              f"   [~113 A = on target]")
        lag = along_path_offset(frames_o, s_anchor, point_at, s_coil[-1])
        print("  along-path offset (+ ahead / - behind its genomic spot), by position along chain:")
        print("    " + "  ".join(f"{int(q*100)}%:{lag[int(q*(n-1))]:+.0f}" for q in np.linspace(0, 1, 9)) + " A")
    print(f"  distance from path  : median {np.median(drift):.0f} A, "
          f"95th {np.percentile(drift,95):.0f}, max {drift.max():.0f} A   "
          f"[expect ~113 A fibre radius + {args.coil_radius:.0f} A coil]")
    from scipy.spatial import cKDTree
    d, j = cKDTree(frames_o).query(frames_o, k=2)
    seq = np.abs(np.arange(n) - j[:, 1])
    far = seq > 2
    print(f"  min centre-centre   : {d[:,1].min():.0f} A overall; "
          f"{d[far,1].min() if far.any() else float('nan'):.0f} A between "
          f"nucleosomes >2 apart in sequence   [<110 A is a collision]")
    ext, aim = [], []
    for i in range(n - 1):                        # every linker, not a sample
        ex = place((frames_o[i], frames_R[i]), loc_exit)
        en = place((frames_o[i + 1], frames_R[i + 1]), loc_entry)
        sp = float(np.linalg.norm(en[0] - ex[0]))
        ext.append(sp / ((int(nrl[i]) - 147 + 1) * cl.B_RISE))
        dd = (en[0] - ex[0]) / sp
        aim.append(max(np.degrees(np.arccos(np.clip(np.dot(ex[1][:, 2], dd), -1, 1))),
                       np.degrees(np.arccos(np.clip(np.dot(en[1][:, 2], dd), -1, 1)))))
    ext = np.array(ext); aim = np.array(aim)
    print(f"  extension           : median {np.median(ext):.2f}, "
          f"over 1.0 (unreachable): {int((ext>1).sum())}/{len(ext)}")
    print(f"  aim                 : median {np.median(aim):.1f} deg, "
          f"over 28: {int((aim>28).sum())}/{len(aim)}   [b2 reference 12-28]")

    np.savez_compressed(args.out, origins=frames_o, rotations=frames_R,
                        dyad_pos=pos, source=src, nrl=nrl, step_len=Ls,
                        advance=adv, meander=meander,
                        coil_radius=args.coil_radius, damp=args.damp)
    print(f"\nwrote {args.out}  ({n} frames)")
    lagA = along_path_offset(frames_o, S_REF[0], PA[0], SE[0], window=3000.0)
    lagW = along_path_offset(frames_o, S_REF[0], PA[0], SE[0], window=60000.0, step=20.0)
    reg = np.linalg.norm(frames_o - np.array([PA[0](v) for v in S_REF[0]]), axis=1)
    np.savez(args.out.replace('.npz','_s1fig.npz'), corr=corr, lag_local=lagA, lag_wide=lagW, reg=reg, s_ref=S_REF[0], origins=frames_o, pos=pos, src=src)
    print('S1 Fig extras saved')


if __name__ == "__main__":
    main()
