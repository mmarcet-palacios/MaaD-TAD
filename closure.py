"""
Path-following with closure: build N base pairs that start at one frame and
ARRIVE at another with correct position, orientation and twist phase.

Marching forward from a nucleosome's exit is an initial-value problem and is
trivial. Landing on the next nucleosome's entry frame is a boundary-value
problem in SE(3) -- six constraints, of which the twist phase is the one that
bites, because over N steps you accumulate ~34.3*N degrees of twist and the
arrival register has to match.

Formulation (an elastic-rod ground state, not a simulation):

  * DNA bends and twists readily but does not stretch. So rise, slide and
    shift are held at their relaxed B-DNA values and only roll, tilt and
    twist vary. This is the standard soft-mode/stiff-mode split.

  * The three soft parameters are expressed along the chain as a short sum of
    smooth cosine modes rather than 3N independent values. That reduces the
    problem from ~600 unknowns to ~20, and more importantly it makes the
    solution a SMOOTH bend -- an elastic rod under end constraints adopts a
    smooth shape, not a chain of local kinks. Kinky closure is exactly the
    kind of artifact we are trying to stop producing.

  * Closure is enforced as a heavily weighted residual (position + rotation
    vector to the target frame); elastic strain is the quantity minimised.
    Twist-phase mismatch is therefore absorbed as distributed over/undertwist,
    which is the physical linking-number degree of freedom, rather than being
    dumped into one discontinuity.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

import bpstep as bs

# relaxed B-DNA in solution (~10.5 bp/turn)
B_TWIST = 34.3
B_RISE = 3.34
B_SLIDE = 0.0
B_SHIFT = 0.0

# relative stiffnesses for the soft modes (bend vs twist are comparable in DNA)
K_ROLL = 1.0
K_TILT = 1.0
K_TWIST = 1.2

CLOSURE_W = 3000.0     # large -> closure satisfied to numerical precision


def mode_basis(n_steps, n_modes):
    """Smooth cosine modes over the chain. Mode 0 is a uniform offset (a
    constant bend), higher modes add progressively finer variation."""
    i = np.arange(n_steps) / max(n_steps - 1, 1)
    return np.stack([np.cos(m * np.pi * i) for m in range(n_modes)], axis=1)


def helical_phase(n_steps):
    """Accumulated twist before each step. Roll and tilt are expressed in the
    LOCAL frame, which rotates by ~34.3 deg per step, so a constant roll does
    not bend the duplex coherently -- it traces a helical wobble that averages
    out. Coherent bending in a fixed lab direction requires roll/tilt to be
    modulated at the helical frequency. This is exactly why real nucleosomal
    DNA shows roll oscillating with a ~10.4 bp period, which is what the 1KX5
    extraction measured."""
    return np.radians(B_TWIST) * np.arange(n_steps)


def params_from_coeffs(coeffs, basis_rt, basis_tw, n_steps):
    n_rt = basis_rt.shape[1]
    n_tw = basis_tw.shape[1]
    a = coeffs[:n_rt]
    b = coeffs[n_rt:2 * n_rt]
    c = coeffs[2 * n_rt:2 * n_rt + n_tw]

    # smooth bend components in the LAB frame, then rotated into each step's
    # local frame by the accumulated helical phase
    u = basis_rt @ a
    v = basis_rt @ b
    phi = helical_phase(n_steps)
    cos_p, sin_p = np.cos(phi), np.sin(phi)

    # Rotation of the lab-frame bend vector BY +phi into the local frame.
    # The opposite sign (rotating by -phi) looks equally plausible and is
    # wrong: measured over 60 steps, a constant 1 deg bend deflects the duplex
    # 94.5 A with this convention and only 3.0 A with the other, because the
    # wrong sign counter-rotates against the helix and the bend averages out.
    roll = u * cos_p - v * sin_p
    tilt = u * sin_p + v * cos_p
    twist = B_TWIST + basis_tw @ c

    p = np.zeros((n_steps, 6))
    p[:, 0] = B_SHIFT
    p[:, 1] = B_SLIDE
    p[:, 2] = B_RISE
    p[:, 3] = tilt
    p[:, 4] = roll
    p[:, 5] = twist
    return p


def _rodrigues(v):
    """Fast rotation matrix from a rotation vector. scipy's Rotation objects
    are correct but allocate per call, and the chain composition inside the
    optimiser calls this tens of thousands of times."""
    theta = np.sqrt(v @ v)
    if theta < 1e-12:
        return np.eye(3)
    k = v / theta
    K = np.array([[0.0, -k[2], k[1]],
                  [k[2], 0.0, -k[0]],
                  [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def _rotvec(R):
    """Fast rotation vector from a rotation matrix."""
    c = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    theta = np.arccos(c)
    if theta < 1e-9:
        return np.zeros(3)
    if abs(np.pi - theta) < 1e-6:
        # near-pi: recover the axis from the symmetric part
        A = (R + np.eye(3)) / 2.0
        axis = np.sqrt(np.clip(np.diag(A), 0.0, None))
        if R[2, 1] - R[1, 2] < 0: axis[0] = -axis[0]
        if R[0, 2] - R[2, 0] < 0: axis[1] = -axis[1]
        if R[1, 0] - R[0, 1] < 0: axis[2] = -axis[2]
        return theta * axis / np.linalg.norm(axis)
    s = 2.0 * np.sin(theta)
    return (theta / s) * np.array([R[2, 1] - R[1, 2],
                                   R[0, 2] - R[2, 0],
                                   R[1, 0] - R[0, 1]])


def compose_chain(frame0, params):
    o, R = frame0
    for p in params:
        theta = np.radians(p[3:6])
        R_rel = _rodrigues(theta)
        R_half = _rodrigues(theta / 2.0)
        o = o + R @ (R_half @ p[0:3])
        R = R @ R_rel
    return o, R


def frame_error(frame_a, frame_b):
    """6-vector: position difference, then rotation vector from a to b."""
    o_a, R_a = frame_a
    o_b, R_b = frame_b
    dpos = o_b - o_a
    drot = _rotvec(R_a.T @ R_b)
    return np.concatenate([dpos, np.degrees(drot)])


# ---------------------------------------------------------------------------
# Obstacle awareness
#
# The closure solver used to see only the two endpoint frames. Anything the
# path crossed on the way was invisible to it, which is how a slack linker came
# to loop straight over the next nucleosome: at NRL 190 on b2's zigzag scale,
# 147 A of contour had to fit a 52.6 A gap, and the surplus went wherever the
# minimum-strain solution put it -- measured, bp 22 of 43 passed 41.2 A from
# the next core's centre, inside its 41.8 A DNA gyre, and 4.3 A from its
# histone atoms. The audit passed it because endpoint span and per-step bend
# were both fine. A big lazy loop is geometrically comfortable; it is just in
# the wrong place.
#
# A nucleosome is modelled as a solid DISC rather than a sphere, because the
# difference matters: the core is ~104 A across and ~55 A thick, so a sphere
# large enough to cover the rim needlessly blocks the faces, which is exactly
# where linkers legitimately run in a zigzag.
# ---------------------------------------------------------------------------

NUC_DISC_R = 52.0      # core radius including the wrapped DNA gyre
NUC_DISC_H = 30.0      # half-thickness including the wrapped DNA
DUPLEX_R = 10.0        # linker duplex radius -- clearance is measured to this


def trace_path(frame_start, params):
    """Base-pair origins along a solved path, including both endpoints."""
    o = np.asarray(frame_start[0], float).copy()
    R = np.asarray(frame_start[1], float).copy()
    pts = [o.copy()]
    for q in params:
        th = np.radians(q[3:6])
        o = o + R @ (_rodrigues(th / 2.0) @ q[0:3])
        R = R @ _rodrigues(th)
        pts.append(o.copy())
    return np.array(pts)


def disc_distance(points, centre, R_nuc, disc_r=NUC_DISC_R, disc_h=NUC_DISC_H):
    """Signed distance from each point to a nucleosome disc.

    Positive outside, negative inside. R_nuc is the nucleosome frame whose
    third column is the disc normal.
    """
    v = np.asarray(points, float) - np.asarray(centre, float)
    loc = v @ np.asarray(R_nuc, float)          # into the nucleosome frame
    r = np.hypot(loc[:, 0], loc[:, 1])
    z = np.abs(loc[:, 2])
    dr = r - disc_r
    dz = z - disc_h
    outside = np.hypot(np.maximum(dr, 0.0), np.maximum(dz, 0.0))
    inside = np.maximum(dr, dz)                 # negative when within both
    return np.where((dr > 0) | (dz > 0), outside, inside)


def path_clearance(pts, obstacles, skip_ends=2):
    """Smallest clearance between a path and a set of nucleosome discs.

    obstacles entries are (centre, R) or (centre, R, skip), where skip is a
    per-obstacle number of base pairs to ignore at EACH end of the path. The
    donor and acceptor cores need a larger skip than the rest: the linker
    emerges from the donor's surface, and the exit base pair sits only 4.5 A
    outside the disc by construction, so counting those base pairs would
    register a permanent, unavoidable collision and mask every real one.
    """
    if obstacles is None or len(obstacles) == 0:
        return np.inf
    n = len(pts)
    worst = np.inf
    for ob in obstacles:
        if isinstance(ob[0], str) and ob[0] == "pts":
            # point-cloud obstacle ("pts", cKDTree, need): atoms the duplex AXIS must
            # stay `need` A from (histone tails, a tether partner's built backbone).
            # Reported on the disc scale, so clearance_req=DUPLEX_R still means "clear".
            _, tree, need = ob
            q = pts[2:n - 2] if n > 4 else pts
            if len(q):
                worst = min(worst, float(tree.query(q)[0].min()) - need + DUPLEX_R)
            continue
        c, R = ob[0], ob[1]
        sk = ob[2] if len(ob) > 2 else skip_ends
        q = pts[sk:n - sk] if n > 2 * sk else None
        if q is None or len(q) == 0:
            continue
        worst = min(worst, float(disc_distance(q, c, R).min()))
    return worst


def nucleosome_obstacles(frames, donor=None, acceptor=None,
                         skip_other=1, skip_flank=6):
    """Obstacle list for a linker running from frames[donor] to frames[acceptor]."""
    out = []
    for i, (c, R) in enumerate(frames):
        sk = skip_flank if i in (donor, acceptor) else skip_other
        out.append((c, R, sk))
    return out


def solve_path(frame_start, frame_end, n_steps, n_modes_bend=8, n_modes_twist=4,
               verbose=False):
    """Solve step parameters taking frame_start to frame_end in n_steps."""
    basis_rt = mode_basis(n_steps, n_modes_bend)
    basis_tw = mode_basis(n_steps, n_modes_twist)
    n_coeff = 2 * n_modes_bend + n_modes_twist

    def closure_of(coeffs):
        p = params_from_coeffs(coeffs, basis_rt, basis_tw, n_steps)
        return frame_error(compose_chain_fast(frame_start, p), frame_end)

    # Stiffness weighting for the minimum-norm step: bending and twisting are
    # the soft modes, and a mode's cost scales with how much of the chain it
    # deforms. This makes the Gauss-Newton step prefer gentle, distributed
    # deformation over a sharp local one.
    w_diag = np.concatenate([
        np.full(n_modes_bend, K_ROLL),
        np.full(n_modes_bend, K_TILT),
        np.full(n_modes_twist, K_TWIST),
    ])
    Winv = np.diag(1.0 / w_diag)

    # Constrained Gauss-Newton. The closure error is very nearly linear in the
    # coefficients (a constant 1 deg bend deflects 94.5 A over 60 steps, and
    # deflection scales linearly with amplitude to within a percent), so
    # Newton steps converge in a handful of iterations. At each step take the
    # stiffness-weighted MINIMUM-NORM correction that satisfies the
    # linearised closure -- starting from zero strain, this lands on
    # essentially the minimum-strain solution rather than an arbitrary one.
    # Seed with a small bend. The perfectly straight configuration is a
    # near-singular point of this parametrisation: a straight chain does not
    # move when you twist it, so the twist columns of the Jacobian carry no
    # positional sensitivity and the condition number is ~1e5. Displacing by
    # a fraction of a degree of bend drops it to ~60, and Newton then
    # converges in a handful of steps. Measured, not guessed.
    x = np.zeros(n_coeff)
    x[0] = 0.3
    x[n_modes_bend] = 0.2
    err = closure_of(x)
    for it in range(60):
        if np.linalg.norm(err[:3]) < 1e-9 and np.linalg.norm(err[3:]) < 1e-9:
            break
        J = np.zeros((6, n_coeff))
        h = 1e-5
        for k in range(n_coeff):
            xp = x.copy(); xp[k] += h
            J[:, k] = (closure_of(xp) - err) / h
        JW = J @ Winv @ J.T
        # Levenberg damping scaled to the problem, so a near-singular
        # direction produces a bounded step instead of an explosion.
        damp = 1e-9 * np.trace(JW) / 6.0
        try:
            lam = np.linalg.solve(JW + damp * np.eye(6), -err)
        except np.linalg.LinAlgError:
            lam = np.linalg.lstsq(JW, -err, rcond=None)[0]
        step = Winv @ J.T @ lam
        # damp if the linearisation over-reaches
        for scale in (1.0, 0.5, 0.25, 0.1):
            trial = x + scale * step
            e_new = closure_of(trial)
            if np.linalg.norm(e_new) < np.linalg.norm(err) or scale == 0.1:
                x, err = trial, e_new
                break

    p = params_from_coeffs(x, basis_rt, basis_tw, n_steps)
    err = closure_of(x)

    if verbose:
        print(f"   closure: |dpos| = {np.linalg.norm(err[:3]):.2e} A, "
              f"|drot| = {np.linalg.norm(err[3:]):.2e} deg  ({it+1} iters)")
    return p, err, x


def strain_report(params, label=""):
    roll, tilt, twist = params[:, 4], params[:, 3], params[:, 5]
    print(f"  {label}")
    print(f"    roll  : {roll.mean():+6.2f} +/- {roll.std():5.2f} deg"
          f"   |max| {np.abs(roll).max():5.2f}")
    print(f"    tilt  : {tilt.mean():+6.2f} +/- {tilt.std():5.2f} deg"
          f"   |max| {np.abs(tilt).max():5.2f}")
    print(f"    twist : {twist.mean():6.2f} +/- {twist.std():5.2f} deg"
          f"   (relaxed {B_TWIST})  -> {360/twist.mean():.2f} bp/turn")
    bend_per_step = np.sqrt(roll**2 + tilt**2)
    print(f"    total bend/step : {bend_per_step.mean():.2f} deg  "
          f"|max| {bend_per_step.max():.2f} deg")
    return bend_per_step


def straight_endpoint(frame_start, n_steps):
    """Where relaxed, unbent B-DNA would land -- the no-strain reference."""
    p = np.tile([B_SHIFT, B_SLIDE, B_RISE, 0.0, 0.0, B_TWIST], (n_steps, 1))
    return compose_chain(frame_start, p)


# Elastic limits for real duplex DNA. Nucleosomal DNA -- about the most
# sharply bent DNA in biology -- has roll spanning -20 to +23 deg per step
# (measured in 1KX5). Beyond roughly that, DNA does not bend elastically, it
# kinks or melts. The solver will happily satisfy closure at 265 deg/step if
# asked, so geometry that closes must still be checked for whether it is DNA.
ELASTIC_BEND_PER_STEP = 15.0      # comfortable
MAX_BEND_PER_STEP = 25.0          # nucleosome-grade, the practical ceiling


def feasibility(frame_start, frame_end, n_steps, params=None):
    """Is this linker geometry physically buildable? Two separate ways to fail:
    not enough DNA to reach, or enough DNA but an unphysical bend required."""
    contour = n_steps * B_RISE
    reach = float(np.linalg.norm(frame_end[0] - frame_start[0]))
    out = {
        "contour_A": contour,
        "reach_A": reach,
        "reach_fraction": reach / contour,
        "reachable": reach < contour,
        "reach_deficit_A": max(0.0, reach - contour),
    }
    if params is not None:
        bend = np.sqrt(params[:, 4] ** 2 + params[:, 3] ** 2)
        out["bend_mean_deg"] = float(bend.mean())
        out["bend_max_deg"] = float(bend.max())
        out["twist_mean_deg"] = float(params[:, 5].mean())
        out["elastic"] = bool(bend.max() <= ELASTIC_BEND_PER_STEP)
        out["physical"] = bool(bend.max() <= MAX_BEND_PER_STEP)
        if not out["physical"]:
            out["verdict"] = "CLOSES BUT IS NOT DNA -- bend exceeds what duplex DNA can do"
        elif not out["elastic"]:
            out["verdict"] = "strained, nucleosome-grade bending"
        else:
            out["verdict"] = "comfortable"
    if not out["reachable"]:
        out["verdict"] = f"UNREACHABLE -- short by {out['reach_deficit_A']:.1f} A"
    return out


def compose_chain_fast(frame0, params):
    """Same as compose_chain but with the per-step Rodrigues evaluations
    batched. The sequential accumulation is unavoidable, but building all the
    rotation matrices in one vectorised pass removes most of the cost, which
    matters when auditing ~1750 linkers of a few hundred steps each."""
    th = np.radians(params[:, 3:6])
    n = len(th)
    ang = np.sqrt((th * th).sum(axis=1))
    safe = np.where(ang > 1e-12, ang, 1.0)
    k = th / safe[:, None]

    K = np.zeros((n, 3, 3))
    K[:, 0, 1] = -k[:, 2]; K[:, 0, 2] = k[:, 1]
    K[:, 1, 0] = k[:, 2];  K[:, 1, 2] = -k[:, 0]
    K[:, 2, 0] = -k[:, 1]; K[:, 2, 1] = k[:, 0]
    KK = K @ K
    I = np.broadcast_to(np.eye(3), (n, 3, 3))

    s, c = np.sin(ang), 1.0 - np.cos(ang)
    R_rel = I + s[:, None, None] * K + c[:, None, None] * KK
    s2, c2 = np.sin(ang / 2.0), 1.0 - np.cos(ang / 2.0)
    R_half = I + s2[:, None, None] * K + c2[:, None, None] * KK

    zero = ang <= 1e-12
    if zero.any():
        R_rel[zero] = np.eye(3)
        R_half[zero] = np.eye(3)

    trans = np.einsum('nij,nj->ni', R_half, params[:, 0:3])

    o, R = frame0
    o = np.asarray(o, dtype=float).copy()
    R = np.asarray(R, dtype=float).copy()
    for i in range(n):
        o = o + R @ trans[i]
        R = R @ R_rel[i]
    return o, R


def solve_path_multistart(frame_start, frame_end, n_steps, n_modes_bend=8,
                          n_modes_twist=4, n_seeds=12, good_enough=None,
                          obstacles=None, clearance_req=DUPLEX_R):
    """Multi-start wrapper around solve_path, optionally obstacle-aware.

    The single-seed Gauss-Newton reaches closure but can land in a badly
    strained basin, and it does so exactly where it matters: short SLACK
    linkers, where the straight start is far from any valid solution. Measured
    on a 33 bp linker spanning 59.5 A, a circular arc needs 6.3 deg/step for
    position alone, while the single-seed solver returned 158-308 deg/step at
    every mode count tried. The bad basin was the initialisation, not the
    basis.

    Seeding from arcs of assorted curvature and bend direction and keeping the
    least-strained solution that closes fixes it.

    obstacles: list of (centre, R) nucleosome frames the path must clear by
    clearance_req Angstrom. Without this the solver optimises strain alone and
    will happily route a slack linker through a neighbouring core -- which is
    precisely the artifact this argument exists to stop. Candidates are ranked
    by PENETRATION FIRST and strain second, because a strained path is a bad
    model whereas an interpenetrating one is not a structure at all.
    """
    contour = n_steps * B_RISE
    span = float(np.linalg.norm(np.asarray(frame_end[0]) - np.asarray(frame_start[0])))

    # curvature a circular arc would need for the position constraint alone
    turn = 0.0
    if 0 < span < contour:
        from scipy.optimize import brentq
        try:
            turn = brentq(lambda th: 2 * np.sin(th / 2) / th - span / contour,
                          1e-6, 2 * np.pi - 1e-6)
        except ValueError:
            turn = 0.0
    base = np.degrees(turn) / max(n_steps, 1)

    # Seeding. Mode 0 alone is a single circular ARC, and that is the wrong
    # shape here: the endpoints need a large displacement but almost no net
    # change of helical-axis direction (measured across this fibre, the axis
    # turns only 9.7 deg end to end). The path that does that is an S-CURVE,
    # which is mode 1. Seeding arcs only left the solver hunting for an S from
    # an arc start, and it landed wherever it happened to: bend came out
    # non-monotonic in linker length -- 21.6, 31.4, 41.0, 35.8, 37.0, 33.3 deg
    # per step for 15 to 20 bp across the same gap -- which is the signature of
    # scattered basins, not of a geometric requirement.
    seeds = []
    amps = (0.0, 0.5, 1.0, 1.5, 2.0) if obstacles else (0.0, 0.5, 1.0, 1.5)
    n_phase = max(n_seeds // 2, 8) if obstacles else max(n_seeds // 4, 1)
    for amp in amps:
        for k in range(n_phase):
            phi = 2 * np.pi * k / n_phase
            for shape in ("arc", "s", "both"):
                v = np.zeros(2 * n_modes_bend + n_modes_twist)
                a = base * amp
                if shape in ("arc", "both"):
                    v[0] = a * np.cos(phi)
                    v[n_modes_bend] = a * np.sin(phi)
                if shape in ("s", "both") and n_modes_bend > 1:
                    b = a * (2.0 if shape == "s" else 1.2)
                    v[1] = b * np.cos(phi)
                    v[n_modes_bend + 1] = b * np.sin(phi)
                seeds.append(v)

    # Ranking. Penetration used to dominate the score outright, which threw
    # away good paths for marginal clearance gains -- on a 14 bp linker it
    # swapped a 30.4 deg/step path for a 52.4 deg/step one to move the
    # clearance from 9 A to 10 A. Clearance is scored smoothly instead, so it
    # only outweighs strain when the shortfall is real; an actual overlap
    # (negative clearance) is still effectively disqualifying.
    def score_of(bend, clr):
        if clr < 0.0:
            return 1e4 + abs(clr) * 100.0 + bend
        short = max(0.0, clearance_req - clr)
        return bend + 4.0 * short ** 2

    best = None
    for seed in seeds:
        try:
            p, err, x = _solve_from(frame_start, frame_end, n_steps,
                                    n_modes_bend, n_modes_twist, seed)
        except Exception:
            continue
        if np.linalg.norm(err[:3]) > 1e-3 or np.linalg.norm(err[3:]) > 1e-3:
            continue
        bend = float(np.sqrt(p[:, 4] ** 2 + p[:, 3] ** 2).max())
        clr = (path_clearance(trace_path(frame_start, p), obstacles)
               if obstacles else np.inf)
        sc = score_of(bend, clr)
        if best is None or sc < best[4]:
            best = (p, err, x, bend, sc, clr)
        # stop once a comfortable solution is in hand -- trying the rest of
        # the seed grid only matters when nothing good has been found yet
        if (good_enough is not None and bend <= good_enough
                and (not obstacles or clr >= clearance_req)):
            return best[0], best[1], best[2]
    if best is None:
        return solve_path(frame_start, frame_end, n_steps,
                          n_modes_bend, n_modes_twist)
    return best[0], best[1], best[2]


def path_quality(frame_start, params, obstacles=None):
    """Bend and obstacle clearance of a solved path."""
    pts = trace_path(frame_start, params)
    return {"bend_max_deg": float(np.sqrt(params[:, 4] ** 2 + params[:, 3] ** 2).max()),
            "clearance_A": path_clearance(pts, obstacles) if obstacles else np.inf,
            "trace": pts}


def _solve_from(frame_start, frame_end, n_steps, n_modes_bend, n_modes_twist, seed):
    """solve_path's Gauss-Newton, started from an explicit seed."""
    basis_rt = mode_basis(n_steps, n_modes_bend)
    basis_tw = mode_basis(n_steps, n_modes_twist)
    n_coeff = 2 * n_modes_bend + n_modes_twist

    def closure_of(c):
        p = params_from_coeffs(c, basis_rt, basis_tw, n_steps)
        return frame_error(compose_chain_fast(frame_start, p), frame_end)

    w = np.concatenate([np.full(n_modes_bend, K_ROLL),
                        np.full(n_modes_bend, K_TILT),
                        np.full(n_modes_twist, K_TWIST)])
    Winv = np.diag(1.0 / w)

    x = np.array(seed, dtype=float)
    if np.allclose(x, 0):
        x[0] = 0.3
        x[n_modes_bend] = 0.2
    err = closure_of(x)
    for _ in range(60):
        if np.linalg.norm(err[:3]) < 1e-9 and np.linalg.norm(err[3:]) < 1e-9:
            break
        J = np.zeros((6, n_coeff))
        h = 1e-5
        for k in range(n_coeff):
            xp = x.copy(); xp[k] += h
            J[:, k] = (closure_of(xp) - err) / h
        JW = J @ Winv @ J.T
        damp = 1e-9 * np.trace(JW) / 6.0
        try:
            lam = np.linalg.solve(JW + damp * np.eye(6), -err)
        except np.linalg.LinAlgError:
            lam = np.linalg.lstsq(JW, -err, rcond=None)[0]
        step = Winv @ J.T @ lam
        for s in (1.0, 0.5, 0.25, 0.1):
            trial = x + s * step
            e_new = closure_of(trial)
            if np.linalg.norm(e_new) < np.linalg.norm(err) or s == 0.1:
                x, err = trial, e_new
                break
    p = params_from_coeffs(x, basis_rt, basis_tw, n_steps)
    return p, closure_of(x), x
