"""
Adaptive zigzag: the nucleosome-to-nucleosome STEP LENGTH follows the local
linker length.

The bug this fixes
------------------
The b2-measured zigzag was measured on an array whose own NRL is ~165. Running
it at NRL 190 puts 147 A of linker contour into a 52.6 A exit->entry gap --
extension 0.36, nearly 3x more DNA than the gap needs. The closure solver then
did the only thing available to it and threw a wide smooth loop, which swung
over the NEXT nucleosome: measured, bp 22 of 43 passed 41.2 A from its centre,
inside the 41.8 A DNA gyre, and 4.3 A from its histone atoms. In a viewer that
reads exactly as "the linker aims at the centre of the next nucleosome". It is
not a targeting error -- the endpoints were right to 4.17 A. It is slack.

Why the step transform and not a cylinder
-----------------------------------------
The first attempt parameterised the fibre as a cylinder and scaled its radius
per nucleosome. That fails on heterogeneous NRL, and not for a fixable reason:
the radius at nucleosome i is SHARED by the linkers on either side of it, so a
41 bp linker and a 14 bp linker that meet at the same nucleosome pull it in
opposite directions. Measured, a joint least-squares over 8 such linkers left 5
outside any reasonable band.

Iterating b2's fitted step transform instead gives a knob per LINKER. Scaling
the step translation between nucleosome i and i+1 moves i+1 and everything
downstream rigidly, so the relative geometry of every OTHER linker is
untouched. Verified: with L = [136.4, L1, 136.4], varying L1 over 100-260 A
leaves the flanking spans at 57.40 A to all printed digits. Each linker
therefore solves independently, in one dimension, with no optimiser.

The step rotation is held at b2's measured 138.69 deg, so the zigzag register
is preserved and only the spacing breathes.

A hard lower bound falls out
----------------------------
Span is not monotonic in step length -- it bottoms out near L = 100 A at about
40.3 A and rises on both sides. The short-L branch is unusable because the step
length IS the centre-to-centre distance, so it collides. On the usable branch
the smallest reachable span is ~43 A, which needs about 13 bp of linker at 0.95
extension. So this zigzag cannot represent NRL below ~159 at all. That is a
statement about the geometry, not a solver limitation: short-linker chromatin
adopts a stacked arrangement rather than an open zigzag, and should be modelled
as one.
"""
import numpy as np

import fibre_fit as ff

TARGET_EXT = 0.90
EXT_HI = 0.95            # reachability; above this the duplex must stretch
B_RISE = 3.34

L_MIN = 110.0            # centre-to-centre floor: below this the cores collide
L_MAX = 600.0

_STEP_CACHE = {}


def contour(link_bp):
    """Contour available to a linker of link_bp interior base pairs.

    The path runs from the donor core's terminal bp to the acceptor core's
    terminal bp, so it spans link_bp + 1 steps.
    """
    return (np.asarray(link_bp) + 1) * B_RISE


def b2_step(src_frames=None):
    """b2's fitted nucleosome->nucleosome step: unit direction and rotation."""
    if "fit" not in _STEP_CACHE:
        if src_frames is None:
            from build_zigzag import b2_frames
            src_frames = b2_frames()
        fit = ff.fit_step(src_frames)
        t = np.asarray(fit["t"], float)
        _STEP_CACHE["fit"] = {"t_hat": t / np.linalg.norm(t),
                              "L0": float(np.linalg.norm(t)),
                              "R": np.asarray(fit["R"], float),
                              "rot_deg": float(fit["rot_deg"])}
    return _STEP_CACHE["fit"]


def build_chain(lengths, step=None, origin=None, R0=None):
    """Place nucleosome frames by iterating the step transform.

    lengths[i] is the centre-to-centre distance from nucleosome i to i+1.
    """
    st = step or b2_step()
    o = np.zeros(3) if origin is None else np.asarray(origin, float).copy()
    R = np.eye(3) if R0 is None else np.asarray(R0, float).copy()
    frames = [(o.copy(), R.copy())]
    for L in lengths:
        o = o + R @ (st["t_hat"] * float(L))
        R = R @ st["R"]
        frames.append((o.copy(), R.copy()))
    return frames


def span_of(L, place_fn, loc_entry, loc_exit, step=None):
    """Exit->entry span produced by a single step of length L."""
    fr = build_chain([L], step=step)
    return float(np.linalg.norm(place_fn(fr[1], loc_entry)[0]
                                - place_fn(fr[0], loc_exit)[0]))


def solve_length(link_bp, place_fn, loc_entry, loc_exit,
                 target_ext=TARGET_EXT, step=None, tol=1e-4):
    """Step length putting one linker at target_ext, on the usable branch.

    Returns (L, span, extension, note). note is "" when the target was met,
    or a short reason when it could not be.
    """
    C = float(contour(link_bp))
    want = target_ext * C

    lo, hi = L_MIN, L_MAX
    s_lo = span_of(lo, place_fn, loc_entry, loc_exit, step)
    s_hi = span_of(hi, place_fn, loc_entry, loc_exit, step)

    if want <= s_lo:
        # already at the steric floor and the span is still too big
        ext = s_lo / C
        note = "" if ext <= EXT_HI else f"UNREACHABLE (ext {ext:.2f} > {EXT_HI})"
        return lo, s_lo, ext, note
    if want >= s_hi:
        return hi, s_hi, s_hi / C, "step length clipped at L_MAX"

    # span is monotonically increasing on [L_MIN, L_MAX]; plain bisection
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        s = span_of(mid, place_fn, loc_entry, loc_exit, step)
        if abs(s - want) < tol:
            break
        if s < want:
            lo = mid
        else:
            hi = mid
    L = 0.5 * (lo + hi)
    s = span_of(L, place_fn, loc_entry, loc_exit, step)
    return L, s, s / C, ""


AIM_MAX = 28.0          # b2's own exit/entry aim spans 12-28 deg


def aim_angles(frame_i, frame_j, place_fn, loc_entry, loc_exit):
    """How well the DNA leaves and arrives ALONG its direction of travel.

    This is the constraint the extension target misses. A linker can be
    perfectly reachable and still need a violent bend if the duplex exits the
    donor pointing 43 deg away from where it has to go, because it must then
    turn away and turn back. Measured on b2's own array the exit axis sits
    12-25 deg off the travel direction and the entry axis 15-28 deg, so that
    band is the reference.
    """
    ex = place_fn(frame_i, loc_exit)
    en = place_fn(frame_j, loc_entry)
    d = en[0] - ex[0]
    s = np.linalg.norm(d)
    d = d / s
    out = np.degrees(np.arccos(np.clip(np.dot(ex[1][:, 2], d), -1, 1)))
    inn = np.degrees(np.arccos(np.clip(np.dot(en[1][:, 2], d), -1, 1)))
    return float(out), float(inn), float(s)


def supported_nrl(place_fn, loc_entry, loc_exit, target_ext=TARGET_EXT,
                  aim_max=AIM_MAX, step=None):
    """The NRL window this zigzag actually supports.

    Aim is a function of step length alone, and it is good only over a middle
    band: too short and the exit axis swings away from the travel direction
    (43 deg at the steric floor), too long and the entry axis does. Converting
    that band of step lengths into linker lengths at the working extension
    gives the supported NRL range. Measured: about 165 to 195, with b2's own
    NRL sitting near its centre -- the fibre geometry selects a repeat length.
    """
    Ls = np.linspace(L_MIN, 420.0, 300)
    ok = []
    for L in Ls:
        fr = build_chain([L], step=step)
        out, inn, s = aim_angles(fr[0], fr[1], place_fn, loc_entry, loc_exit)
        if max(out, inn) <= aim_max:
            ok.append(147 + s / (target_ext * B_RISE) - 1)
    if not ok:
        return None
    return int(np.ceil(min(ok))), int(np.floor(max(ok)))


def solve_fibre(link_bp, place_fn, loc_entry, loc_exit,
                target_ext=TARGET_EXT, step=None, verbose=False):
    """Solve every step length, then build the frames.

    Each linker is solved on its own because the step parameterisation
    decouples them exactly.
    """
    link_bp = np.asarray(link_bp, int)
    rows = [solve_length(b, place_fn, loc_entry, loc_exit, target_ext, step)
            for b in link_bp]
    L = np.array([r[0] for r in rows])
    spans = np.array([r[1] for r in rows])
    ext = np.array([r[2] for r in rows])
    notes = [r[3] for r in rows]
    frames = build_chain(L, step=step)
    aims = [aim_angles(frames[i], frames[i + 1], place_fn, loc_entry, loc_exit)
            for i in range(len(link_bp))]
    for i, (out, inn, _) in enumerate(aims):
        if max(out, inn) > AIM_MAX and not notes[i]:
            notes[i] = f"AIM {max(out, inn):.0f} deg (b2 reference 12-28)"
    if verbose:
        for i, (b, l, s, e, n) in enumerate(zip(link_bp, L, spans, ext, notes)):
            print(f"    linker {i}: {b:2d} bp  contour {contour(b):6.1f} A  "
                  f"span {s:6.1f} A  ext {e:.3f}  step {l:6.1f} A  "
                  f"aim {aims[i][0]:4.1f}/{aims[i][1]:4.1f} deg"
                  + (f"   <-- {n}" if n else ""))
    return frames, L, spans, ext, notes


def min_nrl(place_fn, loc_entry, loc_exit, step=None, ext_hi=EXT_HI):
    """Smallest NRL this zigzag can represent, from the span floor."""
    Ls = np.linspace(L_MIN, 260.0, 200)
    s = np.array([span_of(x, place_fn, loc_entry, loc_exit, step) for x in Ls])
    s_min = float(s.min())
    link = int(np.ceil(s_min / (ext_hi * B_RISE) - 1))
    return 147 + link, s_min, float(Ls[int(np.argmin(s))])
