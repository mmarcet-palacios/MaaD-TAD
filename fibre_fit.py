"""
Fit a repeating nucleosome step transform from a real array, correctly.

Two bugs in the first attempt, both fixed here:

  1. FRAME SIGN. The disc normal came from an SVD, whose sign is arbitrary, and
     nothing pinned it. A frame that flips between neighbours scrambles every
     step transform derived from it. Fixed by signing z from the HANDEDNESS of
     the DNA wrap: sum(r_i x r_{i+1}) over the wrapped path is an
     angular-momentum-like vector that points along the superhelical axis in a
     direction set by how the DNA actually winds, which is a property of the
     structure rather than of the numerics.

  2. ROTATION AVERAGING. Averaging rotation vectors componentwise is only valid
     for small angles. These steps are ~140 deg apart, where it is badly wrong:
     it turned a set of ~160 deg rotations into a 111 deg "mean". Fixed with the
     chordal (Frobenius) mean: sum the rotation matrices and project back onto
     SO(3) by SVD.
"""
import numpy as np


def wrap_handedness(pts, centre):
    """Direction of the superhelical wrap, from the DNA path itself."""
    r = pts - centre
    L = np.cross(r[:-1], r[1:]).sum(axis=0)
    n = np.linalg.norm(L)
    return L / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])


def nucleosome_frame(wrapped_pts, dyad_pt):
    """centre / z (disc normal, signed by wrap handedness) / y (toward dyad)."""
    centre = wrapped_pts.mean(axis=0)
    _, _, vt = np.linalg.svd(wrapped_pts - centre)
    z = vt[2]
    if np.dot(z, wrap_handedness(wrapped_pts, centre)) < 0:
        z = -z
    y = dyad_pt - centre
    y = y - np.dot(y, z) * z
    ny = np.linalg.norm(y)
    if ny < 1e-9:
        raise ValueError("degenerate dyad direction")
    y /= ny
    x = np.cross(y, z)
    return centre, np.column_stack([x, y, z])


def mean_rotation(Rs):
    """Chordal mean on SO(3): sum, then project by SVD. Valid at any angle."""
    M = np.sum(Rs, axis=0)
    U, _, Vt = np.linalg.svd(M)
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1.0, 1.0, d]) @ Vt


def step_transforms(frames):
    out = []
    for (o_a, R_a), (o_b, R_b) in zip(frames[:-1], frames[1:]):
        out.append((R_a.T @ (o_b - o_a), R_a.T @ R_b))
    return out


def rot_angle(R):
    return np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))


def fit_step(frames):
    """Mean step transform plus a regularity report."""
    steps = step_transforms(frames)
    T = np.array([t for t, _ in steps])
    Rs = [R for _, R in steps]
    t_mean = T.mean(axis=0)
    R_mean = mean_rotation(Rs)

    # residuals of each step against the fitted mean
    dt = np.linalg.norm(T - t_mean, axis=1)
    dR = np.array([rot_angle(R.T @ R_mean) for R in Rs])
    return {
        "t": t_mean, "R": R_mean,
        "t_norm": float(np.linalg.norm(t_mean)),
        "rot_deg": float(rot_angle(R_mean)),
        "step_dist": np.linalg.norm(T, axis=1),
        "resid_t": dt, "resid_R": dR,
    }


def generate(n, t, R, start=None):
    if start is None:
        start = (np.zeros(3), np.eye(3))
    out = [start]
    for _ in range(n - 1):
        o, Rc = out[-1]
        out.append((o + Rc @ t, Rc @ R))
    return out


def helix_params(frames):
    C = np.array([f[0] for f in frames])
    n = len(C)
    cen = C.mean(axis=0)
    _, _, vt = np.linalg.svd(C - cen)
    axis = vt[0]
    proj = (C - cen) @ axis
    if proj[-1] < proj[0]:
        axis, proj = -axis, -proj
    rvec = (C - cen) - np.outer(proj, axis)
    radius = np.linalg.norm(rvec, axis=1)
    e1 = rvec[0] / np.linalg.norm(rvec[0])
    e2 = np.cross(axis, e1)
    ang = np.unwrap(np.arctan2(rvec @ e2, rvec @ e1))
    dphi = np.degrees(np.diff(ang))
    rise = (proj.max() - proj.min()) / (n - 1)
    return {
        "radius_A": float(radius.mean()), "radius_sd": float(radius.std()),
        "rise_A": float(rise),
        "deg_per_nuc": float(np.median(dphi)),
        "nuc_per_turn": float(abs(360 / np.median(dphi))),
        "pitch_A": float(abs(360 / np.median(dphi)) * rise),
        "diameter_nm": float(2 * radius.mean() / 10),
    }
