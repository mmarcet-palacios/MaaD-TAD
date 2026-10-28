"""
Three-level chromatin construction.

  level 1  bp        -- the DNA helix (bpstep / closure)
  level 2  nucleosome -- the zigzag fibre, measured from b2
  level 3  fibre      -- the fibre axis MEANDERS, because a straight fibre at
                        biological NRL is 1.5-2x too long for the contact-map
                        spline (9.13 um). The old pipeline had no level 3, so
                        the 6.3x excess landed on the linkers, where it was
                        hopeless. Here it lands on the fibre axis, where a
                        ~1.8x coil absorbs it comfortably.

Level 2 is built in cylindrical coordinates about the (steerable) fibre axis
rather than by blindly iterating a step transform, so that steering and the
zigzag stay independent. The cylindrical parameters ARE the b2 measurements:

    radius 68.6 A, 139.1 deg per nucleosome, axial rise 45.6 A,
    dyad pointing at the axis (measured 11.9 +/- 6.5 deg),
    disc normal ~perpendicular to radial (82.8 deg) and 44.0 deg to the axis.
"""
import numpy as np

# measured in b2_onelayer
R_FIBRE = 68.6
DEG_PER_NUC = 139.1
RISE_PER_NUC = 45.6
NORMAL_VS_AXIS = 44.0


def parallel_frames(points):
    """Rotation-minimising (parallel transport) normal frame along a curve.
    A naive Frenet frame spins wildly where curvature is low, which would
    inject spurious rotation into the zigzag phase."""
    T = np.gradient(points, axis=0)
    T /= np.linalg.norm(T, axis=1)[:, None]

    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(T[0], ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    n1 = ref - np.dot(ref, T[0]) * T[0]
    n1 /= np.linalg.norm(n1)

    N1 = [n1]
    for i in range(1, len(T)):
        v = N1[-1]
        v = v - np.dot(v, T[i]) * T[i]
        nv = np.linalg.norm(v)
        N1.append(v / nv if nv > 1e-9 else N1[-1])
    N1 = np.array(N1)
    N2 = np.cross(T, N1)
    return T, N1, N2


def resample_by_arclength(points, spacing):
    """Evenly spaced samples along a polyline, plus the arc length."""
    d = np.linalg.norm(np.diff(points, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(d)])
    total = s[-1]
    n = max(int(total // spacing), 1)
    targets = np.arange(n) * spacing
    out = np.empty((n, 3))
    for k in range(3):
        out[:, k] = np.interp(targets, s, points[:, k])
    return out, total


def coil_axis(centre_points, coil_radius, coil_pitch, samples=20000):
    """Wrap a helical meander around a guide curve so the fibre axis is longer
    than the curve it follows. arc/progress = sqrt(1 + (2*pi*r/p)^2)."""
    pts, total = resample_by_arclength(centre_points, total_spacing(centre_points, samples))
    T, N1, N2 = parallel_frames(pts)
    s = np.arange(len(pts)) * (total / max(len(pts) - 1, 1))
    theta = 2 * np.pi * s / coil_pitch
    return pts + coil_radius * (np.cos(theta)[:, None] * N1
                                + np.sin(theta)[:, None] * N2)


def total_spacing(points, samples):
    d = np.linalg.norm(np.diff(points, axis=0), axis=1).sum()
    return d / samples


def coil_factor(coil_radius, coil_pitch):
    return float(np.sqrt(1.0 + (2 * np.pi * coil_radius / coil_pitch) ** 2))


def solve_coil(required_factor, coil_radius):
    """Pitch that produces the required arc/progress ratio at a given radius."""
    if required_factor <= 1.0:
        return np.inf
    return 2 * np.pi * coil_radius / np.sqrt(required_factor ** 2 - 1.0)


def build_nucleosomes(axis_points, n_nuc, radius=R_FIBRE, deg=DEG_PER_NUC,
                      rise=RISE_PER_NUC, normal_vs_axis=NORMAL_VS_AXIS,
                      phase0=0.0, flip_normal=True):
    """Place nucleosome frames around a (possibly meandering) fibre axis.

    Returns a list of (centre, R) with columns [x, y, z] where
       y = dyad, pointing IN at the fibre axis
       z = disc normal, perpendicular to y and inclined to the axis
    """
    pts, total = resample_by_arclength(axis_points, rise)
    if len(pts) < 2:
        raise ValueError("axis too short for even one step")
    n_nuc = min(n_nuc, len(pts))
    T, N1, N2 = parallel_frames(pts)

    beta = np.radians(normal_vs_axis)
    frames = []
    for i in range(n_nuc):
        phi = np.radians(phase0 + i * deg)
        radial = np.cos(phi) * N1[i] + np.sin(phi) * N2[i]
        radial /= np.linalg.norm(radial)
        centre = pts[i] + radius * radial

        y = -radial                       # dyad points at the axis
        t = T[i]
        perp = np.cross(t, radial)
        np_ = np.linalg.norm(perp)
        if np_ < 1e-9:
            continue
        perp /= np_
        z = np.cos(beta) * t + np.sin(beta) * perp
        z -= np.dot(z, y) * y
        z /= np.linalg.norm(z)
        # Rotating the nucleosome 180 deg about its dyad leaves the dyad
        # pointing inward either way, but puts the DNA ends on opposite sides.
        # The four sign combinations give exit->entry spans of 31.5 / 52.6 /
        # 91.3 / 117.1 A for the same 136.4 A centre spacing. Placing 1KX5 on
        # b2's OWN measured frames gives 54.6 A, so this is the branch that
        # matches the reference array; the default (unflipped) was the worst.
        if flip_normal:
            z = -z
        x = np.cross(y, z)
        frames.append((centre, np.column_stack([x, y, z])))
    return frames, total
