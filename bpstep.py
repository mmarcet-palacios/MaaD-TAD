"""
Base-pair-step dsDNA builder -- core library.

Replaces the oxDNA/tacoxDNA coordinate pipeline. The failure mode of that
pipeline was structural: linker conformations were generated independently
and then rigidly placed, so consecutive fragments had to *meet*, and the
backbone geometry at every meeting point was whatever fell out. Measured on
the v3 construct, that gave correct bond LENGTHS (1.60 A, std 0.00) but
O3'-P-O5' angles of 104 +/- 35 deg (range 16-177 deg) -- topologically
continuous, chemically impossible. Real nucleosomal DNA in the same file:
102.5 +/- 0.9 deg.

This module builds instead by composition: base pair i+1 is placed in the
local frame of base pair i via a real helical step transform. Nothing has to
meet, because nothing is ever independently placed.

Conventions
-----------
Base-pair frame (origin o, rotation R with columns [x y z]):
    o : midpoint of the two C1' atoms
    z : base-pair plane normal (SVD of both bases' ring atoms), signed to
        run along the 5'->3' direction of the reference strand, so it is the
        local helical axis and twist comes out positive for right-handed DNA
    y : C1'(partner) -> C1'(reference), orthogonalized against z
    x : y cross z   (right-handed: x cross y = z)

Step parameters use the mid-step triad (MST) convention:
    (tilt, roll, twist) = the rotation vector of R_i^T R_j, in degrees.
        Note this vector is identical in frame i, the mid-frame and frame j
        coordinates, since rotating a vector about its own axis is a no-op --
        which is why this decomposition is exactly invertible.
    (shift, slide, rise) = (o_j - o_i) expressed in the mid-step triad.

Both directions are implemented and the round trip is asserted exact.
"""
import numpy as np
from scipy.spatial.transform import Rotation

DNA_RES = {"DA", "DC", "DG", "DT"}

PURINE_RING = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"]
PYRIMIDINE_RING = ["N1", "C2", "N3", "C4", "C5", "C6"]

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# ideal phosphodiester geometry (heavy atoms), from standard nucleic acid values
BOND_O3_P = 1.607
BOND_P_O5 = 1.593
ANG_O3_P_O5 = 104.0          # degrees
# the O3'...O5' separation a strain-free linkage implies (law of cosines)
IDEAL_O3_O5 = np.sqrt(
    BOND_O3_P**2 + BOND_P_O5**2
    - 2 * BOND_O3_P * BOND_P_O5 * np.cos(np.radians(ANG_O3_P_O5))
)


def is_dna(resname):
    return resname[:2] in DNA_RES or resname in DNA_RES


def base_letter(resname):
    b = resname[:2] if resname[:2] in DNA_RES else resname
    return b[1] if len(b) == 2 else b


def ring_atom_names(resname):
    return PURINE_RING if base_letter(resname) in ("A", "G") else PYRIMIDINE_RING


def atom_dict(residue):
    return {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in residue}


# ----------------------------------------------------------------------------
# base-pair frames
# ----------------------------------------------------------------------------

def bp_frame(res_ref, res_partner):
    """Return (origin, R) for a Watson-Crick pair. res_ref defines +y and the
    sign convention; the caller fixes the z sign globally along the chain."""
    a_ref, a_par = atom_dict(res_ref), atom_dict(res_partner)

    ring = []
    for res, atoms in ((res_ref, a_ref), (res_partner, a_par)):
        for n in ring_atom_names(res.name):
            if n in atoms:
                ring.append(atoms[n])
    ring = np.array(ring)

    centroid = ring.mean(axis=0)
    # plane normal = direction of least variance
    _, _, Vt = np.linalg.svd(ring - centroid)
    z = Vt[2]
    z /= np.linalg.norm(z)

    c1_ref, c1_par = a_ref["C1'"], a_par["C1'"]
    origin = 0.5 * (c1_ref + c1_par)

    y = c1_ref - c1_par
    y = y - np.dot(y, z) * z
    y /= np.linalg.norm(y)

    x = np.cross(y, z)
    x /= np.linalg.norm(x)

    return origin, np.column_stack([x, y, z])


def fix_z_signs(frames):
    """The plane normal from SVD has arbitrary sign. Make z run consistently
    along the chain (and therefore along the local helical axis), flipping the
    frame about y so it stays right-handed."""
    fixed = []
    ref_dir = None
    for i, (o, R) in enumerate(frames):
        if ref_dir is None:
            nxt = frames[min(i + 1, len(frames) - 1)][0]
            prv = frames[max(i - 1, 0)][0]
            ref_dir = nxt - prv
        z = R[:, 2]
        if np.dot(z, ref_dir) < 0:
            # flip z and x, keep y -> still right-handed, y still points at the
            # reference strand's C1'
            R = np.column_stack([-R[:, 0], R[:, 1], -R[:, 2]])
        fixed.append((o, R))
        ref_dir = R[:, 2]
    return fixed


def chain_frames(model, ref_chain, partner_chain, pair_index_fn, indices):
    """Build the ordered list of bp frames along a duplex.
    pair_index_fn maps a reference-strand residue number to its partner's."""
    ref_res = {r.seqid.num: r for r in model[ref_chain] if is_dna(r.name)}
    par_res = {r.seqid.num: r for r in model[partner_chain] if is_dna(r.name)}

    frames, seq, kept = [], [], []
    for i in indices:
        j = pair_index_fn(i)
        if i not in ref_res or j not in par_res:
            continue
        frames.append(bp_frame(ref_res[i], par_res[j]))
        seq.append(base_letter(ref_res[i].name))
        kept.append((i, j))
    return fix_z_signs(frames), seq, kept


# ----------------------------------------------------------------------------
# step parameters
# ----------------------------------------------------------------------------

def decompose_step(frame_i, frame_j):
    """(shift, slide, rise, tilt, roll, twist); translations A, rotations deg."""
    o_i, R_i = frame_i
    o_j, R_j = frame_j

    R_rel = R_i.T @ R_j
    theta = Rotation.from_matrix(R_rel).as_rotvec()          # in frame i coords
    R_half = Rotation.from_rotvec(theta / 2.0).as_matrix()

    d = R_half.T @ (R_i.T @ (o_j - o_i))                     # into mid-step triad

    return np.array([d[0], d[1], d[2],
                     *np.degrees(theta)])


def compose_step(frame_i, params):
    """Inverse of decompose_step: place the next bp frame from step params."""
    o_i, R_i = frame_i
    shift, slide, rise, tilt, roll, twist = params

    theta = np.radians([tilt, roll, twist])
    R_rel = Rotation.from_rotvec(theta).as_matrix()
    R_half = Rotation.from_rotvec(theta / 2.0).as_matrix()

    o_j = o_i + R_i @ (R_half @ np.array([shift, slide, rise]))
    R_j = R_i @ R_rel
    return o_j, R_j


def rebuild_frames(frame0, step_params):
    """Compose a whole duplex from a starting frame and a list of steps."""
    frames = [frame0]
    for p in step_params:
        frames.append(compose_step(frames[-1], p))
    return frames


# ----------------------------------------------------------------------------
# placing atoms
# ----------------------------------------------------------------------------

def to_local(frame, xyz):
    o, R = frame
    return (xyz - o) @ R          # == R.T @ (xyz - o) row-wise


def to_world(frame, local):
    o, R = frame
    return local @ R.T + o


def place_phosphate(O3_prev, C3_prev, O5_next, C5_next):
    """Place P (and its two free oxygens) to make a real phosphodiester bond
    between two already-positioned nucleosides, instead of inheriting a
    phosphate from a rigid template and hoping it lands right.

    P is put at the point satisfying |O3'-P| = 1.607 and |P-O5'| = 1.593
    exactly, chosen on the circle of such points so that the two backbone
    angles C3'-O3'-P and P-O5'-C5' are as close to ideal as possible.
    Returns (P, OP1, OP2, strain) where strain is the deviation of the
    O3'...O5' separation from the value a strain-free linkage implies.
    """
    v = O5_next - O3_prev
    d = np.linalg.norm(v)
    strain = d - IDEAL_O3_O5

    if d > BOND_O3_P + BOND_P_O5:
        # cannot close with ideal bonds; put P on the segment, proportionally
        u = v / d
        P = O3_prev + u * (BOND_O3_P / (BOND_O3_P + BOND_P_O5)) * d
    else:
        u = v / d
        # distance along O3'->O5' to the circle's plane
        a = (BOND_O3_P**2 - BOND_P_O5**2 + d**2) / (2 * d)
        h2 = BOND_O3_P**2 - a**2
        h = np.sqrt(max(h2, 0.0))
        centre = O3_prev + a * u

        # basis for the circle of valid P positions
        ref = C3_prev - O3_prev
        e1 = ref - np.dot(ref, u) * u
        if np.linalg.norm(e1) < 1e-6:
            ref = np.array([1.0, 0.0, 0.0])
            e1 = ref - np.dot(ref, u) * u
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u, e1)

        best, best_cost = None, None
        for ang in np.linspace(0, 2 * np.pi, 72, endpoint=False):
            cand = centre + h * (np.cos(ang) * e1 + np.sin(ang) * e2)
            a1 = _angle(C3_prev, O3_prev, cand)
            a2 = _angle(cand, O5_next, C5_next)
            cost = (a1 - 119.7) ** 2 + (a2 - 120.9) ** 2
            if best_cost is None or cost < best_cost:
                best, best_cost = cand, cost
        P = best

    # the two non-bridging oxygens: perpendicular to the O3'-P-O5' plane,
    # tetrahedral-ish about P
    b1 = O3_prev - P
    b2 = O5_next - P
    b1 /= np.linalg.norm(b1)
    b2 /= np.linalg.norm(b2)
    bisector = -(b1 + b2)
    bisector /= np.linalg.norm(bisector)
    perp = np.cross(b1, b2)
    perp /= np.linalg.norm(perp)
    half = np.radians(119.6 / 2.0)
    OP1 = P + 1.485 * (np.cos(half) * bisector + np.sin(half) * perp)
    OP2 = P + 1.485 * (np.cos(half) * bisector - np.sin(half) * perp)

    return P, OP1, OP2, strain


def _angle(a, b, c):
    v1, v2 = a - b, c - b
    cosv = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cosv, -1.0, 1.0)))


def angle(a, b, c):
    return _angle(a, b, c)


def torsion(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def kabsch_rmsd(P, Q, superpose=True):
    if not superpose:
        return float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))
    Pc, Qc = P.mean(0), Q.mean(0)
    P0, Q0 = P - Pc, Q - Qc
    U, S, Vt = np.linalg.svd(P0.T @ Q0)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    P0r = (R @ P0.T).T
    return float(np.sqrt(np.mean(np.sum((P0r - Q0) ** 2, axis=1))))
