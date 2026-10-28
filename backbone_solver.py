"""
Backbone construction with geometry as an INPUT, not an outcome.

Measured on 1KX5: bond lengths vary by +/-0.005 A and bond angles by
+/-0.4-1.3 deg, while torsions swing 22-36 deg. All the conformational
freedom is torsional. And the linkage tolerance is brutal -- real DNA holds
O3'...O5' at 2.487 +/- 0.018 A, and the O3'-P-O5' angle is so sensitive to
that separation (2.4 A -> 97 deg, 3.0 A -> 139 deg) that a rigid-template
placement error of +/-0.44 A destroys it. That is 25x too coarse, and it is
why the rigid-nucleoside rebuild failed even at 0.47 A atom RMSD.

So: never place backbone atoms and hope the bonds come out right.

    base      -- fixed by the base-pair frame
    C1'       -- fixed by the base (it sits ON the glycosidic axis)
    sugar     -- a REAL conformer drawn from a 1KX5-derived library, so ring
                 geometry and pucker are always physically valid by
                 construction; no idealisation, no ring-closure error
    chi       -- continuous rotation of that sugar about the glycosidic axis
    gamma     -- continuous rotation placing O5' off C5'
    phosphate -- placed analytically to make an exact phosphodiester

The only free quantities are torsions, which is exactly where real DNA puts
its freedom. Every bond length and bond angle in the output is an input
value, so the 16-deg O3'-P-O5' artifact is not merely unlikely -- it is
unreachable.
"""
import numpy as np
import bpstep as bs

SUGAR_ATOMS = ["C2'", "C3'", "C4'", "O4'", "O3'", "C5'"]
GLYC_N = {"A": "N9", "G": "N9", "C": "N1", "T": "N1"}
RING_REF = {"A": "C4", "G": "C4", "C": "C2", "T": "C2"}

TARGET_O3_O5 = 2.487      # measured in 1KX5
TOL_O3_O5 = 0.018
# orientation proxies, also measured in 1KX5. Constraining these pins the two
# sugars' relative orientation, which the separation alone does not do.
PROXY_C3_O3_O5 = 129.33
PROXY_TOL_IN = 15.25
PROXY_O3_O5_C5 = 94.14
PROXY_TOL_OUT = 7.17
ORIENT_W = 0.08          # tie-breaker weight; see solve_strand


def glyc_frame(atoms, letter):
    """Frame anchored on the glycosidic bond: origin C1', z along C1'->N.
    Rotating about z IS the chi torsion, and C1' sits on the axis so it
    never moves."""
    C1 = atoms["C1'"]
    N = atoms[GLYC_N[letter]]
    ref = atoms[RING_REF[letter]]

    z = N - C1
    z /= np.linalg.norm(z)
    x = ref - C1
    x = x - np.dot(x, z) * z
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return C1, np.column_stack([x, y, z])


def extract_library(model, chains=("I", "J")):
    """Every real sugar in the structure, expressed in its own glycosidic
    frame. These are real conformers -- correct pucker, correct ring."""
    lib = []
    for cname in chains:
        for r in model[cname]:
            if not bs.is_dna(r.name):
                continue
            a = bs.atom_dict(r)
            letter = bs.base_letter(r.name)
            if "C1'" not in a or GLYC_N[letter] not in a or RING_REF[letter] not in a:
                continue
            if not all(n in a for n in SUGAR_ATOMS):
                continue
            o, R = glyc_frame(a, letter)
            local = {n: R.T @ (a[n] - o) for n in SUGAR_ATOMS}
            entry = {"letter": letter, "local": local}
            if "O5'" in a:
                entry["O5_local"] = R.T @ (a["O5'"] - o)
            lib.append(entry)
    return lib


def place_sugar(frame_o, frame_R, entry, chi_offset):
    """Put a library sugar onto a fixed base, rotated by chi_offset about the
    glycosidic axis."""
    c, s = np.cos(chi_offset), np.sin(chi_offset)
    Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    M = frame_R @ Rz
    return {n: frame_o + M @ v for n, v in entry["local"].items()}


def place_O5(sugar, entry, frame_o, frame_R, chi_offset, gamma_offset):
    """O5' hangs off C5' with fixed bond and angle; gamma rotates it about
    the C4'-C5' bond. Bond length and angle are inputs, so they are exact."""
    C5, C4 = sugar["C5'"], sugar["C4'"]
    axis = C5 - C4
    axis /= np.linalg.norm(axis)

    c, s = np.cos(chi_offset), np.sin(chi_offset)
    Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    O5_ref = frame_o + (frame_R @ Rz) @ entry["O5_local"]

    v = O5_ref - C5
    v_par = np.dot(v, axis) * axis
    v_perp = v - v_par
    cg, sg = np.cos(gamma_offset), np.sin(gamma_offset)
    v_rot = v_par + cg * v_perp + sg * np.cross(axis, v_perp)
    return C5 + v_rot


def solve_strand(base_frames, letters, library, lookahead=True, verbose=False):
    """Walk the strand, choosing (conformer, chi, gamma) per nucleotide so each
    O3'...O5' linkage lands on the real-DNA target separation.

    base_frames: list of (origin, R) glycosidic frames, one per nucleotide,
    in 5'->3' order. These come from the base placement and are FIXED.
    """
    n = len(base_frames)
    cand_by_letter = {}
    for e in library:
        cand_by_letter.setdefault(e["letter"], []).append(e)

    CHI_GRID = np.radians(np.arange(-25, 25.1, 2.5))
    GAM_GRID = np.radians(np.arange(-40, 40.1, 5.0))

    chosen = []
    prev_O3 = None
    prev_C3 = None
    for i in range(n):
        o, R = base_frames[i]
        letter = letters[i]
        cands = cand_by_letter.get(letter, library)

        best = None
        for entry in cands:
            if "O5_local" not in entry:
                continue
            for chi in CHI_GRID:
                sugar = place_sugar(o, R, entry, chi)
                if prev_O3 is None:
                    cost = abs(chi) * 0.5
                    rec = (cost, entry, chi, 0.0, sugar, None)
                    if best is None or rec[0] < best[0]:
                        best = rec
                    continue
                for gam in GAM_GRID:
                    O5 = place_O5(sugar, entry, o, R, chi, gam)
                    d = np.linalg.norm(prev_O3 - O5)
                    # primary: hit the real linkage separation
                    cost = ((d - TARGET_O3_O5) / TOL_O3_O5) ** 2
                    # also constrain the two sugars' relative ORIENTATION, not
                    # just their separation. Without this the phosphate is free
                    # on a circle and the flanking angles drift ~10 deg. These
                    # proxies need no phosphate placement, so they are cheap.
                    # Weighted LOW on purpose: the separation is a hard
                    # constraint (the O3'-P-O5' angle is hypersensitive to it),
                    # whereas orientation is a preference. At equal weight the
                    # solver trades away distance accuracy and the primary angle
                    # degrades from 102.1 +/- 0.6 to 103.2 +/- 5.8. As a
                    # tie-breaker among candidates that already satisfy the
                    # distance, it improves orientation for free.
                    ang_in = bs.angle(prev_C3, prev_O3, O5)
                    ang_out = bs.angle(prev_O3, O5, sugar["C5'"])
                    cost += ORIENT_W * ((ang_in - PROXY_C3_O3_O5) / PROXY_TOL_IN) ** 2
                    cost += ORIENT_W * ((ang_out - PROXY_O3_O5_C5) / PROXY_TOL_OUT) ** 2
                    # secondary: stay near the conformer's own real torsions
                    cost += 0.15 * (np.degrees(chi) / 10.0) ** 2
                    cost += 0.05 * (np.degrees(gam) / 20.0) ** 2
                    if best is None or cost < best[0]:
                        best = (cost, entry, chi, gam, sugar, O5)
        cost, entry, chi, gam, sugar, O5 = best
        chosen.append({"entry": entry, "chi": chi, "gamma": gam,
                       "sugar": sugar, "O5": O5, "cost": cost})
        prev_O3 = sugar["O3'"]
        prev_C3 = sugar["C3'"]
        if verbose and i % 25 == 0:
            print(f"   nt {i:3d}: cost {cost:8.2f}")
    return chosen


def assemble(chosen):
    """Turn the solved per-nucleotide choices into a full backbone, placing
    each phosphate analytically. Returns per-link geometry measurements."""
    links = []
    for i in range(len(chosen) - 1):
        a, b = chosen[i], chosen[i + 1]
        O3, C3 = a["sugar"]["O3'"], a["sugar"]["C3'"]
        O5, C5 = b["O5"], b["sugar"]["C5'"]
        if O5 is None:
            continue
        P, OP1, OP2, strain = bs.place_phosphate(O3, C3, O5, C5)
        links.append({
            "P": P, "OP1": OP1, "OP2": OP2,
            "d_O3_O5": float(np.linalg.norm(O3 - O5)),
            "strain": strain,
            "O3_P_O5": bs.angle(O3, P, O5),
            "C3_O3_P": bs.angle(C3, O3, P),
            "P_O5_C5": bs.angle(P, O5, C5),
            "bond_O3_P": float(np.linalg.norm(O3 - P)),
            "bond_P_O5": float(np.linalg.norm(P - O5)),
        })
    return links


def report(links, label):
    d = np.array([l["d_O3_O5"] for l in links])
    a1 = np.array([l["O3_P_O5"] for l in links])
    a2 = np.array([l["C3_O3_P"] for l in links])
    a3 = np.array([l["P_O5_C5"] for l in links])
    b1 = np.array([l["bond_O3_P"] for l in links])
    b2 = np.array([l["bond_P_O5"] for l in links])
    print(f"\n{label}   n={len(links)}")
    print(f"  O3'...O5' separation : {d.mean():.3f} +/- {d.std():.3f} A"
          f"   [real DNA 2.487 +/- 0.018]")
    print(f"  O3'-P bond           : {b1.mean():.3f} +/- {b1.std():.3f} A   [1.601 +/- 0.005]")
    print(f"  P-O5' bond           : {b2.mean():.3f} +/- {b2.std():.3f} A   [1.589 +/- 0.004]")
    print(f"  O3'-P-O5' angle      : {a1.mean():6.1f} +/- {a1.std():4.1f} deg"
          f"   range [{a1.min():.1f}, {a1.max():.1f}]   [102.5 +/- 0.9]")
    print(f"  C3'-O3'-P angle      : {a2.mean():6.1f} +/- {a2.std():4.1f} deg   [119.9 +/- 0.4]")
    print(f"  P-O5'-C5' angle      : {a3.mean():6.1f} +/- {a3.std():4.1f} deg   [120.0 +/- 0.6]")
    bad = int(np.sum((a1 < 95) | (a1 > 115)))
    print(f"  junctions outside 95-115 deg: {bad} / {len(links)}")
    return a1


def prune_library(library, max_per_letter=36):
    """Cluster the conformer library down to representatives.

    294 real sugars carry a lot of redundancy, and solve_strand's cost is
    linear in library size. Clustering to medoids keeps genuine pucker
    diversity while making a full array tractable.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    out = []
    by_letter = {}
    for e in library:
        by_letter.setdefault(e["letter"], []).append(e)
    for letter, entries in by_letter.items():
        if len(entries) <= max_per_letter:
            out.extend(entries); continue
        X = np.array([np.concatenate([e["local"][n] for n in SUGAR_ATOMS])
                      for e in entries])
        lab = fcluster(linkage(X, "ward"), max_per_letter, criterion="maxclust")
        for k in range(1, max_per_letter + 1):
            idx = np.where(lab == k)[0]
            if not len(idx):
                continue
            sub = X[idx]
            medoid = idx[int(np.argmin(np.linalg.norm(sub - sub.mean(0), axis=1)))]
            out.append(entries[medoid])
    return out


def solve_strand_fast(base_frames, letters, library, chi_step=5.0, gam_step=10.0):
    """solve_strand on coarser torsion grids. The fine grid buys little once
    the conformer library supplies most of the geometric variety."""
    return solve_strand_bounded(base_frames, letters, library,
                                chi_step=chi_step, gam_step=gam_step)


def _link_cost(O3, C3, O5, C5):
    """Cost of an O3'...O5' linkage: separation (hard) plus orientation (soft).

    Separation is weighted hard because the O3'-P-O5' angle is hypersensitive
    to it -- 2.40 A -> 97.2 deg, 3.00 A -> 139.3 deg, 3.19 A -> 170.9 deg.
    """
    d = np.linalg.norm(O3 - O5)
    cost = ((d - TARGET_O3_O5) / TOL_O3_O5) ** 2
    ai = bs.angle(C3, O3, O5)
    ao = bs.angle(O3, O5, C5)
    cost += ORIENT_W * ((ai - PROXY_C3_O3_O5) / PROXY_TOL_IN) ** 2
    cost += ORIENT_W * ((ao - PROXY_O3_O5_C5) / PROXY_TOL_OUT) ** 2
    return cost


def rotate_O5(O5_ref, C5, C4, gamma):
    """Rotate an O5' about the C4'-C5' bond. This IS the gamma torsion."""
    axis = C5 - C4
    axis = axis / np.linalg.norm(axis)
    v = O5_ref - C5
    v_par = np.dot(v, axis) * axis
    v_perp = v - v_par
    return C5 + v_par + np.cos(gamma) * v_perp + np.sin(gamma) * np.cross(axis, v_perp)


def _link_cost_many(O3, C3, O5s, C5):
    """_link_cost against many candidate O5' positions at once."""
    d = np.linalg.norm(O5s - O3, axis=1)
    cost = ((d - TARGET_O3_O5) / TOL_O3_O5) ** 2

    def ang(p, q, r):
        u = p - q
        v = r - q
        u = u / np.linalg.norm(u, axis=-1, keepdims=True)
        v = v / np.linalg.norm(v, axis=-1, keepdims=True)
        return np.degrees(np.arccos(np.clip(np.sum(u * v, axis=-1), -1.0, 1.0)))

    ai = ang(C3[None, :], O3[None, :], O5s)
    ao = ang(O3[None, :], O5s, C5[None, :])
    cost = cost + ORIENT_W * ((ai - PROXY_C3_O3_O5) / PROXY_TOL_IN) ** 2
    cost = cost + ORIENT_W * ((ao - PROXY_O3_O5_C5) / PROXY_TOL_OUT) ** 2
    return cost


def solve_strand_bounded(base_frames, letters, library,
                         prev_O3=None, prev_C3=None,
                         next_O5=None, next_C5=None, next_C4=None,
                         chi_step=5.0, gam_step=10.0, terminal_w=1.0):
    """solve_strand_fast with FIXED FLANKING NUCLEOTIDES.

    Why this exists
    ---------------
    The unbounded solver has no boundary condition at either end of a segment.
    The first nucleotide therefore gets no O5' at all, and the last one's O3'
    is placed by whatever conformer its own cost happened to pick, with nothing
    pulling it onto the next residue. In the 9-nucleosome array that left all
    16 core<->linker junction phosphates unbuilt, so every linker duplex was
    covalently FREE AT BOTH ENDS -- the "floating fragments". The terminal
    O3'...O5' separations came out at 0.66 A on the reference strand where real
    DNA holds 2.487 +/- 0.018, i.e. an overlap, not a bond.

    prev_O3 / prev_C3 : the upstream (5') neighbour's O3' and C3'. Given these,
        nucleotide 0 is solved like any interior one and gets a real phosphate.
    next_O5 / next_C5 : the downstream (3') neighbour's O5' and C5'. Given
        these, the LAST nucleotide additionally pays for its own O3' landing on
        that O5'.

    The two terminal constraints are nearly decoupled in the search: gamma
    moves O5' (the upstream link) while conformer and chi move O3' (the
    downstream link), so one nucleotide can satisfy both.
    """
    n = len(base_frames)
    cand = {}
    for e in library:
        cand.setdefault(e["letter"], []).append(e)
    CHI = np.radians(np.arange(-25, 25.1, chi_step))
    GAM = np.radians(np.arange(-40, 40.1, gam_step))

    # The downstream core nucleotide's O5' is the 5'-terminal hydroxyl of
    # 1KX5's chain -- a dangling OH whose gamma torsion carries no structural
    # information, because in the real fibre the DNA does not end there. Left
    # fixed it makes the last junction a one-sided problem, with only the
    # linker's conformer and chi to bridge a frozen target; measured, that gave
    # terminal O3'...O5' of 2.78 +/- 0.71 A and O3'-P-O5' angles up to 180 deg.
    # Rotating it is free and restores the two-sided freedom every interior
    # linkage already has.
    GAM_TERM = np.radians(np.arange(0.0, 360.0, 2.0))
    free_next = next_O5 is not None and next_C4 is not None
    # The ring of reachable acceptor O5' positions depends only on the
    # acceptor's own atoms, so build it once instead of per candidate.
    next_ring = (np.array([rotate_O5(next_O5, next_C5, next_C4, g)
                           for g in GAM_TERM]) if free_next else None)

    chosen = []
    for i in range(n):
        o, R = base_frames[i]
        is_last = (i == n - 1)
        term = is_last and next_O5 is not None
        best = None
        best_next_O5 = None
        for entry in cand.get(letters[i], library):
            if "O5_local" not in entry:
                continue
            for chi in CHI:
                sugar = place_sugar(o, R, entry, chi)
                base_cost = 0.15 * (np.degrees(chi) / 10.0) ** 2
                this_next_O5 = None
                if term and free_next:
                    costs = _link_cost_many(sugar["O3'"], sugar["C3'"],
                                            next_ring, next_C5)
                    k = int(np.argmin(costs))
                    this_next_O5 = next_ring[k]
                    base_cost += terminal_w * float(costs[k])
                elif term:
                    base_cost += terminal_w * _link_cost(
                        sugar["O3'"], sugar["C3'"], next_O5, next_C5)
                if prev_O3 is None:
                    c = base_cost + abs(chi) * 0.5
                    if best is None or c < best[0]:
                        best = (c, entry, chi, 0.0, sugar, None)
                        best_next_O5 = this_next_O5
                    continue
                for gam in GAM:
                    O5 = place_O5(sugar, entry, o, R, chi, gam)
                    cost = base_cost + _link_cost(
                        prev_O3, prev_C3, O5, sugar["C5'"])
                    cost += 0.05 * (np.degrees(gam) / 20.0) ** 2
                    if best is None or cost < best[0]:
                        best = (cost, entry, chi, gam, sugar, O5)
                        best_next_O5 = this_next_O5
        cost, entry, chi, gam, sugar, O5 = best
        rec = {"entry": entry, "chi": chi, "gamma": gam,
               "sugar": sugar, "O5": O5, "cost": cost}
        if is_last and best_next_O5 is not None:
            # the caller must write this back onto the downstream core residue
            rec["next_O5"] = best_next_O5
        chosen.append(rec)
        prev_O3, prev_C3 = sugar["O3'"], sugar["C3'"]
    return chosen


def close_junction(O3, C3, O5, C5):
    """Place the phosphate of one junction and report its geometry."""
    P, OP1, OP2, strain = bs.place_phosphate(O3, C3, O5, C5)
    return {"P": P, "OP1": OP1, "OP2": OP2,
            "d_O3_O5": float(np.linalg.norm(O3 - O5)), "strain": strain,
            "O3_P_O5": bs.angle(O3, P, O5), "C3_O3_P": bs.angle(C3, O3, P),
            "P_O5_C5": bs.angle(P, O5, C5),
            "bond_O3_P": float(np.linalg.norm(O3 - P)),
            "bond_P_O5": float(np.linalg.norm(P - O5))}
