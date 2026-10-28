"""
Full-atom 9-nucleosome array, v2 -- four faults fixed.

  1. ADAPTIVE FIBRE SCALE. v1 ran NRL 190 on b2's own zigzag scale, giving
     extension 0.36 (147 A of contour into a 52.6 A gap). The closure solver
     looped the slack over the next core -- bp 22 of 43 passed 41.2 A from its
     centre, inside the 41.8 A DNA gyre, and 4.3 A from its histones. That is
     the "linker aims at the centre of the next nucleosome" artifact. Here the
     zigzag radius and rise follow the LOCAL linker length (fibre_adaptive).

  2. O5' IS EMITTED. v1 computed O5' to place each phosphate and then never
     wrote it, so 3,300 nucleotides carried P with nothing bridging P to C5'.

  3. JUNCTION PHOSPHATES. v1's solver had no boundary conditions, so all 16
     core<->linker junctions were unbuilt and every linker duplex floated free
     at both ends. Now the flanking core nucleotides are passed in as fixed
     boundaries (backbone_solver.solve_strand_bounded).

  4. CHAIN B ORDER. 1KX5 chain J runs 5'->3' from the exit side to the entry
     side, so the partner strand runs BACKWARDS along the fibre. v1 wrote its
     blocks in nucleosome order, implying 64 A jumps between consecutive
     residue numbers. Correct order is nuc8, linker7, nuc7, ..., linker0, nuc0.

Nucleosome cores are 1KX5 placed rigidly, so histone/DNA geometry inside a core
is real by construction. Everything between cores is built.
"""
import sys
import numpy as np
import gemmi

sys.path.insert(0, "/tmp/dna_poc/bpstep")
import bpstep as bs
import fibre_fit as ff
import fibre_adaptive as fa
import closure as cl
import backbone_solver as solver
from build_zigzag import canonical_1kx5, place
from rebuild_test import build_templates, collect, load

KX5 = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'data', '1KX5.cif')
AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
OUT = "/mnt/user-data/outputs/array9_full_atom_v2.cif"

# A realistic spread rather than a single NRL, so the adaptive scaling is
# actually exercised the way it will be on nucleosome-called data. All are at
# inside the window this zigzag actually supports, NRL 163-205
# (fibre_adaptive.supported_nrl). The window is set by AIM, not reach: the
# hard reachability floor is NRL 160, but at that spacing the duplex leaves
# the core 43 deg off its direction of travel and has to turn away and back.
# Outside the window is not a solver failure, it is a different local
# geometry -- stacked rather than open zigzag.
NRL = [167, 182, 172, 193, 165, 177, 170, 188]
WRAP = 147
LINK_BP = [n - WRAP for n in NRL]
N_NUC = len(NRL) + 1

# relaxed B-DNA, used for the two steps that bridge core and linker
IDEAL_STEP = np.array([0.0, 0.0, cl.B_RISE, 0.0, 0.0, cl.B_TWIST])


def back_step(frame_j, params):
    """Inverse of bpstep.compose_step: the frame one step BEFORE frame_j."""
    o_j = np.asarray(frame_j[0], float)
    R_j = np.asarray(frame_j[1], float)
    th = np.radians(params[3:6])
    R_i = R_j @ cl._rodrigues(th).T
    return (o_j - R_i @ (cl._rodrigues(th / 2.0) @ params[0:3]), R_i)


BASE_ONLY = {"A": ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],
             "G": ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
             "C": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
             "T": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C7"]}


def kx5_parts():
    """1KX5 in its own sign-fixed nucleosome frame, split into protein and the
    two DNA strands, keeping seqid so the terminal bp can be found again."""
    st = gemmi.read_structure(KX5); st.setup_entities(); m = st[0]
    wrapped, dyad = [], None
    for r in m["I"]:
        if bs.is_dna(r.name):
            for a in r:
                if a.name == "C1'":
                    wrapped.append([a.pos.x, a.pos.y, a.pos.z])
                    if r.seqid.num == 0:
                        dyad = np.array([a.pos.x, a.pos.y, a.pos.z])
    o_n, R_n = ff.nucleosome_frame(np.array(wrapped), dyad)

    parts = {"protein": [], "dnaI": [], "dnaJ": []}
    for ch in m:
        for r in ch:
            if r.name in AA:
                key = "protein"
            elif bs.is_dna(r.name) and ch.name in ("I", "J"):
                key = "dnaI" if ch.name == "I" else "dnaJ"
            else:
                continue
            atoms = {a.name: R_n.T @ (np.array([a.pos.x, a.pos.y, a.pos.z]) - o_n)
                     for a in r}
            parts[key].append((ch.name, r.name, r.seqid.num, atoms))
    return parts


def main():
    print(f"9-nucleosome array v2, NRL {NRL}")
    print(f"  linkers: {LINK_BP} bp\n")

    loc_entry, loc_exit = canonical_1kx5()

    # ---- level 2: adaptive zigzag -----------------------------------------
    st = fa.b2_step()
    print(f"  b2 step transform: {st['L0']:.1f} A, {st['rot_deg']:.2f} deg")
    print(f"  adaptive step length (target extension {fa.TARGET_EXT}):")
    frames, step_L, spans, ext, notes = fa.solve_fibre(
        LINK_BP, place, loc_entry, loc_exit, verbose=True)
    if any(notes):
        raise RuntimeError(f"unrepresentable linkers: {notes}")
    C = np.array([f[0] for f in frames])
    D = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=-1)
    iu = np.triu_indices(len(frames), 1)
    print(f"  min centre-centre over all pairs: {D[iu].min():.1f} A\n")

    # ---- templates and conformer library ----------------------------------
    parts = kx5_parts()
    model = load()
    src_frames, seq, kept, ref, par = collect(model)
    templates, _ = build_templates(src_frames, kept, ref, par)
    lib = solver.prune_library(solver.extract_library(model))
    print(f"  sugar conformer library: {len(lib)} (pruned)\n")

    # ---- place the cores ---------------------------------------------------
    coreA, coreB, protein, term = [], [], [], []
    for ni, (o_n, R_n) in enumerate(frames):
        def world(v): return R_n @ v + o_n
        protein.append([(ch, rn, sq, {n: world(v) for n, v in at.items()})
                        for ch, rn, sq, at in parts["protein"]])
        blocks = {}
        for key, dest in (("dnaI", coreA), ("dnaJ", coreB)):
            blk = [[rn, sq, {n: world(v) for n, v in at.items()}]
                   for ch, rn, sq, at in parts[key]]
            dest.append(blk)
            blocks[key] = blk
        # terminal nucleotides, by seqid, for the junction boundary conditions
        byI = {b[1]: b[2] for b in blocks["dnaI"]}
        byJ = {b[1]: b[2] for b in blocks["dnaJ"]}
        term.append({"I_exit": byI[73], "I_entry": byI[-73],
                     "J_entry3": byJ[73], "J_exit5": byJ[-73]})

    # ---- linkers -----------------------------------------------------------
    linkA, linkB = [], []
    jstats = []
    for ni in range(N_NUC - 1):
        nb = LINK_BP[ni]
        ex = place(frames[ni], loc_exit)
        en = place(frames[ni + 1], loc_entry)
        obst = cl.nucleosome_obstacles(frames, donor=ni, acceptor=ni + 1)
        # The two steps that bridge core and linker are pinned to relaxed
        # B-DNA rather than left to the optimiser. Letting the closure solver
        # choose them put arbitrary roll, tilt and twist right at the junction,
        # and no sugar conformer can then bridge O3' to a FIXED core O5':
        # measured, the terminal linkages came out at 3.31 +/- 1.45 A against
        # real DNA's 2.487 +/- 0.018. Free linker DNA relaxes to B form anyway,
        # so pinning the junction steps costs nothing and the interior absorbs
        # all the bending.
        ex1 = bs.compose_step(ex, IDEAL_STEP)          # first linker bp
        en1 = back_step(en, IDEAL_STEP)                # last linker bp
        p, err, _ = cl.solve_path_multistart(
            ex1, en1, nb - 1, good_enough=cl.ELASTIC_BEND_PER_STEP,
            obstacles=obst)
        q = cl.path_quality(ex1, p, obst)
        bpf = bs.rebuild_frames(ex1, p)
        clr = q["clearance_A"]
        print(f"  linker {ni}: closure {np.linalg.norm(err[:3]):.1e} A, "
              f"bend {q['bend_max_deg']:5.1f} deg/step, clearance {clr:6.2f} A, "
              f"{len(bpf)} bp"
              + ("" if clr >= cl.DUPLEX_R else "   <-- TIGHT"))
        if np.linalg.norm(err[:3]) > 1e-3:
            raise RuntimeError(f"linker {ni} failed to close ({err[:3]})")

        lets = (["A", "T", "G", "C"] * 30)[:len(bpf)]
        gl = {"R": [], "P": []}
        bl = {"R": [], "P": []}
        for f, L in zip(bpf, lets):
            t = templates[L]
            wa = {k: bs.to_world(f, v[None, :])[0] for k, v in t.items()}
            for tag in ("R", "P"):
                letter = L if tag == "R" else bs.COMPLEMENT[L]
                keep = {nm: wa[(tg, nm)] for (tg, nm) in wa if tg == tag
                        and (nm in BASE_ONLY[letter] or nm == "C1'")}
                gl[tag].append(solver.glyc_frame(keep, letter))
                bl[tag].append((letter, keep))

        # Reference strand (1KX5 chain I-like) runs FORWARD along the fibre:
        #   upstream = this nucleosome's bp +73, downstream = next one's bp -73.
        # Partner strand (chain J-like) runs BACKWARDS, so its order is the
        # reverse of the bp index and its neighbours are the other way round.
        seg = {}
        for tag, dest in (("R", linkA), ("P", linkB)):
            glf, blist = gl[tag], bl[tag]
            if tag == "P":
                glf = glf[::-1]; blist = blist[::-1]
                up, dn = term[ni + 1]["J_entry3"], term[ni]["J_exit5"]
            else:
                up, dn = term[ni]["I_exit"], term[ni + 1]["I_entry"]
            chosen = solver.solve_strand_bounded(
                glf, [b[0] for b in blist], lib,
                prev_O3=up["O3'"], prev_C3=up["C3'"],
                next_O5=dn["O5'"], next_C5=dn["C5'"], next_C4=dn["C4'"])
            # the solver is allowed to spin the acceptor's dangling 5'-OH;
            # adopt its choice before closing the junction
            if "next_O5" in chosen[-1]:
                dn["O5'"] = chosen[-1]["next_O5"]
            links = solver.assemble(chosen)

            blk = []
            for i, (ci, (letter, batoms)) in enumerate(zip(chosen, blist)):
                at = dict(batoms)
                at.update(ci["sugar"])
                if ci["O5"] is not None:
                    at["O5'"] = ci["O5"]
                if i > 0:
                    lk = links[i - 1]
                    for nm in ("P", "OP1", "OP2"):
                        at[nm] = lk[nm]
                blk.append(["D" + letter, i + 1, at])

            # junction 1: upstream core -> linker nucleotide 0
            j1 = solver.close_junction(up["O3'"], up["C3'"],
                                       chosen[0]["O5"], chosen[0]["sugar"]["C5'"])
            for nm in ("P", "OP1", "OP2"):
                blk[0][2][nm] = j1[nm]
            # junction 2: linker last nucleotide -> downstream core
            last = chosen[-1]["sugar"]
            j2 = solver.close_junction(last["O3'"], last["C3'"],
                                       dn["O5'"], dn["C5'"])
            for nm in ("P", "OP1", "OP2"):
                dn[nm] = j2[nm]
            dest.append(blk)
            jstats += [(ni, tag, "5'", j1), (ni, tag, "3'", j2)]
            seg[tag] = (j1, j2)

    # ---- emit --------------------------------------------------------------
    st = gemmi.Structure(); st.name = "array9_v2"
    gm = gemmi.Model("1")

    def make_chain(name, blocks):
        ch = gemmi.Chain(name)
        sid = 1
        for blk in blocks:
            for rn, _, at in blk:
                res = gemmi.Residue(); res.name = rn
                res.seqid = gemmi.SeqId(sid, " "); sid += 1
                for nm, xyz in at.items():
                    a = gemmi.Atom(); a.name = nm
                    a.element = gemmi.Element(nm[0])
                    a.pos = gemmi.Position(*xyz); a.occ = 1.0; a.b_iso = 30.0
                    res.add_atom(a)
                ch.add_residue(res)
        return ch

    seqA, seqB = [], []
    for ni in range(N_NUC):
        seqA.append(coreA[ni])
        if ni < N_NUC - 1:
            seqA.append(linkA[ni])
    # partner strand 5'->3' runs from the LAST nucleosome back to the first
    for ni in range(N_NUC - 1, -1, -1):
        seqB.append(coreB[ni])
        if ni > 0:
            seqB.append(linkB[ni - 1])

    gm.add_chain(make_chain("A", seqA))
    gm.add_chain(make_chain("B", seqB))
    for ni, plist in enumerate(protein):
        byc = {}
        for ch, rn, sq, at in plist:
            byc.setdefault(ch, []).append((rn, sq, at))
        for ch, rows in byc.items():
            c = gemmi.Chain(f"P{ni}{ch}")
            for rn, sq, at in rows:
                res = gemmi.Residue(); res.name = rn
                res.seqid = gemmi.SeqId(sq, " ")
                for nm, xyz in at.items():
                    a = gemmi.Atom(); a.name = nm
                    a.element = gemmi.Element(nm[0])
                    a.pos = gemmi.Position(*xyz); a.occ = 1.0; a.b_iso = 30.0
                    res.add_atom(a)
                c.add_residue(res)
            gm.add_chain(c)

    st.add_model(gm)
    st.setup_entities()
    st.make_mmcif_document().write_file(OUT)
    n_atoms = sum(len(r) for c in gm for r in c)
    print(f"\nwrote {OUT}  ({len(gm)} chains, {n_atoms} atoms)")

    d = np.array([j[3]["d_O3_O5"] for j in jstats])
    a = np.array([j[3]["O3_P_O5"] for j in jstats])
    print(f"\n  {len(jstats)} core<->linker junctions built")
    print(f"    O3'...O5' : {d.mean():.3f} +/- {d.std():.3f} A  "
          f"[real 2.487 +/- 0.018]  range {d.min():.2f}-{d.max():.2f}")
    print(f"    O3'-P-O5' : {a.mean():.1f} +/- {a.std():.1f} deg  "
          f"[real 102.5 +/- 0.9]  range {a.min():.1f}-{a.max():.1f}")
    return OUT


if __name__ == "__main__":
    main()
