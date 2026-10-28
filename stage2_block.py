"""
Stage 2: full-atom build of one block of the fibre.

Generalises build_segment.py (the validated 15-nucleosome build) to read its
frames from stage 1 instead of generating them.

    python stage2_block.py --frames frames.npz --block 7
    python stage2_block.py --frames frames.npz --block 0 --size 3     # quick test

Block b owns cores  start .. end          (start = b*size, end = start+size-1)
          and linkers start .. end        (linker i joins core i to core i+1;
                                           the model's final core has none)
Core end+1 is placed so the last linker has its acceptor, but it is NOT
written: it belongs to the next block. The one thing the last linker changes on
that core -- the phosphate (P, OP1, OP2) and O5' of chain A's first nucleotide
-- goes to block_XX_seam.json, and stage 3 applies it. Nothing else crosses a
block boundary, so blocks need no communication and can run in any order.

Residue numbering is global: chain A runs core 0 -> core N-1, chain B runs the
other way, and each block computes its offsets from the cumulative repeat
lengths, so both strands are continuous across all blocks.

Outputs (written to a temp name, then renamed, so a killed run leaves nothing
half-written and the driver's skip-if-exists logic stays correct):
  blocks/block_XX.cif        full atom: chains A, B (DNA) and P<nuc><chain> (histones)
  blocks/block_XX_qc.csv     one row per linker
  blocks/block_XX_seam.json  patch for the next block's first core
"""
import argparse, csv, json, os, sys, time
import numpy as np
import gemmi

import bpstep as bs
import closure as cl
import backbone_solver as solver
from build_zigzag import canonical_1kx5, place
from rebuild_test import build_templates, collect, load
from build_array9b import kx5_parts, back_step, IDEAL_STEP, BASE_ONLY

OBSTACLE_WINDOW = 5          # nucleosomes either side passed as obstacles
CORE_NT = 147


CA_NEED = 14.0     # duplex axis to histone CA: P sits ~9 A off-axis, + 4 A QC limit, + 1
P_NEED = 19.0      # duplex axis to a partner P: 9 A + 5 A QC limit + 5 A for the partner's own groove
CA_TEMPLATE = np.array([at["CA"] for ch, rn, sq, at in kx5_parts()["protein"] if "CA" in at])


def partner_P(outdir, size, pj, ni, F):
    """P atoms of tether linker pj, read from its already-built block (None if not built)."""
    b = pj // size
    cif = os.path.join(outdir, f"block_{b:02d}.cif")
    if b == ni // size or not os.path.exists(cif):
        return None
    dy = F["dyad_pos"].astype(int); nb = int(F["nrl"][pj]) - CORE_NT
    g0 = int(dy[pj]) + 74; g1 = g0 + nb - 1                      # genomic span of the linker
    ra = (g0 - 12771345, g1 - 12771345)                          # chain A: r = genomic - 12,771,345
    rb = (13423346 - g1, 13423346 - g0)                          # chain B: r = 13,423,346 - genomic
    out = []
    for ch in gemmi.read_structure(cif)[0]:
        lo_, hi_ = ra if ch.name == "A" else rb if ch.name == "B" else (1, 0)
        for r in ch:
            if lo_ <= r.seqid.num <= hi_:
                a_ = r.find_atom("P", "*")
                if a_:
                    out.append(a_.pos.tolist())
    return np.array(out)


def numbering(nrl):
    """Global first-residue index of every core and linker on both strands.

    Linker i carries NRL_i - 147 nucleotides (the rebuilt path; checked at
    build time). Chain A: core0, link0, core1, ..., coreN-1.
    Chain B: coreN-1, link N-2, ..., link0, core0."""
    n = len(nrl) + 1
    link = np.asarray(nrl, int) - CORE_NT
    a_core = np.zeros(n, int); a_link = np.zeros(n - 1, int)
    pos = 1
    for i in range(n):
        a_core[i] = pos; pos += CORE_NT
        if i < n - 1:
            a_link[i] = pos; pos += link[i]
    b_core = np.zeros(n, int); b_link = np.zeros(n - 1, int)
    pos = 1
    for i in range(n - 1, -1, -1):
        b_core[i] = pos; pos += CORE_NT
        if i > 0:
            b_link[i - 1] = pos; pos += link[i - 1]
    return a_core, a_link, b_core, b_link


def add_residue(chain, name, seq, atoms, b=30.0):
    res = gemmi.Residue(); res.name = name
    res.seqid = gemmi.SeqId(int(seq), " ")
    for nm, xyz in atoms.items():
        at = gemmi.Atom(); at.name = nm
        at.element = gemmi.Element(nm[0])
        at.pos = gemmi.Position(*xyz); at.occ = 1.0; at.b_iso = b
        res.add_atom(at)
    chain.add_residue(res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--block", type=int, required=True)
    ap.add_argument("--size", type=int, default=100)
    ap.add_argument("--outdir", default="blocks")
    ap.add_argument("--force", action="store_true", help="rebuild even if output exists")
    ap.add_argument("--tethers", default="data/tethers.json",
                    help="cohesin tether waypoints from tethers_layer3.py (ignored if absent)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    tag = f"block_{a.block:02d}"
    out_cif = os.path.join(a.outdir, tag + ".cif")
    if os.path.exists(out_cif) and not a.force:
        print(f"{tag}: exists, skipping"); return

    t0 = time.time()
    F = np.load(a.frames)
    O, Rm, nrl = F["origins"], F["rotations"], F["nrl"].astype(int)
    src = F["source"]
    n = len(O)
    start = a.block * a.size
    if start >= n:
        print(f"{tag}: past the end ({n} nucleosomes)"); return
    end = min(start + a.size, n) - 1
    frames = [(O[i], Rm[i]) for i in range(n)]
    link_bp = nrl - CORE_NT
    a_core, a_link, b_core, b_link = numbering(nrl)
    last = min(end, n - 2)            # last linker this block builds
    print(f"{tag}: cores {start}-{end}, linkers {start}-{last} "
          f"(of {n} nucleosomes)", flush=True)

    TW = {}                                   # linker index -> waypoint spec
    if a.tethers and os.path.exists(a.tethers):
        for t in json.load(open(a.tethers)):
            L = t["linkers"]
            for me, other in ((L[0], L[1]), (L[1], L[0])):
                TW[int(me["linker"])] = dict(t=t["id"], k=int(me["split_bp"]), nb=int(me["nb"]),
                                              w=np.array(me["waypoint"]), d=np.array(me["dir"]),
                                              partner=int(other["linker"]))
    loc_entry, loc_exit = canonical_1kx5()
    parts = kx5_parts()
    model = load()
    src_frames, seq, kept, ref, par = collect(model)
    templates, _ = build_templates(src_frames, kept, ref, par)
    lib = solver.prune_library(solver.extract_library(model))

    # ---- cores: start .. end, plus end+1 as acceptor only ----------------
    place_hi = min(end + 1, n - 1)
    coreA, coreB, protein, term = {}, {}, {}, {}
    for ni in range(start, place_hi + 1):
        o_n, R_n = frames[ni]
        def world(v, o_n=o_n, R_n=R_n): return R_n @ v + o_n
        protein[ni] = [(ch, rn, sq, {nm: world(v) for nm, v in at.items()})
                       for ch, rn, sq, at in parts["protein"]]
        blk = {}
        for key, dest in (("dnaI", coreA), ("dnaJ", coreB)):
            b = [[rn, sq, {nm: world(v) for nm, v in at.items()}]
                 for ch, rn, sq, at in parts[key]]
            dest[ni] = b; blk[key] = b
        byI = {b[1]: b[2] for b in blk["dnaI"]}
        byJ = {b[1]: b[2] for b in blk["dnaJ"]}
        term[ni] = {"I_exit": byI[73], "I_entry": byI[-73],
                    "J_entry3": byJ[73], "J_exit5": byJ[-73]}

    assert all(len(coreA[i]) == CORE_NT and len(coreB[i]) == CORE_NT for i in coreA), \
        "1KX5 core is not 147 nt per strand -- numbering would be wrong"

    # ---- linkers ----------------------------------------------------------
    from scipy.spatial import cKDTree
    from stage1_frames import rot_between
    ctree = cKDTree(O)

    def bond_max(up, blk, dn):
        """largest O3'-P distance along this linker strand incl. both junctions"""
        d = [np.linalg.norm(up["O3'"] - blk[0][2]["P"])]
        for x, y in zip(blk[:-1], blk[1:]):
            d.append(np.linalg.norm(x[2]["O3'"] - y[2]["P"]))
        d.append(np.linalg.norm(blk[-1][2]["O3'"] - dn["P"]))
        return float(max(d))

    # retry ladder: more bending modes / seeds and a stricter "good enough"
    LADDER = [dict(n_modes_bend=8, n_seeds=12, good_enough=cl.ELASTIC_BEND_PER_STEP),
              dict(n_modes_bend=12, n_seeds=36, good_enough=10.0),
              dict(n_modes_bend=16, n_seeds=72, good_enough=8.0)]

    linkA, linkB, qc, jstats = {}, {}, [], []
    for ni in range(start, last + 1):
        tl = time.time()
        nb = int(link_bp[ni])
        ex = place(frames[ni], loc_exit)
        en = place(frames[ni + 1], loc_entry)
        lo, hi = max(0, ni - OBSTACLE_WINDOW), min(n - 1, ni + 1 + OBSTACLE_WINDOW)
        obst = cl.nucleosome_obstacles(frames[lo:hi + 1], donor=ni - lo, acceptor=ni + 1 - lo)
        # SPATIAL obstacles: any nucleosome the linker could physically reach,
        # wherever it is in sequence (folded fibre, tether pinches, crossings)
        mid = 0.5 * (ex[0] + en[0])
        reach = 0.5 * np.linalg.norm(en[0] - ex[0]) + 0.35 * nb * cl.B_RISE + 120.0
        seq_near = set(range(lo, hi + 1))
        far = [k for k in ctree.query_ball_point(mid, reach) if k not in seq_near]
        obst = obst + [(O[k], Rm[k], 1) for k in sorted(far)]
        ex1 = bs.compose_step(ex, IDEAL_STEP)
        en1 = back_step(en, IDEAL_STEP)
        if ni in TW:
            tw = TW[ni]
            if tw["nb"] != nb:
                raise RuntimeError(f"linker {ni}: tethers.json expects {tw['nb']} bp, frames give {nb}")
            pj = tw["partner"]
            have = {k for k in range(lo, hi + 1)} | set(far)
            obst = obst + [(O[k], Rm[k], 1) for k in range(max(0, pj - OBSTACLE_WINDOW),
                                                             min(n, pj + 2 + OBSTACLE_WINDOW)) if k not in have]
            # atom-level keep-outs the disc model misses (stage-4 QC found both):
            #  - histone CAs (incl. N-terminal tails) of every core in the obstacle set
            #    more than 3 away in sequence: axis >= CA_NEED from any CA
            #  - the partner tether linker's P atoms, when the partner is built first
            #    (lower linker index; run_stage2.py builds those blocks in an earlier wave)
            cores = sorted({k for k in have | set(range(max(0, pj - OBSTACLE_WINDOW),
                                                           min(n, pj + 2 + OBSTACLE_WINDOW)))
                            if abs(k - ni) > 3 and abs(k - (ni + 1)) > 3})
            ca = np.array([Rm[k] @ v + O[k] for k in cores for v in CA_TEMPLATE]) if cores else np.zeros((0, 3))
            if len(ca):
                obst = obst + [("pts", cKDTree(ca), CA_NEED)]
            pp = partner_P(a.outdir, a.size, pj, ni, F) if pj < ni else None
            if pp is not None and len(pp):
                obst = obst + [("pts", cKDTree(pp), P_NEED)]
                print(f"  linker {ni}: avoiding {len(pp)} P atoms of partner linker {pj} "
                      f"and {len(ca)} histone CAs", flush=True)
            k1 = tw["k"]
            Fk = ex1
            for _ in range(k1):
                Fk = bs.compose_step(Fk, IDEAL_STEP)
            Wf = (tw["w"], rot_between(Fk[1][:, 2], tw["d"]) @ Fk[1])

        # snapshot the flanking-core atoms a strand build edits, so a retry starts clean
        snap = {k: {kk: np.array(vv) for kk, vv in term[k2][k].items()}
                for k2, k in ((ni, "J_exit5"), (ni + 1, "I_entry"))}
        best = None
        for attempt, kw in enumerate(LADDER):
            for (k2, k) in ((ni, "J_exit5"), (ni + 1, "I_entry")):
                term[k2][k].clear(); term[k2][k].update({kk: np.array(vv) for kk, vv in snap[k].items()})
            if ni in TW:
                kwa = dict(kw, n_modes_bend=min(kw["n_modes_bend"], max(2, k1 - 1)))
                kwb = dict(kw, n_modes_bend=min(kw["n_modes_bend"], max(2, nb - 2 - k1)))
                p1, e1, _ = cl.solve_path_multistart(ex1, Wf, k1, obstacles=obst, **kwa)
                p2, e2, _ = cl.solve_path_multistart(Wf, en1, nb - 1 - k1, obstacles=obst, **kwb)
                p = np.vstack([p1, p2])
                err = e1 if np.linalg.norm(e1[:3]) > np.linalg.norm(e2[:3]) else e2
                wp_err = float(np.linalg.norm(bs.rebuild_frames(ex1, p)[k1][0] - tw["w"]))
            else:
                kw2 = dict(kw, n_modes_bend=min(kw["n_modes_bend"], max(2, nb - 2)))
                p, err, _ = cl.solve_path_multistart(ex1, en1, nb - 1, obstacles=obst, **kw2)
                wp_err = None
            q = cl.path_quality(ex1, p, obst)
            bpf = bs.rebuild_frames(ex1, p)
            closure = float(np.linalg.norm(err[:3]))
            if len(bpf) != nb:
                raise RuntimeError(f"linker {ni}: {len(bpf)} bp rebuilt, expected {nb}")

            lets = (["A", "T", "G", "C"] * (len(bpf) // 4 + 1))[:len(bpf)]   # placeholder sequence
            gl = {"R": [], "P": []}; bl = {"R": [], "P": []}
            for fr, L in zip(bpf, lets):
                t = templates[L]
                wa = {k: bs.to_world(fr, v[None, :])[0] for k, v in t.items()}
                for tg in ("R", "P"):
                    letter = L if tg == "R" else bs.COMPLEMENT[L]
                    keep = {nm: wa[(tt, nm)] for (tt, nm) in wa if tt == tg
                            and (nm in BASE_ONLY[letter] or nm == "C1'")}
                    gl[tg].append(solver.glyc_frame(keep, letter))
                    bl[tg].append((letter, keep))

            built, js, worst = {}, [], 0.0
            for tg in ("R", "P"):
                glf, blist = gl[tg], bl[tg]
                if tg == "P":
                    glf = glf[::-1]; blist = blist[::-1]
                    up, dn = term[ni + 1]["J_entry3"], term[ni]["J_exit5"]
                else:
                    up, dn = term[ni]["I_exit"], term[ni + 1]["I_entry"]
                chosen = solver.solve_strand_bounded(
                    glf, [b_[0] for b_ in blist], lib,
                    prev_O3=up["O3'"], prev_C3=up["C3'"],
                    next_O5=dn["O5'"], next_C5=dn["C5'"], next_C4=dn["C4'"])
                if "next_O5" in chosen[-1]:
                    dn["O5'"] = chosen[-1]["next_O5"]
                links = solver.assemble(chosen)
                blk = []
                for i, (ci, (letter, batoms)) in enumerate(zip(chosen, blist)):
                    at = dict(batoms); at.update(ci["sugar"])
                    if ci["O5"] is not None:
                        at["O5'"] = ci["O5"]
                    if i > 0:
                        for nm in ("P", "OP1", "OP2"):
                            at[nm] = links[i - 1][nm]
                    blk.append(["D" + letter, i, at])
                j1 = solver.close_junction(up["O3'"], up["C3'"], chosen[0]["O5"], chosen[0]["sugar"]["C5'"])
                for nm in ("P", "OP1", "OP2"):
                    blk[0][2][nm] = j1[nm]
                lst = chosen[-1]["sugar"]
                j2 = solver.close_junction(lst["O3'"], lst["C3'"], dn["O5'"], dn["C5'"])
                for nm in ("P", "OP1", "OP2"):
                    dn[nm] = j2[nm]
                built[tg] = blk; js += [j1, j2]
                worst = max(worst, bond_max(up, blk, dn))
            cand = dict(attempt=attempt, p=p, q=q, closure=closure, wp_err=wp_err, built=built, js=js,
                        worst=worst, term={k: {kk: np.array(vv) for kk, vv in term[k2][k].items()}
                                           for k2, k in ((ni, "J_exit5"), (ni + 1, "I_entry"))})
            if best is None or (cand["closure"] <= 1e-3 and cand["worst"] < best["worst"]) or best["closure"] > 1e-3:
                best = cand
            if cand["closure"] <= 1e-3 and cand["worst"] <= 1.8:
                break
        # commit the best attempt
        for (k2, k) in ((ni, "J_exit5"), (ni + 1, "I_entry")):
            term[k2][k].clear(); term[k2][k].update(best["term"][k])
        linkA[ni], linkB[ni] = best["built"]["R"], best["built"]["P"]
        jstats += best["js"]
        q, closure, wp_err = best["q"], best["closure"], best["wp_err"]
        row = {"linker": ni, "nrl": int(nrl[ni]), "source": str(src[ni]),
               "closure_A": closure, "bend_max_deg": q["bend_max_deg"],
               "clearance_A": q["clearance_A"], "ok": closure <= 1e-3,
               "tether": TW[ni]["t"] if ni in TW else "",
               "waypoint_err_A": round(wp_err, 4) if wp_err is not None else "",
               "far_obstacles": len(far), "attempts": best["attempt"] + 1,
               "max_O3P_A": round(best["worst"], 3), "backbone_ok": best["worst"] <= 1.8}
        row["seconds"] = round(time.time() - tl, 1)
        qc.append(row)
        print(f"  linker {ni:4d} NRL {nrl[ni]:3d} {src[ni]:>11s}: closure {closure:.1e} A, "
              f"bend {q['bend_max_deg']:5.1f} deg/step, clearance {q['clearance_A']:6.2f} A, "
              f"O3'-P max {best['worst']:.2f} A, far obstacles {len(far)}, attempt {best['attempt']+1}, "
              f"{row['seconds']:.0f}s" + (f"  [tether {TW[ni]['t']}: waypoint hit to {wp_err:.1e} A]" if ni in TW else ""), flush=True)

    # ---- write ------------------------------------------------------------
    st = gemmi.Structure(); st.name = tag
    gm = gemmi.Model("1")
    chA, chB = gemmi.Chain("A"), gemmi.Chain("B")
    for ni in range(start, end + 1):
        for k, (rn, sq, at) in enumerate(coreA[ni]):
            add_residue(chA, rn, a_core[ni] + k, at)
        if ni in linkA:
            for rn, k, at in linkA[ni]:
                add_residue(chA, rn, a_link[ni] + k, at)
    # chain B runs core N-1 -> core 0; linker i's B strand sits between core
    # i+1 and core i. So this block's B order is:
    #   linkB[end], coreB[end], linkB[end-1], coreB[end-1], ..., linkB[start], coreB[start]
    for ni in range(end, start - 1, -1):
        if ni in linkB:
            for rn, k, at in linkB[ni]:
                add_residue(chB, rn, b_link[ni] + k, at)
        for k, (rn, sq, at) in enumerate(coreB[ni]):
            add_residue(chB, rn, b_core[ni] + k, at)
    gm.add_chain(chA); gm.add_chain(chB)
    for ni in range(start, end + 1):
        byc = {}
        for ch, rn, sq, at in protein[ni]:
            byc.setdefault(ch, []).append((rn, sq, at))
        for ch, rws in byc.items():
            c = gemmi.Chain(f"P{ni}{ch}")
            for rn, sq, at in rws:
                add_residue(c, rn, sq, at)
            gm.add_chain(c)
    st.add_model(gm); st.setup_entities()
    tmp = out_cif + ".part"
    st.make_mmcif_document().write_file(tmp)

    seam = {}
    if end + 1 < n:
        e = term[end + 1]["I_entry"]
        seam = {"core": int(end + 1), "chain": "A", "resi": int(a_core[end + 1]),
                "atoms": {nm: [float(x) for x in e[nm]] for nm in ("P", "OP1", "OP2", "O5'") if nm in e}}
    with open(os.path.join(a.outdir, tag + "_seam.json"), "w") as f:
        json.dump(seam, f)
    with open(os.path.join(a.outdir, tag + "_qc.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(qc[0].keys()) if qc else ["linker"])
        w.writeheader(); w.writerows(qc)
    bad = [r["linker"] for r in qc if not r["ok"]]
    if bad:
        # never mark a block done with a broken linker: keep the evidence, exit non-zero
        os.replace(tmp, out_cif.replace(".cif", "_FAILED.cif"))
        print(f"{tag}: FAILED closure at linkers {bad} -- block NOT marked done")
        sys.exit(2)
    os.replace(tmp, out_cif)                      # last: marks the block done

    d = np.array([j["d_O3_O5"] for j in jstats]) if jstats else np.array([np.nan])
    ang = np.array([j["O3_P_O5"] for j in jstats]) if jstats else np.array([np.nan])
    n_atoms = sum(len(r) for c in gm for r in c)
    bad = sum(not r["ok"] for r in qc)
    print(f"{tag}: wrote {out_cif}  ({n_atoms:,} atoms, {len(qc)} linkers, "
          f"{bad} failed closure) in {time.time()-t0:.0f}s")
    print(f"  junctions: O3'...O5' {np.nanmean(d):.3f} +/- {np.nanstd(d):.3f} A, "
          f"O3'-P-O5' {np.nanmean(ang):.1f} +/- {np.nanstd(ang):.1f} deg")


if __name__ == "__main__":
    main()
