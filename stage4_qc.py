"""
Stage 4: QC of the whole built model -- per block, across seams, and long-range.

    python stage4_qc.py --workers 16          (reads blocks/, frames.npz)

Checks
  1. Completeness: every block present; chain A and B each cover residues
     1..N with no numbering gaps.
  2. Covalent continuity: every O3'(i)-P(i+1) bond on both strands, INCLUDING
     the 67 block seams (chain A's seam phosphates come from block_XX_seam.json).
  3. Phosphodiester geometry: O3'...O5' distance and O3'-P-O5' angle; linkages
     outside 95-115 deg counted (the known outside-strand limitation at high
     bend), split by strand and by measured/placeholder region.
  4. Long-range sterics, the check per-block QC cannot do: one KD-tree over all
     DNA phosphorus and histone C-alpha atoms, reporting close approaches
     between parts of the fibre MORE THAN 5 NUCLEOSOMES APART in sequence:
       P-P   < 10 A   contacts (what the contact map is about; counted, not failed)
       P-P   <  5 A   clash
       CA-P  <  4 A   clash
       CA-CA <  3.8 A clash
Writes qc_summary.txt and qc_contacts.csv (every far pair under 10 A).
"""
import argparse, csv, glob, json, os, re, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.spatial import cKDTree

SEP = 5


def read_block(path):
    import gemmi
    st = gemmi.read_structure(path)
    out = {"A": {}, "B": {}, "ca": [], "ca_nuc": []}
    for ch in st[0]:
        if ch.name in ("A", "B"):
            n = len(ch)
            resi = np.empty(n, np.int64)
            xyz = {k: np.full((n, 3), np.nan) for k in ("P", "O3'", "O5'")}
            for j, r in enumerate(ch):
                resi[j] = r.seqid.num
                for a in r:
                    if a.name in xyz:
                        xyz[a.name][j] = (a.pos.x, a.pos.y, a.pos.z)
            out[ch.name] = {"resi": resi, **xyz}
        else:
            m = re.match(r"P(\d+)", ch.name)
            if not m:
                continue
            ni = int(m.group(1))
            for r in ch:
                a = r.find_atom("CA", "*")
                if a:
                    out["ca"].append((a.pos.x, a.pos.y, a.pos.z)); out["ca_nuc"].append(ni)
    out["ca"] = np.array(out["ca"]).reshape(-1, 3)
    out["ca_nuc"] = np.array(out["ca_nuc"], np.int64)
    return path, out


def angle(p, q, r):
    u, v = p - q, r - q
    c = np.einsum("ij,ij->i", u, v) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--blocks", default="blocks")
    ap.add_argument("--size", type=int, default=50)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    t0 = time.time()
    from stage2_block import numbering, CORE_NT

    F = np.load(a.frames)
    nrl = F["nrl"].astype(int); src = F["source"]; n = len(nrl) + 1
    a_core, a_link, b_core, b_link = numbering(nrl)
    NRES = int(nrl.sum()) + CORE_NT
    nblocks = (n + a.size - 1) // a.size
    files = sorted(glob.glob(os.path.join(a.blocks, "block_??.cif")))
    lines = []

    def say(s=""):
        print(s, flush=True); lines.append(s)

    say(f"STAGE 4 QC  {time.strftime('%Y-%m-%d %H:%M')}")
    say(f"blocks found {len(files)}/{nblocks}; expected residues per strand {NRES:,}")
    failed = glob.glob(os.path.join(a.blocks, "*_FAILED.cif"))
    if failed:
        say(f"  WARNING: failed blocks present: {[os.path.basename(f) for f in failed]}")

    # nucleosome index of every residue (linker i -> nucleosome i) and region
    nucA = np.zeros(NRES + 2, np.int64); nucB = np.zeros(NRES + 2, np.int64)
    for i in range(n):
        nucA[a_core[i]:a_core[i] + CORE_NT] = i; nucB[b_core[i]:b_core[i] + CORE_NT] = i
        if i < n - 1:
            L = nrl[i] - CORE_NT
            nucA[a_link[i]:a_link[i] + L] = i; nucB[b_link[i]:b_link[i] + L] = i
    measured = np.array([s == "real" for s in src])

    # ---- read all blocks (parallel) ---------------------------------------
    S = {c: {k: np.full((NRES + 2, 3), np.nan) for k in ("P", "O3'", "O5'")} for c in "AB"}
    seen = {c: np.zeros(NRES + 2, np.int32) for c in "AB"}
    ca, ca_nuc = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, (path, d) in enumerate(ex.map(read_block, files)):
            for c in "AB":
                if len(d[c]):
                    r = d[c]["resi"]; seen[c][r] += 1
                    for nm in ("P", "O3'", "O5'"):
                        S[c][nm][r] = d[c][nm]
            ca.append(d["ca"]); ca_nuc.append(d["ca_nuc"])
            print(f"  read {os.path.basename(path)} ({k+1}/{len(files)})", end="\r", flush=True)
    print()
    for f in sorted(glob.glob(os.path.join(a.blocks, "block_??_seam.json"))):
        s = json.load(open(f))
        if s and os.path.exists(os.path.join(a.blocks, f"block_{(s['core']) // a.size:02d}.cif")):
            for nm, v in s["atoms"].items():
                if nm in S["A"]:
                    S["A"][nm][s["resi"]] = v
    ca = np.vstack(ca); ca_nuc = np.concatenate(ca_nuc)

    # ---- 1. completeness ---------------------------------------------------
    say("\n1. COMPLETENESS")
    for c in "AB":
        present = seen[c][1:NRES + 1]
        say(f"  chain {c}: {int((present > 0).sum()):,}/{NRES:,} residues present, "
            f"missing {int((present == 0).sum()):,}, duplicated {int((present > 1).sum())}")

    # ---- 2 + 3. backbone ---------------------------------------------------
    say("\n2. COVALENT CONTINUITY (O3'-P, all linkages incl. seams)")
    geo_rows = {}
    problems = []
    for c in "AB":
        O3 = S[c]["O3'"][1:NRES]; P = S[c]["P"][2:NRES + 1]; O5 = S[c]["O5'"][2:NRES + 1]
        ok = ~np.isnan(O3[:, 0]) & ~np.isnan(P[:, 0])
        d = np.full(len(O3), np.nan); d[ok] = np.linalg.norm(O3[ok] - P[ok], axis=1)
        breaks = np.where(~ok | (d > 2.0))[0] + 1                  # residue i whose bond to i+1 is broken
        nucX = nucA if c == "A" else nucB
        for r in np.where(ok & (d > 2.0))[0] + 1:            # both atoms present, bond too long
            problems.append(dict(kind="break", chain=c, resi=int(r), linker=int(nucX[r]), nuc_i="", nuc_j="",
                                 dist=round(float(d[r - 1]), 2) if not np.isnan(d[r - 1]) else ""))
        say(f"  chain {c}: {int(ok.sum()):,}/{NRES-1:,} O3'-P bonds, length {np.nanmean(d):.3f} +/- "
            f"{np.nanstd(d):.3f} A, breaks {len(breaks)}  -> {'PASS' if len(breaks)==0 else 'FAIL'}"
            + (f"  first at resi {breaks[:5].tolist()}" if len(breaks) else ""))
        # the seams explicitly: bonds that cross a block boundary
        blk = np.zeros(NRES + 2, np.int64)
        nuc = nucA if c == "A" else nucB
        blk[:] = nuc // a.size
        cross = np.where(blk[1:NRES] != blk[2:NRES + 1])[0]
        cd = d[cross]
        say(f"           seams crossed: {len(cross)}, O3'-P there {np.nanmin(cd) if len(cd) else float('nan'):.3f}"
            f"-{np.nanmax(cd) if len(cd) else float('nan'):.3f} A, broken {int(np.sum(~(cd < 2.0)))}")
        ok2 = ok & ~np.isnan(O5[:, 0])
        ang = np.full(len(O3), np.nan); ang[ok2] = angle(O3[ok2], P[ok2], O5[ok2])
        dd = np.full(len(O3), np.nan); dd[ok2] = np.linalg.norm(O3[ok2] - O5[ok2], axis=1)
        geo_rows[c] = (ang, dd, nuc[1:NRES])

    say("\n3. PHOSPHODIESTER GEOMETRY (O3'-P-O5' outside 95-115 deg = strained)")
    for c in "AB":
        ang, dd, nuc = geo_rows[c]
        bad = (ang < 95) | (ang > 115)
        m = measured[np.clip(nuc, 0, n - 1)]
        say(f"  chain {c}: O3'...O5' {np.nanmean(dd):.3f} +/- {np.nanstd(dd):.3f} A, "
            f"angle {np.nanmean(ang):.1f} +/- {np.nanstd(ang):.1f} deg; strained "
            f"{int(bad.sum()):,}/{int((~np.isnan(ang)).sum()):,} "
            f"({100*bad.sum()/max(1,(~np.isnan(ang)).sum()):.2f}%) "
            f"[measured {int((bad & m).sum())}, placeholder {int((bad & ~m).sum())}]")

    # ---- 4. long-range sterics --------------------------------------------
    say(f"\n4. LONG-RANGE STERICS (parts > {SEP} nucleosomes apart)")
    Ppts, Pnuc = [], []
    for c, nuc in (("A", nucA), ("B", nucB)):
        p = S[c]["P"][1:NRES + 1]; ok = ~np.isnan(p[:, 0])
        Ppts.append(p[ok]); Pnuc.append(nuc[1:NRES + 1][ok])
    Ppts = np.vstack(Ppts); Pnuc = np.concatenate(Pnuc)
    say(f"  {len(Ppts):,} P atoms, {len(ca):,} histone CA atoms")
    tP = cKDTree(Ppts)
    pp = tP.query_pairs(10.0, output_type="ndarray")
    far = np.abs(Pnuc[pp[:, 0]] - Pnuc[pp[:, 1]]) > SEP
    pp = pp[far]
    dpp = np.linalg.norm(Ppts[pp[:, 0]] - Ppts[pp[:, 1]], axis=1) if len(pp) else np.array([])
    say(f"  P-P  < 10 A (contacts):  {len(pp):,} pairs between "
        f"{len({tuple(sorted((Pnuc[i], Pnuc[j]))) for i, j in pp}):,} nucleosome pairs")
    say(f"  P-P  <  5 A (clash):     {int((dpp < 5).sum()):,}  "
        f"-> {'PASS' if (dpp < 5).sum() == 0 else 'FAIL'}   min {dpp.min() if len(dpp) else float('nan'):.1f} A")
    tC = cKDTree(ca)
    sd = tC.sparse_distance_matrix(tP, 4.0, output_type="ndarray")
    sd = sd[np.abs(ca_nuc[sd["i"]] - Pnuc[sd["j"]]) > SEP] if len(sd) else sd
    ncp = len(sd); mincp = float(sd["v"].min()) if ncp else np.inf
    for i, j in pp[dpp < 5] if len(pp) else []:
        problems.append(dict(kind="P-P", chain="", resi="", linker="", nuc_i=int(Pnuc[i]), nuc_j=int(Pnuc[j]),
                             dist=round(float(np.linalg.norm(Ppts[i] - Ppts[j])), 2)))
    for r in sd:
        problems.append(dict(kind="CA-P", chain="", resi="", linker="", nuc_i=int(ca_nuc[r["i"]]),
                             nuc_j=int(Pnuc[r["j"]]), dist=round(float(r["v"]), 2)))
    say(f"  CA-P < 4 A (clash):      {ncp:,}  -> {'PASS' if ncp == 0 else 'FAIL'}"
        + (f"   min {mincp:.1f} A" if ncp else ""))
    cc = tC.query_pairs(3.8, output_type="ndarray")
    cc = cc[np.abs(ca_nuc[cc[:, 0]] - ca_nuc[cc[:, 1]]) > SEP] if len(cc) else cc
    say(f"  CA-CA < 3.8 A (clash):   {len(cc):,}  -> {'PASS' if len(cc) == 0 else 'FAIL'}")
    for i, j in cc:
        problems.append(dict(kind="CA-CA", chain="", resi="", linker="", nuc_i=int(ca_nuc[i]), nuc_j=int(ca_nuc[j]),
                             dist=round(float(np.linalg.norm(ca[i] - ca[j])), 2)))
    with open("qc_problems.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["kind", "chain", "resi", "linker", "nuc_i", "nuc_j", "dist"])
        w.writeheader(); w.writerows(problems)
    say(f"  -> {len(problems)} problems written to qc_problems.csv (input for relax_frames.py --problems)")

    with open("qc_contacts.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["nuc_i", "nuc_j", "P_P_distance_A", "i_measured", "j_measured"])
        for (i, j), dist in zip(pp, dpp):
            ni, nj = sorted((int(Pnuc[i]), int(Pnuc[j])))
            w.writerow([ni, nj, round(float(dist), 2), bool(measured[ni]), bool(measured[nj])])
    say(f"\nwrote qc_summary.txt, qc_contacts.csv   ({time.time()-t0:.0f}s)")
    open("qc_summary.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
