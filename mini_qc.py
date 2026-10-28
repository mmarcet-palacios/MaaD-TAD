"""Stage-4 checks (breaks, far-pair clashes) on a folder of partial blocks."""
import argparse, csv, glob, os, re
import numpy as np, gemmi
from scipy.spatial import cKDTree
from stage2_block import numbering, CORE_NT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", required=True); ap.add_argument("--blocks", required=True)
    a = ap.parse_args()
    F = np.load(a.frames); nrl = F["nrl"].astype(int); n = len(nrl) + 1
    a_core, a_link, b_core, b_link = numbering(nrl)
    NRES = int(nrl.sum()) + CORE_NT
    nucA = np.zeros(NRES + 2, int); nucB = np.zeros(NRES + 2, int)
    for i in range(n):
        nucA[a_core[i]:a_core[i] + CORE_NT] = i; nucB[b_core[i]:b_core[i] + CORE_NT] = i
        if i < n - 1:
            L = nrl[i] - CORE_NT; nucA[a_link[i]:a_link[i] + L] = i; nucB[b_link[i]:b_link[i] + L] = i
    P, Pn, CA, CAn = [], [], [], []
    for f in sorted(glob.glob(os.path.join(a.blocks, "block_*.cif"))):
        if f.endswith("_FAILED.cif"):
            continue
        st = gemmi.read_structure(f)
        for ch in st[0]:
            if ch.name in ("A", "B"):
                nuc = nucA if ch.name == "A" else nucB
                for r in ch:
                    at = r.find_atom("P", "*")
                    if at:
                        P.append(at.pos.tolist()); Pn.append(nuc[r.seqid.num])
            else:
                m = re.match(r"P(\d+)", ch.name)
                for r in ch:
                    at = r.find_atom("CA", "*")
                    if at:
                        CA.append(at.pos.tolist()); CAn.append(int(m.group(1)))
    P, Pn, CA, CAn = map(np.array, (P, Pn, CA, CAn))
    rows = []
    for f in sorted(glob.glob(os.path.join(a.blocks, "block_*_qc.csv"))):
        rows += list(csv.DictReader(open(f)))
    print("LINKERS (backbone):")
    for r in sorted(rows, key=lambda r: int(r["linker"])):
        print(f"  linker {int(r['linker']):5d}: closure {float(r['closure_A']):.0e} A, bend {float(r['bend_max_deg']):5.1f}, "
              f"O3'-P max {float(r['max_O3P_A']):.2f} A {'OK' if r['backbone_ok']=='True' else 'BROKEN'}, "
              f"far obstacles {r['far_obstacles']}, attempts {r['attempts']}")
    tP = cKDTree(P); tC = cKDTree(CA)
    res = {}
    def add(kind, i, j, d):
        k = (kind,) + tuple(sorted((int(i), int(j))))
        res[k] = min(res.get(k, 99), d)
    for i, j in tP.query_pairs(5.0, output_type="ndarray"):
        if abs(Pn[i] - Pn[j]) > 5: add("P-P<5", Pn[i], Pn[j], np.linalg.norm(P[i] - P[j]))
    sd = tC.sparse_distance_matrix(tP, 4.0, output_type="ndarray")
    for r in sd:
        if abs(CAn[r["i"]] - Pn[r["j"]]) > 5: add("CA-P<4", CAn[r["i"]], Pn[r["j"]], r["v"])
    for i, j in tC.query_pairs(3.8, output_type="ndarray"):
        if abs(CAn[i] - CAn[j]) > 5: add("CA-CA<3.8", CAn[i], CAn[j], np.linalg.norm(CA[i] - CA[j]))
    print(f"FAR-PAIR CLASHES among rebuilt pieces: {len(res)}")
    for k, d in sorted(res.items()):
        print(f"  {k[0]:<10} nuc {k[1]} - {k[2]}: {d:.2f} A")


if __name__ == "__main__":
    main()
