"""
Generate a fibre from the step transform fitted to b2, then audit its linkers
with exactly the test that gave the spline model 17.4% OK.

If placement-by-measured-transform is right, feasibility should be near-total
without bending anything, because the geometry is copied from an array that
already exists.
"""
import sys
import numpy as np
import gemmi
from scipy.cluster.hierarchy import fcluster, linkage

sys.path.insert(0, "/tmp/dna_poc/bpstep")
import bpstep as bs
import fibre_fit as ff
import closure as cl
from analyse_b2 import load

KX5 = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'data', '1KX5.cif')


def b2_frames():
    prot, dna = load()
    cents = np.array([v.mean(axis=0) for v in prot.values()])
    n_nuc = len(cents) // 8
    lab = fcluster(linkage(cents, "ward"), n_nuc, criterion="maxclust")
    oct_c = np.array([cents[lab == k].mean(axis=0) for k in range(1, n_nuc + 1)])

    resi, xyz = dna["A"]
    d = np.linalg.norm(xyz[:, None, :] - oct_c[None, :, :], axis=-1)
    near = np.argmin(d, axis=1)
    assign = np.where(d[np.arange(len(xyz)), near] < 60.0, near, -1)

    runs, cur, start = [], assign[0], 0
    for i in range(1, len(assign) + 1):
        if i == len(assign) or assign[i] != cur:
            if cur >= 0 and (i - start) >= 90:
                runs.append((cur, start, i))
            if i < len(assign):
                cur, start = assign[i], i
    return [ff.nucleosome_frame(xyz[s:e], xyz[(s + e) // 2]) for _, s, e in runs]


def canonical_1kx5():
    """1KX5's nucleosome frame in the SAME sign-fixed convention, plus its
    entry(-73) and exit(+73) bp frames expressed in that local frame."""
    st = gemmi.read_structure(KX5)
    st.setup_entities()
    model = st[0]

    wrapped, dyad = [], None
    for r in model["I"]:
        if bs.is_dna(r.name):
            for a in r:
                if a.name == "C1'":
                    wrapped.append([a.pos.x, a.pos.y, a.pos.z])
                    if r.seqid.num == 0:
                        dyad = np.array([a.pos.x, a.pos.y, a.pos.z])
    wrapped = np.array(wrapped)
    nuc = ff.nucleosome_frame(wrapped, dyad)

    frames, _, _ = bs.chain_frames(model, "I", "J", lambda i: -i,
                                   list(range(-73, 74)))
    o_n, R_n = nuc

    def local(fr):
        o, R = fr
        return (R_n.T @ (o - o_n), R_n.T @ R)

    return local(frames[0]), local(frames[-1])


def place(nuc_frame, local_bp):
    o_n, R_n = nuc_frame
    o_l, R_l = local_bp
    return (o_n + R_n @ o_l, R_n @ R_l)


def main():
    src = b2_frames()
    fit = ff.fit_step(src)
    print(f"step fitted to b2 ({len(src)} nucleosomes)")
    print(f"  |t| = {fit['t_norm']:.1f} A   rotation = {fit['rot_deg']:.1f} deg")
    print(f"  residual: {fit['resid_t'].mean():.1f} A / {fit['resid_R'].mean():.1f} deg\n")

    print("SOURCE ARRAY vs FIBRE GENERATED FROM THE FITTED STEP")
    gen = ff.generate(len(src), fit["t"], fit["R"])
    ps, pg = ff.helix_params(src), ff.helix_params(gen)
    print(f"  {'':22s} {'b2 (source)':>13s} {'generated':>12s}")
    for k in ("radius_A", "rise_A", "deg_per_nuc", "nuc_per_turn", "pitch_A", "diameter_nm"):
        print(f"  {k:22s} {ps[k]:13.2f} {pg[k]:12.2f}")

    loc_entry, loc_exit = canonical_1kx5()
    NRL = 165
    linker_bp = NRL - 147
    print(f"\nLINKER FEASIBILITY OF A 60-NUCLEOSOME GENERATED FIBRE")
    print(f"  canonical 1KX5 nucleosomes, NRL {NRL} -> {linker_bp} bp linkers\n")

    fib = ff.generate(60, fit["t"], fit["R"])
    counts, bends, spans = {}, [], []
    for i in range(len(fib) - 1):
        ex = place(fib[i], loc_exit)
        en = place(fib[i + 1], loc_entry)
        reach = float(np.linalg.norm(en[0] - ex[0]))
        contour = linker_bp * cl.B_RISE
        spans.append(reach)
        if reach >= contour:
            counts["UNREACHABLE"] = counts.get("UNREACHABLE", 0) + 1
            continue
        p, err, _ = cl.solve_path(ex, en, linker_bp)
        f = cl.feasibility(ex, en, linker_bp, p)
        if np.linalg.norm(err[:3]) > 1e-3:
            v = "CLOSURE_FAILED"
        elif f["bend_max_deg"] > cl.MAX_BEND_PER_STEP:
            v = "BEND_IMPOSSIBLE"
        elif f["bend_max_deg"] > cl.ELASTIC_BEND_PER_STEP:
            v = "STRAINED"
        else:
            v = "OK"
        counts[v] = counts.get(v, 0) + 1
        bends.append(f["bend_max_deg"])

    spans = np.array(spans)
    print(f"  exit->entry span : {spans.mean():.1f} +/- {spans.std():.1f} A"
          f"   (contour available {linker_bp*cl.B_RISE:.1f} A)")
    n = sum(counts.values())
    print()
    for v, c in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {v:16s} {c:3d} / {n}   ({100*c/n:.0f}%)")
    if bends:
        print(f"\n  bend required : {np.mean(bends):.1f} +/- {np.std(bends):.1f} deg/step"
              f"  (max {np.max(bends):.1f})")
    print(f"\n  spline model, identical test : 17.4% OK")


if __name__ == "__main__":
    main()
