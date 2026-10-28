"""
Local rigid-body relaxation of nucleosome frames.

    python relax_frames.py --frames frames.npz --out frames.npz
    python relax_frames.py --frames frames.npz --out frames.npz --problems qc_problems.csv   (repair mode)

Nucleosomes near any problem move as rigid bodies (3 translations + 3
rotations); everything else is fixed but still counted. Problems found:
  - centre-centre collisions (< --min-clear, > 2 apart)
  - DNA stand-in clashes (dna_proxy.py): DNA axis vs DNA axis < 19 A, DNA axis
    vs histone CA < 11 A, between elements > 2 nucleosomes apart
  - repair mode: linkers whose backbone broke and element pairs that clashed in
    the built model (from stage4_qc.py's qc_problems.csv)
Minimised: the clashes above, plus
  reach    linker exit->entry span <= --max-ext of its B-DNA contour
  turn     extra rotation per step <= max(its value before, 0.6 deg/bp x length,
           capped at --max-step-deg): short linkers are not asked to turn more
  aim      exit/entry aim <= --max-aim
  relief   repair mode: broken linkers pushed toward (aim_in+aim_out)/bp <= 1.2 deg
  stay     small pull back to the original frames
Needs: pip install autograd
"""
import argparse, csv, json, os
import numpy as np
import autograd.numpy as anp
from autograd import grad
from scipy.optimize import minimize
from scipy.spatial import cKDTree

import closure as cl
import fibre_adaptive as fa
import dna_proxy as DP


def check(O, R, nrl, le, lx, st):
    ext, aim, bend, turn = [], [], [], []
    for i in range(len(O) - 1):
        ex_p = R[i] @ lx[0] + O[i]; ex_z = R[i] @ lx[1][:, 2]
        en_p = R[i + 1] @ le[0] + O[i + 1]; en_z = R[i + 1] @ le[1][:, 2]
        d = en_p - ex_p; sp = np.linalg.norm(d)
        a1 = np.degrees(np.arccos(np.clip(ex_z @ d / sp, -1, 1)))
        a2 = np.degrees(np.arccos(np.clip(en_z @ d / sp, -1, 1)))
        ext.append(sp / ((nrl[i] - 147 + 1) * cl.B_RISE)); aim.append(max(a1, a2))
        turn.append((a1 + a2) / (nrl[i] - 147 + 1))
        M = st["R"].T @ R[i].T @ R[i + 1]
        bend.append(np.degrees(np.arccos(np.clip((np.trace(M) - 1) / 2, -1, 1))))
    return np.array(ext), np.array(aim), np.array(bend), np.array(turn)


def centre_pairs(O, d):
    p = cKDTree(O).query_pairs(d, output_type="ndarray")
    return p[np.abs(p[:, 0] - p[:, 1]) > 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--out", default="frames.npz")
    ap.add_argument("--problems", help="qc_problems.csv from stage4_qc.py (repair mode)")
    ap.add_argument("--tethers", default="data/tethers.json")
    ap.add_argument("--min-clear", type=float, default=120.0)
    ap.add_argument("--max-ext", type=float, default=0.93)
    ap.add_argument("--max-step-deg", type=float, default=20.0)
    ap.add_argument("--max-aim", type=float, default=35.0)
    ap.add_argument("--win", type=int, default=6)
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--push", type=float, default=30.0,
                    help="repair mode: flagged nucleosome pairs end up this much further apart (A)")
    ap.add_argument("--hot-margin", type=float, default=12.0,
                    help="extra stand-in clearance (A) around pairs that clashed in the real build")
    a = ap.parse_args()

    F = dict(np.load(a.frames, allow_pickle=True))
    O = F["origins"].copy(); R = F["rotations"].copy(); nrl = F["nrl"].astype(int); n = len(O)
    nb = nrl - 147
    SL = F["step_len"].astype(float)
    tm = DP.local_templates(); AX, CA, le, lx = tm
    st = fa.b2_step()
    TWD, exempt = DP.load_tw(a.tethers)
    ext0, aim0, bend0, turn0 = check(O, R, nrl, le, lx, st)
    cap = np.minimum(a.max_step_deg, np.maximum(bend0, 0.3 * nb))       # per-linker extra-turn cap
    # starting pose deviation of every linker from the ideal b2 step (translation, A)
    SL0 = F["step_len"].astype(float)
    dev0 = np.array([np.linalg.norm(R[i].T @ (O[i + 1] - O[i]) - st["t_hat"] * SL0[i]) for i in range(n - 1)])

    fix_links, fix_pairs = set(), set()
    if a.problems:
        for r in csv.DictReader(open(a.problems)):
            if r["kind"] == "break":
                fix_links.add(int(r["linker"]))
            else:
                fix_pairs.add((int(r["nuc_i"]), int(r["nuc_j"])))
        print(f"repair mode: {len(fix_links)} broken linkers, {len(fix_pairs)} clashing nucleosome pairs from stage 4")
    hot = set()                                   # nucleosomes around real clashes: wider margin
    for i, j in fix_pairs:
        hot |= set(range(i - 2, i + 3)) | set(range(j - 2, j + 3))
    # push targets: each flagged pair (and the nucleosome closing a flagged linker) must end
    # up --push A further apart than it is now; stage 4 found real atoms there, whatever the
    # stand-in says
    PUSH = set()
    for i, j in fix_pairs:
        for x in (i, i + 1):
            for y in (j, j + 1):
                if 0 <= x < n and 0 <= y < n and abs(x - y) > 2:
                    PUSH.add((min(x, y), max(x, y)))
    PUSH = np.array(sorted(PUSH), int).reshape(-1, 2)
    PUSH_T = (np.linalg.norm(O[PUSH[:, 0]] - O[PUSH[:, 1]], axis=1) + a.push) if len(PUSH) else np.array([])

    def report(tag):
        e, am, b, tu = check(O, R, nrl, le, lx, st)
        cp = centre_pairs(O, 110.0); px = DP.clashes(O, R, tm, exempt, TWD)
        fl = sorted(fix_links)
        print(f"{tag}: centre collisions {len(cp)}; DNA stand-in clashes {len(px)}; ext>1 {int((e>1).sum())}; "
              f"aim>35 {int((am>35).sum())}; extra turn max {b.max():.1f} deg"
              + (f"; broken-linker pose dev (A, deg) {[(round(float(np.linalg.norm(R[k].T @ (O[k+1]-O[k]) - st['t_hat']*SL[k])),0), round(float(b[k]),1)) for k in fl]}" if fl else ""), flush=True)
        return cp, px
    cp, px = report("before")
    if a.problems:
        pass

    for rnd in range(a.rounds):
        involved = set(np.unique(centre_pairs(O, a.min_clear)).tolist())
        for ea, eb, _, _ in px:
            involved |= {ea // 2, eb // 2}
            if ea % 2: involved.add(ea // 2 + 1)
            if eb % 2: involved.add(eb // 2 + 1)
        for i, j in fix_pairs:
            involved |= {i, j}
        for k in fix_links:
            involved |= {k, k + 1}
        if not involved:
            break
        mov = np.zeros(n, bool)
        for i in involved:
            mov[max(0, i - a.win):min(n, i + a.win + 1)] = True
        idx = np.where(mov)[0]; m = len(idx)

        # candidate pairs (built now, distances re-evaluated inside E)
        tree = cKDTree(O)
        cpairs = set()
        for i in idx:
            for j in tree.query_ball_point(O[i], a.min_clear + 400):
                if abs(i - j) > 2:
                    cpairs.add((min(i, j), max(i, j)))
        cpairs = np.array(sorted(cpairs))
        P, E, K = DP.build(O, R, tm, TWD)
        # which stand-in point is which: element, kind, local index/param
        loc = []
        for k in range(n):
            loc += [("c", k, v) for v in AX] + [("h", k, v) for v in CA]
            if k < n - 1:
                p0 = R[k] @ lx[0] + O[k]; p3 = R[k + 1] @ le[0] + O[k + 1]
                if k in TWD:
                    W = TWD[k][0]
                    n1 = max(3, int(np.linalg.norm(W - p0) / 10.0) + 1); n2 = max(3, int(np.linalg.norm(p3 - W) / 10.0) + 1)
                    loc += [("a", k, s) for s in np.linspace(0, 1, n1)[1:-1]] + [("w", k, 0.0)] \
                         + [("b", k, s) for s in np.linspace(0, 1, n2)[1:-1]]
                else:
                    span = np.linalg.norm(p3 - p0); nn = max(3, int(span / 10.0) + 1)
                    loc += [("l", k, s) for s in np.linspace(0, 1, nn)[1:-1]]
        assert len(loc) == len(P), "stand-in bookkeeping mismatch"
        elem_mov = np.array([mov[l[1]] or (l[0] in "lab" and mov[min(l[1] + 1, n - 1)]) for l in loc])
        ptree = cKDTree(P)
        pp = ptree.query_pairs(DP.AXIS_MIN + 25.0, output_type="ndarray")
        keep = (np.abs(E[pp[:, 0]] - E[pp[:, 1]]) > 4) & ~((K[pp[:, 0]] == 1) & (K[pp[:, 1]] == 1)) \
            & (elem_mov[pp[:, 0]] | elem_mov[pp[:, 1]])
        exk = {(min(2 * x + 1, 2 * y + 1), max(2 * x + 1, 2 * y + 1)) for x, y in exempt}
        pp = np.array([q for q in pp[keep] if (min(E[q[0]], E[q[1]]), max(E[q[0]], E[q[1]])) not in exk]).reshape(-1, 2)
        thr = np.where((K[pp[:, 0]] == 1) | (K[pp[:, 1]] == 1), DP.AXIS_CA_MIN, DP.AXIS_MIN) if len(pp) else np.array([])
        if len(pp) and hot:
            nod = lambda q: np.array([loc[x][1] for x in q])
            hh = np.isin(nod(pp[:, 0]), list(hot)) & np.isin(nod(pp[:, 1]), list(hot))
            thr = thr + a.hot_margin * hh

        # point evaluators: split the point list into core/histone refs and linker refs
        def refs(ids):
            kind = np.array([loc[i][0] for i in ids])
            node = np.array([loc[i][1] for i in ids])
            vec = np.array([loc[i][2] if loc[i][0] in "ch" else np.zeros(3) for i in ids]).reshape(-1, 3)
            s = np.array([loc[i][2] if loc[i][0] in "lab" else 0.0 for i in ids])
            return kind, node, vec, s
        RA, RB = (refs(pp[:, 0]), refs(pp[:, 1])) if len(pp) else (None, None)

        O0, R0 = O.copy(), R.copy()
        Sel = np.zeros((n, m)); Sel[idx, np.arange(m)] = 1.0
        keepm = 1.0 - mov.astype(float)
        PI, PJ = cpairs[:, 0], cpairs[:, 1]
        use = mov[PI] | mov[PJ]; PI, PJ = PI[use], PJ[use]
        LI = np.array([i for i in range(n - 1) if mov[i] or mov[i + 1]]); LJ = LI + 1
        cont = (nrl[LI] - 147 + 1) * cl.B_RISE
        ca_ = np.cos(np.radians(a.max_aim)); capL = np.cos(np.radians(cap[LI]))
        FL = np.array([k for k in sorted(fix_links) if mov[k]], int)
        # guard: short linkers (< 40 bp) must not drift further from the ideal b2 step
        SH = np.array([k for k in LI if nb[k] < 40 and k not in fix_links], int)
        SHlim = np.maximum(dev0[SH], 10.0) if len(SH) else np.array([])
        RS = st["R"]

        def batch_rot(W):
            th = anp.sqrt(anp.sum(W * W, axis=1) + 1e-12); k = W / th[:, None]; z = anp.zeros(len(W))
            K_ = anp.stack([anp.stack([z, -k[:, 2], k[:, 1]], axis=1), anp.stack([k[:, 2], z, -k[:, 0]], axis=1),
                            anp.stack([-k[:, 1], k[:, 0], z], axis=1)], axis=1)
            return anp.eye(3)[None] + anp.sin(th)[:, None, None] * K_ + (1 - anp.cos(th))[:, None, None] * anp.einsum("mij,mjk->mik", K_, K_)

        def frames_of(v):
            T = anp.reshape(v[:3 * m], (m, 3)); W = anp.reshape(v[3 * m:], (m, 3))
            Of = O0 + anp.dot(Sel, T)
            Rm = anp.einsum("mij,mjk->mik", batch_rot(W), R0[idx])
            return Of, R0 * keepm[:, None, None] + anp.einsum("nm,mij->nij", Sel, Rm)

        def points(Of, Rf, rf):
            kind, node, vec, s = rf
            pc = anp.einsum("nij,nj->ni", Rf[node], vec) + Of[node]
            nx = np.minimum(node + 1, n - 1)
            p0 = anp.einsum("nij,j->ni", Rf[node], lx[0]) + Of[node]; t0 = anp.einsum("nij,j->ni", Rf[node], lx[1][:, 2])
            p3 = anp.einsum("nij,j->ni", Rf[nx], le[0]) + Of[nx]; t3 = anp.einsum("nij,j->ni", Rf[nx], le[1][:, 2])
            h = anp.sqrt(anp.sum((p3 - p0) ** 2, axis=1, keepdims=True) + 1e-9) / 3.0
            ss = s[:, None]
            pl = ((1 - ss) ** 3) * p0 + 3 * ((1 - ss) ** 2) * ss * (p0 + t0 * h) + 3 * (1 - ss) * ss ** 2 * (p3 - t3 * h) + ss ** 3 * p3
            isl = (kind == "l")[:, None].astype(float)
            out = pc * (1 - isl) + pl * isl
            if TWD:
                # tether halves: exit -> waypoint (first half), waypoint -> entry (second half)
                Wp = np.array([TWD[k][0] if k in TWD else np.zeros(3) for k in node])
                Wd = np.array([TWD[k][1] if k in TWD else np.zeros(3) for k in node])
                ha = anp.sqrt(anp.sum((Wp - p0) ** 2, axis=1, keepdims=True) + 1e-9) / 3.0
                hb = anp.sqrt(anp.sum((p3 - Wp) ** 2, axis=1, keepdims=True) + 1e-9) / 3.0
                pa = ((1 - ss) ** 3) * p0 + 3 * ((1 - ss) ** 2) * ss * (p0 + t0 * ha) + 3 * (1 - ss) * ss ** 2 * (Wp - Wd * ha) + ss ** 3 * Wp
                pb = ((1 - ss) ** 3) * Wp + 3 * ((1 - ss) ** 2) * ss * (Wp + Wd * hb) + 3 * (1 - ss) * ss ** 2 * (p3 - t3 * hb) + ss ** 3 * p3
                ia = (kind == "a")[:, None].astype(float); ib = (kind == "b")[:, None].astype(float)
                iw = (kind == "w")[:, None].astype(float)
                out = out * (1 - ia - ib - iw) + pa * ia + pb * ib + Wp * iw
            return out

        def E(v):
            Of, Rf = frames_of(v)
            e = 1e-6 * anp.sum(v[:3 * m] ** 2) + 1e-2 * anp.sum(v[3 * m:] ** 2)
            d = anp.sqrt(anp.sum((Of[PI] - Of[PJ]) ** 2, axis=1) + 1e-9)
            e = e + anp.sum(anp.maximum(a.min_clear - d, 0.0) ** 2) / 100.0
            if RA is not None:
                dd = anp.sqrt(anp.sum((points(Of, Rf, RA) - points(Of, Rf, RB)) ** 2, axis=1) + 1e-9)
                e = e + anp.sum(anp.maximum(thr - dd, 0.0) ** 2)
            Ri, Rj = Rf[LI], Rf[LJ]
            exp = anp.einsum("nij,j->ni", Ri, lx[0]) + Of[LI]; exz = anp.einsum("nij,j->ni", Ri, lx[1][:, 2])
            enp = anp.einsum("nij,j->ni", Rj, le[0]) + Of[LJ]; enz = anp.einsum("nij,j->ni", Rj, le[1][:, 2])
            dv = enp - exp; sp = anp.sqrt(anp.sum(dv ** 2, axis=1) + 1e-9)
            e = e + 100.0 * anp.sum(anp.maximum(sp / cont - a.max_ext, 0.0) ** 2)
            c1 = anp.sum(exz * dv, axis=1) / sp; c2 = anp.sum(enz * dv, axis=1) / sp
            e = e + 10.0 * anp.sum(anp.maximum(ca_ - c1, 0) ** 2 + anp.maximum(ca_ - c2, 0) ** 2)
            M = anp.einsum("ji,njk->nik", RS, anp.einsum("nji,njk->nik", Ri, Rj))
            tr = M[:, 0, 0] + M[:, 1, 1] + M[:, 2, 2]
            e = e + 10.0 * anp.sum(anp.maximum(capL - (tr - 1) / 2, 0) ** 2)
            if len(PUSH):
                dpu = anp.sqrt(anp.sum((Of[PUSH[:, 0]] - Of[PUSH[:, 1]]) ** 2, axis=1) + 1e-9)
                e = e + 0.05 * anp.sum(anp.maximum(PUSH_T - dpu, 0.0) ** 2)
            if len(SH):
                rs = anp.einsum("nji,nj->ni", Rf[SH], Of[SH + 1] - Of[SH])
                dsh = anp.sqrt(anp.sum((rs - st["t_hat"][None] * SL[SH][:, None]) ** 2, axis=1) + 1e-9)
                e = e + 0.5 * anp.sum(anp.maximum(dsh - SHlim, 0.0) ** 2)
            if len(FL):
                # broken linker: pull the two nucleosomes back toward the ideal b2 step
                # (the validated geometry every linker can close)
                Rf_i, Rf_j = Rf[FL], Rf[FL + 1]
                rel_t = anp.einsum("nji,nj->ni", Rf_i, Of[FL + 1] - Of[FL])
                e = e + 0.05 * anp.sum((rel_t - st["t_hat"][None] * SL[FL][:, None]) ** 2)
                Mf = anp.einsum("ji,njk->nik", RS, anp.einsum("nji,njk->nik", Rf_i, Rf_j))
                trf = Mf[:, 0, 0] + Mf[:, 1, 1] + Mf[:, 2, 2]
                e = e + 200.0 * anp.sum(1.0 - (trf - 1) / 2)
            return e

        dE = grad(E)
        res = minimize(lambda v: (E(v), dE(v)), np.zeros(6 * m), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=3000))
        Of, Rf = frames_of(res.x)
        O, R = np.asarray(Of), np.asarray(Rf)
        cp, px = report(f"round {rnd+1} (moved {m}, largest shift {np.abs(res.x[:3*m]).reshape(-1,3).max():.0f} A)")
        met = (not len(PUSH)) or bool(np.all(np.linalg.norm(O[PUSH[:, 0]] - O[PUSH[:, 1]], axis=1) >= PUSH_T - 1.0))
        if not len(cp) and not px and met and rnd >= 1:
            break
        pass

    F["origins"], F["rotations"] = O, R
    np.savez_compressed(a.out, **F)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
