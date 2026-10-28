"""
Path-level refold with homology pairing and cohesin tethers (layer 1 of 3).

Starts from the current guide (data/guide_refolded.xyz) so everything not
involved in a pairing stays where it is, and adds:

  PARALLEL  every homology block's two arms at CONSTANT spacing -- the spacing
            they have now (block median) -- with tangents aligned (direct) or
            anti-aligned (inverted), register-matched: arm-A offset t faces
            arm-C offset t (direct) or 1-t (inverted).
  TETHER    each cohesin site pulled onto its mirror copy (path axes to
            --tether-dist, default 60 A; the naked DNA and ring finish the job
            at nucleosome level).
  PINCH     option (a): the parallel term fades to zero within --pinch-kb of
            a tether, so arms stay parallel away from it and bend in to meet.
  SPACING   the 520 A non-local minimum stays, except between the two
            neighbourhoods of each tether.

Also kept: bead gaps = fibre axis demand (meander 1.00-1.03), bend radius, and
a weak pull toward the starting path.

    python refold_tether.py            (reads data/homology_restraints.json)
Writes data/guide_tethered.xyz and prints a per-block / per-tether scorecard.
Needs:  pip install autograd
"""
import argparse, json
import autograd.numpy as anp
from autograd import grad
import numpy as np
from scipy.optimize import minimize

from stage1_frames import load_guide, ROI_START, BIN, FIRST_SECTION
from refold_guide import demand

FIBRE_R = 110.0


def u_of(g):
    """continuous bead coordinate of genomic position g (bead k at ROI_START + k*BIN)"""
    return (np.asarray(g, float) - ROI_START) / BIN


def interp_idx(u, nb):
    k = np.clip(np.floor(u).astype(int), 0, nb - 2)
    return k, (u - k)


def build_terms(x0, need, R, a):
    nb = len(x0)
    blocks, tethers = R["blocks"], R["tethers"]
    # tether sites (bead coords)
    T = []
    for t in tethers:
        c = 0.5 * (t["cohesin"][0] + t["cohesin"][1])
        T.append((t["id"], t["block"], u_of(c), u_of(t["mirror"])))
    pinch_u = a.pinch_kb * 1000 / BIN
    # parallel samples: every --sample-kb along each block
    S = []
    for b in blocks:
        A0, A1 = u_of(b["armA"]); C0, C1 = u_of(b["armC"])
        n = max(3, int(round((b["armA"][1] - b["armA"][0]) / (a.sample_kb * 1000))) + 1)
        ts = np.linspace(0, 1, n)
        uA = A0 + ts * (A1 - A0)
        uC = C0 + ts * (C1 - C0) if b["type"] == "direct" else C1 - ts * (C1 - C0)
        w = np.ones(n)
        for tid, tb, us, um in T:          # taper near any tether on either arm
            for uu in (uA, uC):
                dmin = np.minimum(np.abs(uu - us), np.abs(uu - um))
                w = np.minimum(w, np.clip(dmin / pinch_u, 0, 1) ** 2)
        kA, fA = interp_idx(uA, nb); kC, fC = interp_idx(uC, nb)
        PA = x0[kA] * (1 - fA)[:, None] + x0[kA + 1] * fA[:, None]
        PC = x0[kC] * (1 - fC)[:, None] + x0[kC + 1] * fC[:, None]
        D = float(np.median(np.linalg.norm(PA - PC, axis=1)))
        sgn = 1.0 if b["type"] == "direct" else -1.0
        S.append(dict(id=b["id"], type=b["type"], kA=kA, fA=fA, kC=kC, fC=fC, w=w, D=D, sgn=sgn,
                      ident=b["ident"]))
    # tether neighbourhoods exempt from the non-local spacing rule
    exempt = np.zeros((nb, nb), bool)
    rad = int(np.ceil(pinch_u))
    for tid, tb, us, um in T:
        i0, j0 = int(round(us)), int(round(um))
        ii = np.arange(max(0, i0 - rad), min(nb, i0 + rad + 1))
        jj = np.arange(max(0, j0 - rad), min(nb, j0 + rad + 1))
        exempt[np.ix_(ii, jj)] = True; exempt[np.ix_(jj, ii)] = True
    return S, T, exempt


def make_energy(x0, need, S, T, exempt, a):
    nb = len(x0)
    I, J = np.triu_indices(nb, 2)
    D0 = np.linalg.norm(x0[I] - x0[J], axis=1)
    ref = D0 < 3000.0
    evm = (np.abs(I - J) >= 3) & ~exempt[I, J]
    Iev, Jev = I[evm], J[evm]
    c_hi = need / a.m_min; c_lo = need / a.m_max
    Tk = [(interp_idx(np.array([us]), nb), interp_idx(np.array([um]), nb)) for _, _, us, um in T]
    # bend limit per bead: relaxed inside tether neighbourhoods
    Rb = np.full(nb - 2, a.bend_radius)
    for _, _, us, um in T:
        for u in (us, um):
            k = int(round(u)); r = int(np.ceil(a.pinch_kb * 1000 / BIN))
            Rb[max(0, k - r - 1):k + r] = a.pinch_bend_radius

    def P(x, k, f):
        return x[k] * (1 - f)[:, None] + x[k + 1] * f[:, None]

    def tang(x, k):
        v = x[k + 1] - x[k]
        return v / anp.sqrt(anp.sum(v * v, axis=1, keepdims=True))

    def E(flat):
        x = anp.reshape(flat, (nb, 3))
        e = 0.0
        d = anp.sqrt(anp.sum((x[I[ref]] - x[J[ref]]) ** 2, axis=1))
        e = e + a.w_ref * anp.sum(((d - D0[ref]) / D0[ref]) ** 2)
        dv = anp.sqrt(anp.sum((x[Iev] - x[Jev]) ** 2, axis=1) + 1e-9)
        e = e + anp.sum(anp.maximum(a.d_min - dv, 0.0) ** 2) / a.d_min ** 2
        b = anp.sqrt(anp.sum((x[1:] - x[:-1]) ** 2, axis=1))
        t = (x[1:] - x[:-1]) / b[:, None]
        cosang = anp.sum(t[1:] * t[:-1], axis=1)
        # curved (spline) length of each segment ~ chord * (1 + turn^2/24), turn = mean of its
        # two end-vertex angles (theta^2 ~ 2(1-cos)); stage 1 walks the spline, so THIS must
        # equal the fibre's axis demand for genomic registration to hold
        th2 = 2.0 * (1.0 - cosang)
        th2_seg = anp.concatenate([th2[:1], 0.5 * (th2[1:] + th2[:-1]), th2[-1:]])
        arc = b * (1.0 + th2_seg / 24.0)
        e = e + 10 * anp.sum((anp.maximum(arc - c_hi, 0) / 100) ** 2 + (anp.maximum(c_lo - arc, 0) / 100) ** 2)
        cmin = anp.cos(anp.minimum(0.5 * (b[1:] + b[:-1]) / Rb, np.pi / 2))
        e = e + 10 * anp.sum(anp.maximum(cmin - cosang, 0) ** 2) * 100
        for s in S:
            PA, PC = P(x, s["kA"], s["fA"]), P(x, s["kC"], s["fC"])
            dd = anp.sqrt(anp.sum((PA - PC) ** 2, axis=1))
            e = e + a.w_par * anp.sum(s["w"] * ((dd - s["D"]) / s["D"]) ** 2)
            ta, tc = tang(x, s["kA"]), tang(x, s["kC"])
            e = e + a.w_dir * anp.sum(s["w"] * (1 - s["sgn"] * anp.sum(ta * tc, axis=1)) ** 2)
            # register: the rung joining homologous positions must be perpendicular to both
            # arms (like a ladder); otherwise parallel arms can slide past each other
            rung = (PC - PA) / s["D"]
            e = e + a.w_reg * anp.sum(s["w"] * (anp.sum(rung * ta, axis=1) ** 2 + anp.sum(rung * tc, axis=1) ** 2))
        for (k1, f1), (k2, f2) in Tk:
            dt = anp.sqrt(anp.sum((P(x, k1, f1) - P(x, k2, f2)) ** 2) + 1e-9)
            e = e + a.w_teth * ((dt - a.tether_dist) / 100.0) ** 2
        return e
    return E, grad(E)


EXEMPT = None


def scorecard(x, S, T, need, a, label):
    nb = len(x)
    P = lambda k, f: x[k] * (1 - f)[:, None] + x[k + 1] * f[:, None]
    print(f"\n{label}")
    print(f"  {'block':<5}{'type':<9}{'ident':>6}{'target':>8}{'spacing CV':>11}{'dir err':>9}")
    for s in S:
        m = s["w"] > 0.99
        if m.sum() < 2:
            m = s["w"] > 0.5
        dd = np.linalg.norm(P(s["kA"], s["fA"]) - P(s["kC"], s["fC"]), axis=1)[m]
        ta = x[s["kA"] + 1] - x[s["kA"]]; tc = x[s["kC"] + 1] - x[s["kC"]]
        ta /= np.linalg.norm(ta, axis=1)[:, None]; tc /= np.linalg.norm(tc, axis=1)[:, None]
        ang = np.degrees(np.arccos(np.clip(s["sgn"] * np.sum(ta * tc, axis=1), -1, 1)))[m]
        rung = P(s["kC"], s["fC"]) - P(s["kA"], s["fA"])
        slide = np.abs(np.sum(rung * ta, axis=1))[m]           # A along the arm
        print(f"  {s['id']:<5}{s['type']:<9}{s['ident']:6.1f}{s['D']:7.0f}A"
              f"{dd.std() / dd.mean():10.1%}{np.median(ang):8.1f}d{np.median(slide):8.0f}A slide")
    for tid, tb, us, um in T:
        (k1, f1), (k2, f2) = interp_idx(np.array([us]), nb), interp_idx(np.array([um]), nb)
        print(f"  tether {tid}: path-axis distance {np.linalg.norm(P(k1, f1) - P(k2, f2)):.0f} A "
              f"(target {a.tether_dist:.0f})")
    D = np.linalg.norm(x[:, None] - x[None], axis=-1); i = np.arange(nb)
    far = (np.abs(i[:, None] - i[None]) >= 3) & ~EXEMPT
    b = np.linalg.norm(np.diff(x, axis=0), axis=1); ok = need > 0
    print(f"  non-local min spacing (outside tether zones) {D[far].min():.0f} A; meander {np.median(need[ok]/b[ok]):.2f} "
          f"({(need[ok]/b[ok]).min():.2f}-{(need[ok]/b[ok]).max():.2f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default="data/nucleosomes_GRCh38_resolved.bed")
    ap.add_argument("--guide", default="data/guide_refolded.xyz")
    ap.add_argument("--gaps", default="data/gap_sections.txt")
    ap.add_argument("--domain", type=float, default=963.5)
    ap.add_argument("--restraints", default="data/homology_restraints.json")
    ap.add_argument("--out", default="data/guide_tethered.xyz")
    ap.add_argument("--sample-kb", type=float, default=1.0)
    ap.add_argument("--pinch-kb", type=float, default=6.0)
    ap.add_argument("--tether-dist", type=float, default=60.0)
    ap.add_argument("--d-min", type=float, default=520.0)
    ap.add_argument("--bend-radius", type=float, default=1300.0)
    ap.add_argument("--pinch-bend-radius", type=float, default=900.0)
    ap.add_argument("--w-ref", type=float, default=0.002)
    ap.add_argument("--w-par", type=float, default=2.0)
    ap.add_argument("--w-dir", type=float, default=10.0)
    ap.add_argument("--w-teth", type=float, default=5.0)
    ap.add_argument("--w-reg", type=float, default=10.0)
    ap.add_argument("--iters", type=int, default=3)
    ap.add_argument("--m-min", type=float, default=1.0,
                    help="fibre demand / curved path length, lower bound (1.0 = path exactly as long as the fibre)")
    ap.add_argument("--m-max", type=float, default=1.005)
    a = ap.parse_args()

    R = json.load(open(a.restraints))
    x0 = load_guide(a.guide, a.gaps, a.domain)
    need = demand(a.bed, len(x0) - 1)
    S, T, exempt = build_terms(x0, need, R, a)
    global EXEMPT; EXEMPT = exempt
    scorecard(x0, S, T, need, a, "BEFORE (current path)")
    x = x0.copy()
    for it in range(a.iters):
        E, dE = make_energy(x0, need, S, T, exempt, a)
        res = minimize(lambda v: (E(v), dE(v)), x.ravel(), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=4000))
        x = res.x.reshape(-1, 3)
        print(f"  pass {it+1}: E = {res.fun:.3f}", flush=True)
    scorecard(x, S, T, need, a, "AFTER")
    A = x0 - x0.mean(0); B = x - x.mean(0)
    print(f"  shift from current path: RMS {np.sqrt(np.mean(np.sum((A-B)**2,1)))/10:.0f} nm")
    ids = [int(l.split()[0]) for l in open(a.guide) if len(l.split()) >= 4]
    with open(a.out, "w") as f:
        for k, p in zip(ids, x):
            f.write(f"{k}\t{p[0]:.3f}\t{p[1]:.3f}\t{p[2]:.3f}\n")
    gaps = {int(l) for l in open(a.gaps) if l.strip()}
    keep = ~np.isin(np.arange(len(x)) + FIRST_SECTION, list(gaps))
    dom = np.linalg.norm(x[keep][:, None] - x[keep][None], axis=-1).max() / 10
    print(f"wrote {a.out}  ->  stage1_frames.py --domain {dom:.1f}")


if __name__ == "__main__":
    main()
