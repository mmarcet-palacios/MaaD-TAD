"""
Option 1: refine the CSynth guide so a physical fibre can thread it.

CSynth folds a bead chain with no finite fibre thickness and no stiffness.
This takes its equilibrium structure and minimises

    E = w_ref * (keep CSynth's non-local distances, i.e. the fold)
      + reach   (bead gap <= fibre axis the section's nucleosomes supply)
      + bend    (guide bend radius >= R_bend)
      + excl    (non-local beads >= d_min apart: 2 x (coil radius + fibre radius))

All lengths are absolute Angstrom: a section's axis demand comes from its
nucleosomes, not from the domain size. --domain only sets the starting scale.

Writes the refined guide in the same format as gap-filled.xyz (Angstrom), and
prints the --domain value to give stage1_frames.py so it is not rescaled.
"""
import argparse
import numpy as np
from scipy.optimize import minimize

import fibre_adaptive as fa
from build_zigzag import canonical_1kx5, place
from stage1_frames import load_dyads, load_guide, ROI_START, BIN, ADV_PER_STEP, FIRST_SECTION

FIBRE_R = 110.0
ADV_NUC = 66.9


def demand(bed, nsec):
    pos, _ = load_dyads(bed)
    nrl = np.diff(pos)
    le, lx = canonical_1kx5()
    cache = {}

    def sl(b):
        b = int(np.clip(b, 10, 200))
        if b not in cache:
            cache[b] = fa.solve_length(b, place, le, lx)[0]
        return cache[b]
    Ls = np.array([sl(int(x) - 147) for x in nrl])
    adv = Ls * ADV_PER_STEP                        # exact b2 steps: advance = L * t_hat.a
    sec_of = np.clip((pos[:-1] - ROI_START) // BIN, 0, nsec - 1)
    return np.array([adv[sec_of == k].sum() for k in range(nsec)])


def energy_factory(x0, need, coil_r, R_bend, w_ref, ref_cut, m_max, c_bend, sep_ev=3,
                   m_min=1.0 / 0.95, ref_min_sep=2):
    n = len(x0)
    I, J = np.triu_indices(n, 2)
    D0 = np.linalg.norm(x0[I] - x0[J], axis=1)
    ref = (np.abs(I - J) >= ref_min_sep) & (D0 < ref_cut)
    ev = np.abs(I - J) >= sep_ev
    d_min = 2 * (coil_r + FIBRE_R)
    bonds = np.arange(n - 1)
    c_hi = need / m_min                       # reachable (m >= m_min)
    c_lo = need / m_max                       # not so short the coil must be huge
    c_lo = np.where(need > 0, c_lo, 0.0)

    def E(flat):
        x = flat.reshape(n, 3)
        g = np.zeros_like(x); e = 0.0

        def pair(ii, jj, d, dE_dd):
            u = (x[ii] - x[jj]) / np.maximum(d, 1e-9)[:, None]
            np.add.at(g, ii, dE_dd[:, None] * u)
            np.add.at(g, jj, -dE_dd[:, None] * u)

        d = np.linalg.norm(x[I] - x[J], axis=1)
        # fold preservation (relative error)
        r = (d[ref] - D0[ref]) / D0[ref]
        e += w_ref * np.sum(r ** 2)
        pair(I[ref], J[ref], d[ref], w_ref * 2 * r / D0[ref])
        # excluded volume
        v = ev & (d < d_min)
        h = (d_min - d[v]) / d_min
        e += np.sum(h ** 2)
        pair(I[v], J[v], d[v], -2 * h / d_min)
        # bond window
        b = np.linalg.norm(x[bonds + 1] - x[bonds], axis=1)
        over = np.maximum(b - c_hi, 0) / 100.0
        under = np.maximum(c_lo - b, 0) / 100.0
        e += 10 * np.sum(over ** 2 + under ** 2)
        pair(bonds + 1, bonds, b, 10 * 2 * (over - under) / 100.0)
        # bend: d(i,i+2) >= 2c cos(c / 2R)
        i2 = np.arange(n - 2)
        c = c_bend                            # frozen per pass (see main)
        tgt = 2 * c * np.cos(np.minimum(c / (2 * R_bend), np.pi / 2))
        d2 = np.linalg.norm(x[i2 + 2] - x[i2], axis=1)
        s = np.maximum(tgt - d2, 0) / 100.0
        e += 10 * np.sum(s ** 2)
        pair(i2 + 2, i2, d2, -10 * 2 * s / 100.0)
        return e, g.ravel()
    return E, d_min


def report(x, need, d_min, R_bend, label):
    n = len(x)
    b = np.linalg.norm(np.diff(x, axis=0), axis=1)
    t = np.diff(x, axis=0) / b[:, None]
    ang = np.arccos(np.clip(np.sum(t[1:] * t[:-1], 1), -1, 1))
    Rb = 0.5 * (b[1:] + b[:-1]) / np.maximum(ang, 1e-9)
    D = np.linalg.norm(x[:, None] - x[None], axis=-1)
    i = np.arange(n)
    clr = D[np.abs(i[:, None] - i[None]) >= 3].min()
    ok = need > 0
    m = need[ok] / b[ok]
    ext = D.max()
    print(f"{label:>9}: extent {ext/10:5.0f} nm | unreachable {int((m<1).sum())} | "
          f"meander med {np.median(m):.2f} max {m.max():.2f} | "
          f"self-approach {clr:4.0f} A (need {d_min:.0f}) | "
          f"bend radius min {Rb.min():4.0f} A, <{R_bend:.0f}: {int((Rb<R_bend).sum())}")
    return ext


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True); ap.add_argument("--guide", required=True)
    ap.add_argument("--gaps", required=True)
    ap.add_argument("--domain", type=float, default=800.0, help="starting scale, nm")
    ap.add_argument("--coil-radius", type=float, default=300.0, help="A")
    ap.add_argument("--bend-radius", type=float, default=700.0, help="min guide bend radius, A")
    ap.add_argument("--w-ref", type=float, default=0.05)
    ap.add_argument("--ref-cut", type=float, default=3000.0,
                    help="bead pairs closer than this (A) keep their CSynth distance")
    ap.add_argument("--m-max", type=float, default=3.0)
    ap.add_argument("--m-min", type=float, default=1.0 / 0.95)
    ap.add_argument("--ref-min-sep", type=int, default=2,
                    help="only keep CSynth distances for beads at least this far apart in sequence")
    ap.add_argument("--no-coil", action="store_true",
                    help="guide becomes the fibre axis itself (m in [1.0, 1.03], run stage 1 with --coil-radius 0); --coil-radius then acts as clearance margin")
    ap.add_argument("--out", default="data/guide_refolded.xyz")
    a = ap.parse_args()

    if a.no_coil:
        a.m_min, a.m_max = 1.0, 1.03
    x0 = load_guide(a.guide, a.gaps, a.domain)
    need = demand(a.bed, len(x0) - 1)
    x = x0.copy()
    for it in range(8):
        # bend target uses the current bond lengths, frozen for this pass so
        # the gradient is exact; re-frozen between passes
        b = np.linalg.norm(np.diff(x, axis=0), axis=1)
        E, d_min = energy_factory(x0, need, a.coil_radius, a.bend_radius,
                                  a.w_ref, a.ref_cut, a.m_max,
                                  0.5 * (b[:-1] + b[1:]),
                                  m_min=a.m_min, ref_min_sep=a.ref_min_sep)
        if it == 0:
            report(x0, need, d_min, a.bend_radius, "CSynth")
        res = minimize(E, x.ravel(), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=20000, maxfun=40000))
        x = res.x.reshape(-1, 3)
    ext = report(x, need, d_min, a.bend_radius, "refolded")
    # how much did the fold change? RMSD after optimal superposition
    A = x0 - x0.mean(0); B = x - x.mean(0)
    U, _, Vt = np.linalg.svd(B.T @ A)
    Rm = (U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt)
    rmsd = np.sqrt(np.mean(np.sum((B @ Rm - A) ** 2, 1)))
    I, J = np.triu_indices(len(x), 3)
    d0 = np.linalg.norm(x0[I] - x0[J], axis=1); d1 = np.linalg.norm(x[I] - x[J], axis=1)
    rho = np.corrcoef(np.argsort(np.argsort(d0)), np.argsort(np.argsort(d1)))[0, 1]
    print(f"fold change: RMSD {rmsd:.0f} A ({rmsd/ (ext/2):.0%} of radius), "
          f"distance rank correlation {rho:.3f}")

    # write, same format as gap-filled.xyz; report the domain value stage 1 needs
    ids = [int(l.split()[0]) for l in open(a.guide) if len(l.split()) >= 4]
    with open(a.out, "w") as f:
        for k, p in zip(ids, x):
            f.write(f"{k}\t{p[0]:.3f}\t{p[1]:.3f}\t{p[2]:.3f}\n")
    gaps = {int(l) for l in open(a.gaps) if l.strip()}
    keep = ~np.isin(np.arange(len(x)) + FIRST_SECTION, list(gaps))
    dom = np.linalg.norm(x[keep][:, None] - x[keep][None], axis=-1).max() / 10
    print(f"wrote {a.out}   ->  give stage1_frames.py  --domain {dom:.1f}")


if __name__ == "__main__":
    main()
