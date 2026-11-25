"""register_path.py -- register error of each homology block on a 2 kb path (Methods, 'Register error'; Fig 6B).
The path is resampled at 20 points per bead; for every fifth point on arm A, the nearest point on arm C is found and
its genomic offset from the homologous position is the register error; the median over the arm is the block value."""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.abspath(__file__))
LO, BIN = 12771273, 2000
R = json.load(open(os.path.join(ROOT, 'data', 'homology_restraints.json')))


def register(xyz_path):
    X = np.loadtxt(xyz_path)[:, 1:]
    n = len(X); t = np.linspace(0, n - 1, (n - 1) * 20 + 1)
    P = np.column_stack([np.interp(t, np.arange(n), X[:, k]) for k in range(3)])
    g = LO + BIN * t                                     # bead k at ROI_START + k*BIN, as in refold_tether.u_of
    out = {}
    for b in R['blocks']:
        A0, A1 = b['armA']; C0, C1 = b['armC']
        ia = np.where((g >= A0) & (g <= A1))[0][::5]; ic = np.where((g >= C0) & (g <= C1))[0]
        near = g[ic][cKDTree(P[ic]).query(P[ia])[1]]
        s = (g[ia] - A0) / (A1 - A0)
        match = C0 + s * (C1 - C0) if b['type'] == 'direct' else C1 - s * (C1 - C0)
        out[b['id']] = float(np.median(np.abs(near - match)))
    return out


if __name__ == '__main__':
    for p in sys.argv[1:]:
        r = register(p); print(os.path.basename(p), 'median over blocks %.0f bp' % np.median(list(r.values())), {k: round(v) for k, v in r.items()})


def pairing(xyz_path, domain_nm, gaps=os.path.join(ROOT, 'data', 'gap_sections.txt')):
    """Pairing of each block on a path, as refold_tether.scorecard measures it (all three tether zones tapered, so every
    path is sampled at the same points): median arm spacing (nm), spacing CV, direction error (deg), and slide (nm),
    the component of the rung joining homologous positions along the arm (0 = in register)."""
    import types
    import refold_tether as RT
    from stage1_frames import load_guide
    x = load_guide(xyz_path, gaps, domain_nm)
    a = types.SimpleNamespace(pinch_kb=6.0, sample_kb=1.0)
    S, T, _ = RT.build_terms(x, None, R, a)
    P = lambda k, f: x[k] * (1 - f)[:, None] + x[k + 1] * f[:, None]
    out = {}
    for s in S:
        m = s['w'] > 0.99
        if m.sum() < 2:
            m = s['w'] > 0.5
        dd = np.linalg.norm(P(s['kA'], s['fA']) - P(s['kC'], s['fC']), axis=1)[m]
        ta = x[s['kA'] + 1] - x[s['kA']]; tc = x[s['kC'] + 1] - x[s['kC']]
        ta /= np.linalg.norm(ta, axis=1)[:, None]; tc /= np.linalg.norm(tc, axis=1)[:, None]
        ang = np.degrees(np.arccos(np.clip(s['sgn'] * np.sum(ta * tc, axis=1), -1, 1)))[m]
        slide = np.abs(np.sum((P(s['kC'], s['fC']) - P(s['kA'], s['fA'])) * ta, axis=1))[m]
        out[s['id']] = dict(spacing_nm=float(np.median(dd)) / 10, cv=float(dd.std() / dd.mean()),
                            dir_deg=float(np.median(ang)), slide_nm=float(np.median(slide)) / 10)
    return out
