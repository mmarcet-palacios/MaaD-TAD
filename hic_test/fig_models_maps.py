"""fig_models_maps.py -- supplementary figure, panel A for each gap-free test model.

H1 Micro-C (above the diagonal) against contacts of the full-atom model (below the diagonal), 2 kb bins.
Model contacts: distance^-3 between the mean nucleosome centre of each 2 kb bin (frames.npz of the build).
Both triangles: log scale, clipped at the 5th and 99.3rd percentiles of their own values, Fig 1B colour ramp.
Track above the map: IGHC -- IG constant genes (green) and IG pseudogenes (dark grey);
                     PWRN -- the two arms of line L4 (light and dark blue, as in Fig 10).

Run in MaaD-TAD\\hic_test (needs numpy, matplotlib):
    python fig_models_maps.py      -> figS_models_A_IGHC.png/.svg, figS_models_A_PWRN.png/.svg
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = 2000
RAMP = LinearSegmentedColormap.from_list('fig1b', ['#ffffff', '#f6b3a3', '#e0402b', '#8e1016'])
RAMP.set_bad('#e9e5da')                                   # discarded bins, as Fig 1B
INK, MUTED, RULE = '#3d3a33', '#6b675c', '#8a8578'
GENE, PSEUDO = '#2f8f5b', '#4d4d47'                       # g_pramef and g_pseudo of fig7c.pml
ARM_A, ARM_C = '#87c7eb', '#1f549e'                       # armA / armC of fig10_d3d4.pml

# GRCh38, Ensembl (overlap/region, 14 Oct 2026), all on the minus strand
IGHC_GENES = [('IGHA2', 105583731, 105588395, 'gene'), ('IGHE', 105597691, 105601728, 'gene'),
              ('IGHG4', 105620506, 105626066, 'gene'), ('IGHG2', 105639559, 105644790, 'gene'),
              ('IGHGP', 105664633, 105669843, 'pseudo'), ('IGHA1', 105703995, 105708665, 'gene'),
              ('IGHEP1', 105721794, 105722527, 'pseudo'), ('IGHG1', 105736343, 105743071, 'gene')]
# line L4 of the PWRN SEQSIM map (longest of its 7 lines: 26 bin pairs, SEQSIM > 1.1), 1-based
PWRN_L4 = [('L4 arm 1', 24196001, 24236000, 'A'), ('L4 arm 2', 24294001, 24342000, 'C')]
LOCI = {'IGHC': ('IGHC_chr14', 'chr14', 105_540_000, 105_750_000, IGHC_GENES),
        'PWRN': ('PWRN_chr15', 'chr15', 24_130_000, 24_440_000, PWRN_L4)}


def model_contacts(frames, start, n):
    F = np.load(frames); O, pos = F['origins'], F['dyad_pos']
    k = np.clip((pos - start) // BIN, 0, n - 1)
    X = np.array([O[k == b].mean(0) if (k == b).any() else [np.nan] * 3 for b in range(n)])
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1)) / 10.0        # nm
    with np.errstate(divide='ignore'):
        C = D ** -3.0
    np.fill_diagonal(C, np.nan)
    return C


def clip_log(M):
    """log10, clipped at the 5th/99.3rd percentiles of the non-zero values; zero counts -> bottom of the ramp"""
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.log10(np.where(M > 0, M, np.nan))
    lo, hi = np.nanpercentile(L, [5, 99.3])
    out = np.clip((L - lo) / (hi - lo), 0, 1)
    out[(M == 0) & ~np.isnan(M)] = 0.0
    return out


def panel(tag):
    name, chrom, start, end, feats = LOCI[tag]
    n = (end - start) // BIN
    z = np.load(os.path.join(HERE, f'{name}_microc.npz'))
    mc = z['H1_bal_2000'].astype(float); disc = np.isnan(z['H1_w_2000'])
    mo = model_contacts(os.path.join(HERE, f'{name}_H1', 'frames.npz'), start, n)
    up, lo = clip_log(mc), clip_log(mo)
    up[disc, :] = np.nan; up[:, disc] = np.nan                      # discarded bins: grey (Micro-C half only)
    A = np.full((n, n), np.nan)
    iu, il = np.triu_indices(n, 1), np.tril_indices(n, -1)
    A[iu] = up[iu]; A[il] = lo[il]
    np.fill_diagonal(A, 1.0)

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'svg.fonttype': 'none'})
    W, H = 3.2, 3.45
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.16, 0.11, 0.78, 0.78 * W / H])
    tr = fig.add_axes([0.16, 0.11 + 0.78 * W / H + 0.005, 0.78, 0.09], sharex=ax)
    ext = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]
    ax.imshow(A, cmap=RAMP, vmin=0, vmax=1, extent=ext, interpolation='nearest')
    ax.text(0.97, 0.97, 'H1 Micro-C', transform=ax.transAxes, ha='right', va='top', color=INK)
    ax.text(0.03, 0.03, 'model', transform=ax.transAxes, ha='left', va='bottom', color=INK)
    t = np.arange(np.ceil(start / 5e4) * 0.05, end / 1e6 + 1e-9, 0.05)
    t = t[(t > start / 1e6 + 0.01) & (t < end / 1e6 - 0.01)]
    ax.set_xticks(t); ax.set_yticks(t)
    ax.set_xticklabels([f'{v:.2f}' for v in t], color=MUTED); ax.set_yticklabels([f'{v:.2f}' for v in t], color=MUTED)
    ax.set_xlabel(f'{chrom} (Mb)', color=MUTED)
    for s in ax.spines.values():
        s.set_color(RULE); s.set_linewidth(0.6)
    ax.tick_params(colors=RULE, width=0.6, length=2.5)
    # feature track
    tr.set_xlim(start / 1e6, end / 1e6); tr.set_ylim(0, 1); tr.axis('off')
    tr.plot([start / 1e6, end / 1e6], [0.25, 0.25], color=RULE, lw=0.6)
    for k, (lab, g0, g1, kind) in enumerate(feats):
        col = {'gene': GENE, 'pseudo': PSEUDO, 'A': ARM_A, 'C': ARM_C}[kind]
        tr.add_patch(Rectangle((g0 / 1e6, 0.1), (g1 - g0) / 1e6, 0.3, color=col, lw=0))
        tr.text((g0 + g1) / 2e6, 0.5 + 0.3 * (k % 2 if tag == 'IGHC' else 0), lab, ha='center', va='bottom',
                fontsize=6, color=INK, style='italic' if tag == 'IGHC' else 'normal')
    for ext_ in ('png', 'svg'):
        fig.savefig(os.path.join(HERE, f'figS_models_A_{tag}.{ext_}'), dpi=600)
    plt.close(fig)
    print(f'wrote figS_models_A_{tag}.png/.svg  ({n} bins, {int(disc.sum())} discarded in H1)')


if __name__ == '__main__':
    for tag in LOCI:
        panel(tag)
