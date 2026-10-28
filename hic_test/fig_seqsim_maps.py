"""fig_seqsim_maps.py -- Fig 10A: H1 Micro-C (above the diagonal) against SEQSIM (below), 2 kb bins, for IGHC and PWRN.

Same layout, Micro-C half, feature track and axes as fig_models_maps.py (Fig 10B), so the two panels line up.
SEQSIM half: <name>_SEQSIM.csv (max of direct and inverted, as used by the hybrid folds), log scale fixed at 0.4-100,
             with the SEQSIM colour ramp of Fig 9A.

Run in MaaD-TAD\\hic_test (needs numpy, matplotlib; imports fig_models_maps.py from the same folder):
    python fig_seqsim_maps.py      -> figS_seqsim_A_IGHC.png/.svg, figS_seqsim_A_PWRN.png/.svg
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from fig_models_maps import HERE, BIN, RAMP, INK, MUTED, RULE, GENE, PSEUDO, ARM_A, ARM_C, LOCI, clip_log

SEQ_STOPS = ['#fffaf2', '#fbebd5', '#f8dbb7', '#f2c493', '#e9a466', '#e0833a',
             '#ca6a2d', '#b34f20', '#953a16', '#6f280e', '#4a1606']           # Fig 9A SEQSIM ramp
SEQ_LO, SEQ_HI = 0.4, 100.0


def read_csynth(path):
    rows = [l.rstrip('\n').split(',') for l in open(path, encoding='utf-8')][1:]
    return np.array([[float(v) for v in r[1:]] for r in rows])


def rgba(norm, ramp):
    out = ramp(np.nan_to_num(norm, nan=0.0))
    out[np.isnan(norm)] = matplotlib.colors.to_rgba('#e9e5da')               # discarded bins, as Fig 1B
    return out


def panel(tag):
    name, chrom, start, end, feats = LOCI[tag]
    n = (end - start) // BIN
    z = np.load(os.path.join(HERE, f'{name}_microc.npz'))
    mc = z['H1_bal_2000'].astype(float); disc = np.isnan(z['H1_w_2000'])
    up = clip_log(mc); up[disc, :] = np.nan; up[:, disc] = np.nan
    S = read_csynth(os.path.join(HERE, f'{name}_SEQSIM.csv'))
    assert S.shape == (n, n), S.shape
    lo = np.clip((np.log10(np.clip(S, SEQ_LO, SEQ_HI)) - np.log10(SEQ_LO)) / (np.log10(SEQ_HI) - np.log10(SEQ_LO)), 0, 1)

    seq_ramp = LinearSegmentedColormap.from_list('seqsim', SEQ_STOPS)
    img = np.ones((n, n, 4))
    iu, il = np.triu_indices(n, 1), np.tril_indices(n, -1)
    U, L = rgba(up, RAMP), seq_ramp(lo)
    img[iu] = U[iu]; img[il] = L[il]
    img[np.arange(n), np.arange(n)] = matplotlib.colors.to_rgba('#8e1016')   # diagonal: top of the Micro-C ramp, as Fig 10B

    plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 7, 'svg.fonttype': 'none'})
    W, H = 3.2, 3.45
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.16, 0.11, 0.78, 0.78 * W / H])
    tr = fig.add_axes([0.16, 0.11 + 0.78 * W / H + 0.005, 0.78, 0.09], sharex=ax)
    ext = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]
    ax.imshow(img, extent=ext, interpolation='nearest')
    ax.text(0.97, 0.97, 'H1 Micro-C', transform=ax.transAxes, ha='right', va='top', color=INK)
    ax.text(0.03, 0.03, 'SEQSIM', transform=ax.transAxes, ha='left', va='bottom', color=INK)
    t = np.arange(np.ceil(start / 5e4) * 0.05, end / 1e6 + 1e-9, 0.05)
    t = t[(t > start / 1e6 + 0.01) & (t < end / 1e6 - 0.01)]
    ax.set_xticks(t); ax.set_yticks(t)
    ax.set_xticklabels([f'{v:.2f}' for v in t], color=MUTED); ax.set_yticklabels([f'{v:.2f}' for v in t], color=MUTED)
    ax.set_xlabel(f'{chrom} (Mb)', color=MUTED)
    for s in ax.spines.values():
        s.set_color(RULE); s.set_linewidth(0.6)
    ax.tick_params(colors=RULE, width=0.6, length=2.5)
    tr.set_xlim(start / 1e6, end / 1e6); tr.set_ylim(0, 1); tr.axis('off')
    tr.plot([start / 1e6, end / 1e6], [0.25, 0.25], color=RULE, lw=0.6)
    for k, (lab, g0, g1, kind) in enumerate(feats):
        col = {'gene': GENE, 'pseudo': PSEUDO, 'A': ARM_A, 'C': ARM_C}[kind]
        tr.add_patch(Rectangle((g0 / 1e6, 0.1), (g1 - g0) / 1e6, 0.3, color=col, lw=0))
        tr.text((g0 + g1) / 2e6, 0.5 + 0.3 * (k % 2 if tag == 'IGHC' else 0), lab, ha='center', va='bottom',
                fontsize=6, color=INK, style='italic' if tag == 'IGHC' else 'normal')
    for ext_ in ('png', 'svg'):
        fig.savefig(os.path.join(HERE, f'figS_seqsim_A_{tag}.{ext_}'), dpi=600)
    plt.close(fig)
    print(f'wrote figS_seqsim_A_{tag}.png/.svg  ({n} bins; SEQSIM pairs > 1.1: {int((S[il] > 1.1).sum())})')


if __name__ == '__main__':
    for tag in LOCI:
        panel(tag)
