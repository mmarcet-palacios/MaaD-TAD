"""fig6b.py -- Fig 6B: register error of each homology block on the CSynth path and on the path of the final model.

Run from the repository root:
    python figures/fig6b.py                          # -> figures/final/fig6B.png (+ register_fig6b.csv)
    python figures/fig6b.py --into Fig6.tif          # also paste the panel into the Fig 6 composite
Register error as in register_path.py (Methods): for every fifth point on arm A (20 points per 2 kb bead), the distance
along arm C between its nearest point and the homologous position; median per block.
"""
import argparse, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from register_path import register, R

ORANGE, BLUE, STEM, GRID, CIRC, TXT = '#d9731f', '#3f6fa8', '#d9d4c7', '#eee9dd', '#b8b3a6', '#5f5c57'
NRL = 194
BOX = (1560, 175, 4500, 1300)          # where panel B's plot sits in Fig6.tif (600 dpi), below its title


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--into'); a = ap.parse_args()
    D = lambda n: os.path.join(ROOT, 'data', n)
    cs, fin = register(D('gap-filled.xyz')), register(D('guide_tethered.xyz'))
    ref = register(D('guide_refolded.xyz'))
    ids = [b['id'] for b in R['blocks']]; typ = {b['id']: b['type'] for b in R['blocks']}
    out = os.path.join(ROOT, 'figures', 'final'); os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'register_fig6b.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['block', 'type', 'csynth_bp', 'refolded_bp', 'final_bp'])
        for i in ids:
            w.writerow([i, typ[i], round(cs[i]), round(ref[i]), round(fin[i])])
    mc, mr, mf = (np.median(list(d.values())) for d in (cs, ref, fin))
    print(f'median register error: CSynth {mc:.0f} bp, refolded {mr:.0f} bp, final {mf:.0f} bp')

    plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 6.5, 'text.color': TXT,
                         'axes.labelcolor': TXT, 'xtick.color': TXT, 'ytick.color': TXT})
    W, H = (BOX[2] - BOX[0]) / 600, (BOX[3] - BOX[1]) / 600
    fig = plt.figure(figsize=(W, H), dpi=600)
    ax = fig.add_axes([0.11, 0.2, 0.86, 0.77])
    x = np.arange(len(ids))
    for k, i in enumerate(ids):
        ax.plot([k, k], [fin[i], cs[i]], color=STEM, lw=1.2, zorder=1)
        ax.scatter(k, cs[i], s=22, facecolor='white', edgecolor=CIRC, lw=1.1, zorder=3)
        ax.scatter(k, fin[i], s=22, color=ORANGE if typ[i] == 'inverted' else BLUE, lw=0, zorder=4)
    ax.axhline(NRL, color='#333333', lw=0.7, ls=(0, (4, 3)), zorder=2)
    ax.set_yscale('log'); ax.set_ylim(2, 12000); ax.set_xlim(-0.5, len(ids) - 0.5)
    ax.set_yticks([10, 100, 1000, 10000]); ax.set_yticklabels(['10', '100', '1 kb', '10 kb'])
    ax.minorticks_off()
    ax.set_xticks(x); ax.set_xticklabels(ids)
    ax.grid(axis='y', color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0, pad=3)
    ax.set_ylabel('register error (bp, log)', labelpad=4)
    k1 = ids.index('D1'); ax.annotate('shared arm', (k1, fin['D1']), xytext=(0, -7), textcoords='offset points',
                                      ha='center', va='top', fontsize=6, color='#222222')
    lg = fig.add_axes([0, 0, 1, 0.09]); lg.set_xlim(0, 1); lg.set_ylim(0, 1); lg.axis('off'); y = 0.45
    lg.scatter(0.1, y, s=22, facecolor='white', edgecolor=CIRC, lw=1.1)
    lg.text(0.112, y, f'path from the map: median {mc:,.0f} bp', va='center')
    lg.scatter([0.43, 0.445], [y, y], s=22, color=[ORANGE, BLUE], lw=0)
    lg.text(0.457, y, f'final model: median {mf:,.0f} bp', va='center')
    lg.plot([0.74, 0.78], [y, y], color='#333333', lw=0.7, ls=(0, (4, 3)))
    lg.text(0.788, y, 'one nucleosome repeat', va='center')
    png = os.path.join(out, 'fig6B.png'); fig.savefig(png, dpi=600, facecolor='white'); plt.close(fig)
    if a.into:
        from PIL import Image
        im = Image.open(a.into).convert('RGB')
        p = Image.open(png).convert('RGB').resize((BOX[2] - BOX[0], BOX[3] - BOX[1]))
        im.paste(p, BOX[:2])
        im.save(a.into, compression='tiff_lzw', dpi=(600, 600))
        print('pasted into', a.into)


if __name__ == '__main__':
    main()
