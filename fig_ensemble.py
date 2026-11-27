"""fig_ensemble.py -- Fig 11: the ensemble of MaaD-TAD paths.
Reads ensemble/ (k = 0.007, calibrated) and ensemble_k04/ (k = 0.04, the strength of the CSynth path), each analysed by
analyse_ensemble.py, and writes figures/final/Fig11_ensemble.png and .tif."""
import json, os, glob, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.argv = sys.argv[:1]
import analyse_ensemble as A

BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#1f1f1e', '#73726c', '#d9d8d2'
POS = LinearSegmentedColormap.from_list('pos', ['#86b6ef', '#2a78d6', '#104281', '#0d366b'])
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
# GraphPad Prism style: black axes, outward ticks, no gridlines, Arial
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 8, 'axes.edgecolor': 'black',
                     'axes.labelcolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black',
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.9,
                     'xtick.major.width': 0.9, 'ytick.major.width': 0.9, 'xtick.minor.width': 0.6, 'ytick.minor.width': 0.6,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'xtick.major.size': 3.5, 'ytick.major.size': 3.5,
                     'ytick.minor.size': 2, 'axes.labelweight': 'bold', 'legend.fontsize': 6.5})


def load(d):
    S = json.load(open(os.path.join(ROOT, d, 'summary.json'))) if os.path.exists(os.path.join(ROOT, d, 'summary.json')) else None
    M = np.load(os.path.join(ROOT, d, 'maps.npz'), allow_pickle=True) if S else None
    rows = [l.strip().split(',') for l in open(os.path.join(ROOT, d, 'members.csv'))] if S else None
    blocks = [l.strip().split(',') for l in open(os.path.join(ROOT, d, 'blocks.csv'))][1:] if S else None
    return S, M, rows, blocks


S1, M1, R1, B1 = load('ensemble')
S2, M2, R2, B2 = load('ensemble_k04')
col = {k: i for i, k in enumerate(R1[0])}
mem1 = [r for r in R1[1:] if r[0] != 'published']; pub = [r for r in R1[1:] if r[0] == 'published'][0]
mem2 = [r for r in R2[1:] if r[0] != 'published'] if R2 else []


def beads(d, frames):
    F = np.load(os.path.join(ROOT, d, frames))
    O = F['origins']; dy = F['dyad_pos'].astype(float)
    g = A.LO + 1000 + 2000 * np.arange(A.N)
    return np.column_stack([np.interp(g, dy, O[:, k]) for k in range(3)]) / 10.0   # nm


def proj(X):
    c = X - X.mean(0); _, _, Vt = np.linalg.svd(c, full_matrices=False)
    return c @ Vt[:2].T


fig = plt.figure(figsize=(7.5, 8.4), dpi=300)
gs = fig.add_gridspec(3, 4, height_ratios=[0.5, 1, 1], hspace=0.42, wspace=0.62, left=0.09, right=0.985, top=0.965, bottom=0.05)

# ---- A: seven paths, same scale
show = [('published', os.path.join('ensemble', '_published'), 'frames_stage1.npz')] + \
       [(r[0], os.path.join('ensemble', r[0]), r[1]) for r in mem1[:5]] + \
       ([('k = 0.04 member', os.path.join('ensemble_k04', mem2[0][0]), mem2[0][1])] if mem2 else [])
sub = gs[0, :].subgridspec(1, len(show), wspace=0.08)
P = [proj(beads(d, f)) for _, d, f in show]
lim = max(np.abs(p).max() for p in P) * 1.05
for i, ((name, _, _), p) in enumerate(zip(show, P)):
    ax = fig.add_subplot(sub[0, i])
    seg = np.stack([p[:-1], p[1:]], 1); real = A.real[:-1] & A.real[1:]
    lc = LineCollection(seg[real], cmap=POS, linewidths=1.0, capstyle='round'); lc.set_array(np.arange(A.N - 1)[real]); ax.add_collection(lc)
    ax.add_collection(LineCollection(seg[~real], colors=GRID, linewidths=0.8, linestyles=(0, (1, 1))))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal'); ax.axis('off')
    tt = name if not name.startswith('m0') else f'member {int(name[1:]) + 1}'
    ax.set_title(tt if name != 'published' else 'full-atom model', fontsize=7, color=INK, pad=2)
    if i == 0:
        ax.plot([-lim, -lim + 500], [-lim * 0.98] * 2, color=INK, lw=1.2, solid_capstyle='butt')
        ax.text(-lim + 250, -lim * 0.9, '500 nm', ha='center', va='bottom', fontsize=6.5, color=INK)
fig.text(0.025, 0.972, 'A', fontsize=10, weight='bold')
cax = fig.add_axes([0.80, 0.815, 0.15, 0.007])
cb = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=POS, norm=matplotlib.colors.Normalize(A.LO / 1e6, (A.LO + 2000 * A.N) / 1e6)), cax=cax, orientation='horizontal')
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6, length=2); cb.set_label('GRCh38 chr1 (Mb); dotted, assembly gap', fontsize=6, color=INK, labelpad=1, weight='normal')


def panel(ax, letter):
    p = ax.get_position(); fig.text(p.x0 - 0.065, p.y1 + 0.012, letter, fontsize=10, weight='bold', va='bottom')


def strip(ax, x, v, c, w=0.28, s=9, alpha=0.85, rng=np.random.default_rng(1)):
    v = np.asarray(v, float); ax.scatter(x + rng.uniform(-w, w, len(v)), v, s=s, color=c, alpha=alpha, lw=0.3, edgecolor='white', zorder=3)
    if len(v):
        ax.plot([x - w - 0.05, x + w + 0.05], [np.median(v)] * 2, color=INK, lw=1.2, zorder=4)


# ---- B: path-to-path similarity
ax = fig.add_subplot(gs[1, 0]); panel(ax, 'B')
rho1 = M1['rho']; n1 = len(mem1)
mm1 = rho1[:n1, :n1][np.triu_indices(n1, 1)]; mp1 = rho1[-1, :n1]
data = [(mm1, BLUE, 'between\nmembers'), (mp1, BLUE, 'with\nmodel')]
if M2 is not None:
    rho2 = M2['rho']; n2 = len(mem2)
    data.append((rho2[:n2, :n2][np.triu_indices(n2, 1)], ORANGE, 'k = 0.04'))
for i, (v, c, lab) in enumerate(data):
    vp = ax.violinplot([v], [i], widths=0.7, showextrema=False)
    for b in vp['bodies']:
        b.set_facecolor(c); b.set_alpha(0.35); b.set_edgecolor('none')
    ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=INK, lw=1.2)
ax.set_xticks(range(len(data))); ax.set_xticklabels([d[2] for d in data], fontsize=5.6); ax.set_xlim(-0.55, len(data) - 0.45)
ax.set_ylabel('path similarity (Spearman ρ)'); ax.set_ylim(-0.1, 0.6); ax.axhline(0, color=MUTED, lw=0.5, ls=':')

# ---- C: homologous-pair enrichment by stage
ax = fig.add_subplot(gs[1, 1:4]); panel(ax, 'C')
stages = [('enr_fold', 'fold'), ('enr_refolded', 'refolded'), ('enr_tethered', 'paired +\ntethered'), ('enr_nucleosomes', 'nucleosomes')]
for i, (k, lab) in enumerate(stages):
    strip(ax, i - 0.2, [float(r[col[k]]) for r in mem1], BLUE, w=0.14, s=7)
    if mem2:
        strip(ax, i + 0.2, [float(r[col[k]]) for r in mem2], ORANGE, w=0.14, s=7)
    ax.scatter([i], [float(pub[col[k]])], marker='D', s=22, color=INK, zorder=5)
ax.axhline(3.0, color=BLUE, lw=0.8, ls='--'); ax.text(3.42, 3.0, ' Micro-C\n calibration', va='center', fontsize=6, color=MUTED)
ax.set_xticks(range(4)); ax.set_xticklabels([s[1] for s in stages], fontsize=7); ax.set_xlim(-0.5, 3.75)
ax.set_ylabel('homologous-pair enrichment'); ax.set_yscale('log'); ax.set_yticks([1, 2, 3, 5, 10]); ax.set_yticklabels(['1', '2', '3', '5', '10'])
h = [plt.Line2D([], [], marker='o', ls='', color=BLUE, ms=4, label='k = 0.007 (calibrated)'),
     plt.Line2D([], [], marker='o', ls='', color=ORANGE, ms=4, label='k = 0.04'),
     plt.Line2D([], [], marker='D', ls='', color=INK, ms=4, label='full-atom model')]
ax.legend(handles=h, fontsize=6, frameon=False, loc='upper left', ncol=3, bbox_to_anchor=(0, 1.13), handletextpad=0.2, columnspacing=0.8)

# ---- E: register of each block before and after pairing (slide along the arms)
ax = fig.add_subplot(gs[2, 0:3]); panel(ax, 'D')
PB1 = [l.strip().split(',') for l in open(os.path.join(ROOT, 'ensemble', 'pairing_blocks.csv'))][1:]
order = [b['id'] for b in A.R_ALL['blocks']]
for i, bid in enumerate(order):
    pre = [float(r[6]) for r in PB1 if r[2] == bid and r[1] == 'refolded' and r[0] != 'published']
    post = [float(r[6]) for r in PB1 if r[2] == bid and r[1] == 'tethered' and r[0] != 'published']
    rng = np.random.default_rng(i)
    ax.scatter(i - 0.18 + rng.uniform(-0.12, 0.12, len(pre)), pre, s=6, facecolor='none', edgecolor=BLUE, lw=0.5, alpha=0.8, zorder=3)
    ax.scatter(i + 0.18 + rng.uniform(-0.12, 0.12, len(post)), post, s=6, color=BLUE, lw=0, alpha=0.85, zorder=3)
    for x, v in ((i - 0.18, pre), (i + 0.18, post)):
        ax.plot([x - 0.15, x + 0.15], [np.median(v)] * 2, color=INK, lw=1.1, zorder=4)
    pv = [float(r[6]) for r in PB1 if r[2] == bid and r[1] == 'tethered' and r[0] == 'published']
    ax.scatter([i + 0.18], pv, marker='D', s=14, color=INK, zorder=5)
ax.set_yscale('log'); ax.set_ylim(0.1, 1000)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=6.5); ax.set_xlim(-0.6, len(order) - 0.4)
ax.set_ylabel('slide along the arms (nm)')
h = [plt.Line2D([], [], marker='o', ls='', mfc='none', mec=BLUE, ms=4, label='refolded, before pairing'),
     plt.Line2D([], [], marker='o', ls='', color=BLUE, ms=4, label='after pairing and tethers'),
     plt.Line2D([], [], marker='D', ls='', color=INK, ms=4, label='full-atom model')]
ax.legend(handles=h, fontsize=6, frameon=False, loc='upper left', ncol=3, bbox_to_anchor=(0, 1.12), handletextpad=0.2, columnspacing=0.8)

# ---- F: tethers on / off
ax = fig.add_subplot(gs[2, 3]); panel(ax, 'E')
for i, t in enumerate(A.TETH):
    on = [float(r[col[t + '_nm']]) for r in mem1 if t in r[col['tethers']].split('+')]
    off = [float(r[col[t + '_nm']]) for r in mem1 if t not in r[col['tethers']].split('+')]
    strip(ax, i - 0.18, on, BLUE, w=0.12, s=6)
    strip(ax, i + 0.18, off, '#9ec5f4', w=0.12, s=6)
    ax.scatter([i - 0.18], [float(pub[col[t + '_nm']])], marker='D', s=14, color=INK, zorder=5)
ax.set_xticks(range(3)); ax.set_xticklabels([t.replace('T_', '') for t in A.TETH], fontsize=6.5)
ax.set_ylabel('cohesin site to mirror copy (nm)'); ax.set_yscale('log')
h = [plt.Line2D([], [], marker='o', ls='', color=BLUE, ms=4, label='tethered'), plt.Line2D([], [], marker='o', ls='', color='#9ec5f4', ms=4, label='not tethered'), plt.Line2D([], [], marker='D', ls='', color=INK, ms=4, label='model')]
ax.legend(handles=h, frameon=False, loc='lower left', bbox_to_anchor=(-0.12, 0.99), ncol=2, handletextpad=0.05, columnspacing=0.5, fontsize=6)

out = os.path.join(ROOT, 'figures', 'final'); os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, 'Fig11_ensemble.png'), dpi=300)
fig.savefig(os.path.join(out, 'Fig11_ensemble.tif'), dpi=300, pil_kwargs={'compression': 'tiff_lzw'})

# ---- S4 Fig: domain size (radius of gyration of the nucleosome centres)
fig = plt.figure(figsize=(2.6, 2.6), dpi=300)
ax = fig.add_axes([0.3, 0.16, 0.65, 0.78])
panel = lambda ax, letter: None
# ---- D: domain size
strip(ax, 0, [float(r[col['rg_nm']]) for r in mem1], BLUE)
if mem2:
    strip(ax, 1, [float(r[col['rg_nm']]) for r in mem2], ORANGE)
ax.scatter([0.5], [float(pub[col['rg_nm']])], marker='D', s=22, color=INK, zorder=5)
ax.set_xticks([0, 1]); ax.set_xticklabels(['k = 0.007', 'k = 0.04'], fontsize=6.5); ax.set_xlim(-0.6, 1.6)
ax.set_ylabel('radius of gyration (nm)'); ax.set_ylim(0, None)

fig.savefig(os.path.join(out, 'S4_Fig_rg.png'), dpi=300)
fig.savefig(os.path.join(out, 'S4_Fig_rg.tif'), dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
print('wrote', out)
