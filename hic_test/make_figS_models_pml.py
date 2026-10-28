"""make_figS_models_pml.py -- write the PyMOL scripts for panels B and C of the gap-free test-model figure.

Writes, for each locus, figS_<tag>.pml (tag = IGHC, PWRN). Each script loads the finished model
(<name>_H1/blocks_final), colours it and stores two views; it renders nothing.
  view figS_B   whole model, grey; IGHC: IG constant genes green, IG pseudogenes dark grey;
                PWRN: arm 1 of line L4 light blue, arm 2 dark blue (DNA and the histones of nucleosomes in them)
  view figS_C   close-up of 8 consecutive nucleosomes in the middle of IGHG1 (IGHC) or of L4 arm 1 (PWRN),
                cartoon (style it with  @figS_<tag>.pml  then  figS_closeup)
Both loci use the same camera distance in figS_B, so the two whole-model images share one scale.
Views come from the nucleosome centres in frames.npz: the image plane holds the two longest principal axes.

Run in MaaD-TAD\\hic_test:   python make_figS_models_pml.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FOV = 20.0                                   # degrees, orthoscopic (as fig4_stage5.pml)
SIZE = 3000                                  # suggested square render, px
COLS = {'hp_dna': [0.62, 0.60, 0.56], 'hp_hist': [0.84, 0.82, 0.78], 'g_gene': [0.184, 0.561, 0.357],
        'g_pseudo': [0.30, 0.30, 0.28], 'armA': [0.53, 0.78, 0.92], 'armC': [0.12, 0.33, 0.62]}
IGHC = [('IGHA2', 105583731, 105588395, 'g_gene'), ('IGHE', 105597691, 105601728, 'g_gene'),
        ('IGHG4', 105620506, 105626066, 'g_gene'), ('IGHG2', 105639559, 105644790, 'g_gene'),
        ('IGHGP', 105664633, 105669843, 'g_pseudo'), ('IGHA1', 105703995, 105708665, 'g_gene'),
        ('IGHEP1', 105721794, 105722527, 'g_pseudo'), ('IGHG1', 105736343, 105743071, 'g_gene')]
PWRN = [('L4 arm 1', 24196001, 24236000, 'armA'), ('L4 arm 2', 24294001, 24342000, 'armC')]
LOCI = {'IGHC': ('IGHC_chr14', 105_540_000, 105_750_000, IGHC, 'IGHG1'),
        'PWRN': ('PWRN_chr15', 24_130_000, 24_440_000, PWRN, 'L4 arm 1')}


def frame(P):
    """camera axes (rows) = principal axes of P; centre = middle of P's bounding box in that frame"""
    c = P.mean(0)
    U, S, Vt = np.linalg.svd(P - c, full_matrices=False)
    x, y = Vt[0], Vt[1]; R = np.array([x, y, np.cross(x, y)])
    Q = (P - c) @ R.T
    return c + R.T @ ((Q.max(0) + Q.min(0)) / 2), R


def view(R, c, d, front, rear):
    m = [R[0, 0], R[1, 0], R[2, 0], R[0, 1], R[1, 1], R[2, 1], R[0, 2], R[1, 2], R[2, 2]]
    v = m + [0.0, 0.0, -d] + list(c) + [front, rear, -FOV]
    return 'set_view (' + ', '.join(f'{x:.6f}' for x in v) + ')'


def chains(idx):
    return '+'.join(f'P{i}*' for i in idx)


def main():
    data = {}
    for tag, (name, start, end, feats, close) in LOCI.items():
        F = np.load(os.path.join(HERE, f'{name}_H1', 'frames.npz'))
        data[tag] = (F['origins'], F['dyad_pos'])
    # one camera distance for both whole-model views (same scale)
    half = 0.0
    for tag, (O, pos) in data.items():
        c, R = frame(O); Q = (O - c) @ R.T
        half = max(half, np.abs(Q[:, :2]).max())
    H = 2 * half * 1.15                                         # visible height, A
    d = H / (2 * np.tan(np.radians(FOV / 2)))
    px100 = SIZE / H * 1000                                     # px per 100 nm
    for tag, (name, start, end, feats, close) in LOCI.items():
        O, pos = data[tag]; n = len(O); nb = (n + 49) // 50
        c, R = frame(O); r = np.linalg.norm(O - c, axis=1).max() + 300
        L = [f'# figS_{tag}.pml -- supplementary figure, gap-free test model {name} (panels B and C)',
             f'# 1. Fresh PyMOL:  cd MaaD-TAD/hic_test   then   @figS_{tag}.pml',
             f'#    (loads {nb} blocks, a few minutes; renders nothing)',
             '# 2. Panel B:  view figS_B   quick look:  png figS_B_test.png, width=600, height=600, ray=0',
             f'#    final:   png figS_{tag}_B.png, width={SIZE}, height={SIZE}, ray=1',
             f'#    scale (same for both loci): at {SIZE} px, 100 nm = {px100:.1f} px; 500 nm = {5 * px100:.0f} px',
             '#    labels are the object labels_B (disable labels_B to render without them)',
             '# 3. Panel C:  figS_closeup   (cartoon of the close-up nucleosomes only, then view figS_C)',
             f'#    final:   png figS_{tag}_C.png, width={SIZE}, height={SIZE}, ray=1',
             '#    Do not touch the mouse between a view command and its png.',
             'reinitialize', 'set auto_dss, 0', 'set suspend_updates, on', 'python',
             'from pymol import cmd', 'import os, time', 't = time.time()',
             f'for i in range({nb}):',
             f'    cmd.load(os.path.join(r"{name}_H1/blocks_final", "block_%02d.cif" % i), "block_%02d" % i)',
             f'print("{name}: %d atoms, %.0f s" % (cmd.count_atoms("all"), time.time() - t))', 'python end',
             'show_as lines, all', 'bg_color white', 'set ray_opaque_background, 1', 'set depth_cue, 1',
             'set orthoscopic, on', f'set field_of_view, {FOV:g}']
        L += [f'set_color {k}, [{v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f}]' for k, v in COLS.items()]
        L += ['color hp_dna, chain A+B', 'color hp_hist, chain P*']
        lab = []
        for name_, g0, g1, col in feats:
            a0, a1 = g0 - start, g1 - start; b0, b1 = end + 1 - g1, end + 1 - g0
            idx = np.where((pos >= g0) & (pos <= g1))[0]
            L.append(f'# {name_} {g0:,}-{g1:,}: {len(idx)} nucleosomes')
            L.append(f'color {col}, (chain A and resi {a0}-{a1}) or (chain B and resi {b0}-{b1})')
            if len(idx):
                L.append(f'color {col}, chain {chains(idx)}')
                p = O[idx].mean(0) + R[1] * 900                   # label 90 nm "above" in the image
                lab.append(f'pseudoatom labels_B, pos=[{p[0]:.1f}, {p[1]:.1f}, {p[2]:.1f}], label="{name_}"')
        L += lab + ['set label_size, -450', 'set label_color, black', 'set label_font_id, 7', 'hide nonbonded, labels_B',
                    view(R, c, d, d - r, d + r), 'view figS_B, store']
        # close-up: 8 consecutive nucleosomes in the middle of the chosen feature
        g0, g1 = next((f[1], f[2]) for f in feats if f[0] == close)
        mid = np.where((pos >= g0) & (pos <= g1))[0]; k = mid[len(mid) // 2]
        idx = np.arange(k - 4, k + 4)
        cc, RR = frame(O[idx]); Q = (O[idx] - cc) @ RR.T
        h = 1.35 * max(np.ptp(Q[:, 0]), np.ptp(Q[:, 1])) + 250
        dd = h / (2 * np.tan(np.radians(FOV / 2))); slab = np.ptp(Q[:, 2]) / 2 + 160
        gA0, gA1 = int(pos[idx[0]] - 73), int(pos[idx[-1]] + 73)
        L += [f'# close-up: nucleosomes {idx[0]}-{idx[-1]} (dyads {pos[idx[0]]:,}-{pos[idx[-1]]:,}) in {close}',
              f'select closeup, (chain A and resi {gA0 - start}-{gA1 - start}) or (chain B and resi {end + 1 - gA1}-{end + 1 - gA0}) or chain {chains(idx)}',
              view(RR, cc, dd, dd - slab, dd + slab), 'view figS_C, store',
              'python', 'def figS_closeup():',
              '    cmd.hide("everything"); cmd.disable("labels_B")',
              '    cmd.set("cartoon_ring_mode", 3); cmd.set("cartoon_ring_finder", 2)',
              '    cmd.show("cartoon", "closeup"); cmd.view("figS_C", "recall")',
              'cmd.extend("figS_closeup", figS_closeup)', 'python end',
              'deselect', 'view figS_B', 'set suspend_updates, off']
        open(os.path.join(HERE, f'figS_{tag}.pml'), 'w').write('\n'.join(L) + '\n')
        print(f'wrote figS_{tag}.pml: {n} nucleosomes, {nb} blocks; whole view {H / 10:.0f} nm across; '
              f'close-up {h / 10:.0f} nm across (nucleosomes {idx[0]}-{idx[-1]})')


if __name__ == '__main__':
    main()
