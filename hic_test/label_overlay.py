"""label_overlay.py -- labels and scale bar for a PyMOL render of a test model, drawn as an overlay.

PyMOL labels are hidden behind the chain and faded by depth cue, so the labels are drawn on top of the
finished render instead. Anchor points (the middle nucleosome of each gene or line arm) are projected
with PyMOL itself, using the exact camera of the render (get_view), so they land on the right pixels.
Each label is placed near its anchor where the render is emptiest, joined to it by a thin leader line.

Run (needs pymol-open-source, numpy, Pillow), from MaaD-TAD\\hic_test:
    python label_overlay.py IGHC figS_IGHC_B.png view_IGHC.txt
  view_IGHC.txt = the set_view (...) text printed by get_view right after the render
Writes  figS_IGHC_B_labelled.png  (flattened)  and  figS_IGHC_B_labels.svg  (render + editable label layer)
"""
import os, re, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED = (61, 58, 51), (107, 103, 92)
FEATS = {
    'IGHC': ('IGHC_chr14', True, [('IGHA2', 105583731, 105588395), ('IGHE', 105597691, 105601728),
                                  ('IGHG4', 105620506, 105626066), ('IGHG2', 105639559, 105644790),
                                  ('IGHGP', 105664633, 105669843), ('IGHA1', 105703995, 105708665),
                                  ('IGHEP1', 105721794, 105722527), ('IGHG1', 105736343, 105743071)]),
    'PWRN': ('PWRN_chr15', False, [('L4 arm 1', 24196001, 24236000), ('L4 arm 2', 24294001, 24342000)])}
BAR_NM = 200
FONT_PX = 160          # about 6 pt when the panel is printed 3.3 inches wide


def read_view(path):
    txt = open(path).read()
    v = [float(x) for x in re.findall(r'-?\d+\.\d+(?:[eE][-+]?\d+)?', txt.split('set_view')[-1])]
    assert len(v) == 18, f'{path}: expected 18 numbers after set_view, found {len(v)}'
    return v


def project(points, view, W, H):
    """pixel (x, y) of each 3D point, rendered by PyMOL with this view at W x H (one marker at a time)"""
    from pymol import cmd
    cmd.reinitialize()
    cmd.bg_color('black'); cmd.set('orthoscopic', 1); cmd.set('field_of_view', abs(view[17]))
    cmd.set('ray_opaque_background', 1); cmd.set('depth_cue', 0); cmd.set('ambient', 1); cmd.set('specular', 0)
    out = []
    for k, p in enumerate(points):
        cmd.delete('m'); cmd.pseudoatom('m', pos=[float(x) for x in p]); cmd.show_as('spheres', 'm')
        cmd.set('sphere_scale', 1.0); cmd.alter('m', 'vdw=20'); cmd.color('white', 'm')
        cmd.set_view(view)
        f = f'/tmp/_mark{k}.png'; cmd.png(f, width=W, height=H, ray=1)
        A = np.asarray(Image.open(f).convert('L')); ys, xs = np.where(A > 128)
        out.append((xs.mean(), ys.mean()) if len(xs) else (np.nan, np.nan))
        os.remove(f)
    return np.array(out)


def font(italic):
    for f in (['ariali.ttf', 'Arial Italic.ttf', 'DejaVuSans-Oblique.ttf'] if italic else ['arial.ttf', 'Arial.ttf', 'DejaVuSans.ttf']):
        try:
            return ImageFont.truetype(f, FONT_PX)
        except OSError:
            pass
    return ImageFont.load_default()


def place(anchors, labels, fg, fnt, W, H):
    """for each anchor, the label box (around a ring of candidate spots) that covers least of the render
    and no earlier label"""
    taken, boxes = np.zeros_like(fg), []
    draw = ImageDraw.Draw(Image.new('L', (10, 10)))
    for (ax, ay), lab in zip(anchors, labels):
        w, h = draw.textbbox((0, 0), lab, font=fnt)[2:]
        best = None
        for r in (220, 320, 440, 580, 740, 920):
            for ang in np.radians(np.arange(0, 360, 15)):
                cx, cy = ax + r * np.cos(ang), ay - r * np.sin(ang)
                x0, y0 = int(cx - w / 2), int(cy - h / 2)
                if x0 < 20 or y0 < 20 or x0 + w > W - 20 or y0 + h > H - 20:
                    continue
                pad = fg[max(0, y0 - 30):y0 + h + 30, max(0, x0 - 30):x0 + w + 30]
                cost = 4 * pad.mean() + 10 * taken[y0:y0 + h, x0:x0 + w].mean() + r / 6000
                if best is None or cost < best[0]:
                    best = (cost, x0, y0, w, h)
        _, x0, y0, w, h = best
        taken[max(0, y0 - 40):y0 + h + 40, max(0, x0 - 40):x0 + w + 40] = 1
        boxes.append((x0, y0, w, h))
    return boxes


def main():
    tag, png, view_txt = sys.argv[1], sys.argv[2], sys.argv[3]
    name, italic, feats = FEATS[tag]
    view = read_view(view_txt)
    img = Image.open(png).convert('RGB'); W, H = img.size
    F = np.load(os.path.join(HERE, f'{name}_H1', 'frames.npz')); O, pos = F['origins'], F['dyad_pos']
    pts, labels = [], []
    for lab, g0, g1 in feats:
        idx = np.where((pos >= g0) & (pos <= g1))[0]
        pts.append(O[idx[len(idx) // 2]]); labels.append(lab)
    # scale: two points 1 um apart along the camera x axis, through the view centre
    R = np.array(view[:9]).reshape(3, 3).T; c = np.array(view[12:15])
    pts += [c - R[0] * 5000, c + R[0] * 5000]
    xy = project(pts, view, W, H)
    px_per_nm = np.linalg.norm(xy[-1] - xy[-2]) / 1000.0
    anchors = xy[:-2]
    fg = (np.asarray(img.convert('L')) < 240).astype(float)
    fnt = font(italic)
    boxes = place(anchors, labels, fg, fnt, W, H)
    d = ImageDraw.Draw(img)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
           f'<image xlink:href="{os.path.basename(png)}" x="0" y="0" width="{W}" height="{H}"/>', '<g id="labels">']
    for (ax, ay), lab, (x0, y0, w, h) in zip(anchors, labels, boxes):
        lx, ly = min(max(ax, x0), x0 + w), min(max(ay, y0), y0 + h)        # nearest point of the box
        d.line([(ax, ay), (lx, ly)], fill=MUTED, width=8)
        d.ellipse([ax - 14, ay - 14, ax + 14, ay + 14], fill=MUTED)
        d.text((x0, y0), lab, font=fnt, fill=INK)
        svg += [f'<line x1="{ax:.0f}" y1="{ay:.0f}" x2="{lx:.0f}" y2="{ly:.0f}" stroke="#6b675c" stroke-width="8"/>',
                f'<circle cx="{ax:.0f}" cy="{ay:.0f}" r="14" fill="#6b675c"/>',
                f'<text x="{x0}" y="{y0 + h * 0.8:.0f}" font-family="Arial" font-size="{FONT_PX}" '
                f'{"font-style=" + chr(34) + "italic" + chr(34) + " " if italic else ""}fill="#3d3a33">{lab}</text>']
    bar = BAR_NM * px_per_nm; bx, by = W - 250 - bar, H - 380
    d.rectangle([bx, by, bx + bar, by + 22], fill=INK)
    bf = font(False); d.text((bx + bar / 2 - d.textlength(f'{BAR_NM} nm', font=bf) / 2, by + 50), f'{BAR_NM} nm', font=bf, fill=INK)
    svg += ['</g>', '<g id="scale_bar">', f'<rect x="{bx:.0f}" y="{by}" width="{bar:.0f}" height="22" fill="#3d3a33"/>',
            f'<text x="{bx + bar / 2:.0f}" y="{by + 50 + FONT_PX * 0.8:.0f}" text-anchor="middle" font-family="Arial" '
            f'font-size="{FONT_PX}" fill="#3d3a33">{BAR_NM} nm</text>', '</g>', '</svg>']
    stem = os.path.splitext(png)[0]
    img.save(stem + '_labelled.png'); open(stem + '_labels.svg', 'w').write('\n'.join(svg) + '\n')
    print(f'{tag}: {px_per_nm:.3f} px per nm ({BAR_NM} nm = {bar:.0f} px); wrote {stem}_labelled.png, {stem}_labels.svg')
    for lab, (ax, ay) in zip(labels, anchors):
        print(f'  {lab:9s} anchor ({ax:.0f}, {ay:.0f})')


if __name__ == '__main__':
    main()
