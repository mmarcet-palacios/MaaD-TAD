"""add_scalebar.py -- draw a scale bar on a render, in the emptiest corner.

    python add_scalebar.py figS_IGHC_C.png 69.33 10
      (image, pixels per nm in that image, bar length in nm)
Writes <image>_bar.png (flattened) and <image>_bar.svg (render + editable bar layer).
Pixels per nm come from the camera of the render (label_overlay.py prints it for whole-model views; for the
IGHC close-up it was measured by re-rendering the close-up with the same view and registering the two images).
"""
import os, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

INK = (61, 58, 51)


def font(px):
    for f in ('arial.ttf', 'Arial.ttf', 'DejaVuSans.ttf'):
        try:
            return ImageFont.truetype(f, px)
        except OSError:
            pass
    return ImageFont.load_default()


def main():
    png, ppn, nm = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    img = Image.open(png).convert('RGB'); W, H = img.size
    fg = np.asarray(img.convert('L')) < 240
    bar, th, fpx, m = nm * ppn, max(6, round(H / 160)), round(H / 22), round(H / 25)
    bw, bh = int(bar), int(th + fpx * 1.6)
    corners = {'bottom right': (W - m - bw, H - m - bh), 'bottom left': (m, H - m - bh),
               'top right': (W - m - bw, m), 'top left': (m, m)}
    name, (x0, y0) = min(corners.items(), key=lambda kv: fg[kv[1][1]:kv[1][1] + bh, kv[1][0]:kv[1][0] + bw].mean())
    d = ImageDraw.Draw(img); f = font(fpx); lab = f'{nm:g} nm'
    ty = y0; by = y0 + fpx * 1.3                                   # label above the bar
    d.text((x0 + bw / 2 - d.textlength(lab, font=f) / 2, ty), lab, font=f, fill=INK)
    d.rectangle([x0, by, x0 + bw, by + th], fill=INK)
    stem = os.path.splitext(png)[0]
    img.save(stem + '_bar.png')
    open(stem + '_bar.svg', 'w').write(
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
        f'<image xlink:href="{os.path.basename(png)}" x="0" y="0" width="{W}" height="{H}"/>\n<g id="scale_bar">\n'
        f'<rect x="{x0}" y="{by:.0f}" width="{bw}" height="{th}" fill="#3d3a33"/>\n'
        f'<text x="{x0 + bw / 2:.0f}" y="{ty + fpx * 0.9:.0f}" text-anchor="middle" font-family="Arial" font-size="{fpx}" fill="#3d3a33">{lab}</text>\n'
        '</g>\n</svg>\n')
    print(f'{png}: {lab} = {bar:.0f} px, placed {name}; wrote {stem}_bar.png and {stem}_bar.svg')


if __name__ == '__main__':
    main()
