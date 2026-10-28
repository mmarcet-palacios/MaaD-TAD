"""
Layer 3, step 2: ring geometry for each cohesin tether.

    python tethers_layer3.py --frames frames.npz      -> data/tethers.json, cohesin_rings.py

For each tether, the naked linker crossing the cohesin site and the one crossing
its mirror are found. A ring centre is placed between them; each duplex gets a
waypoint 17.5 A either side of the centre (duplex axes 35 A apart, i.e. ~15 A
between the two 20 A-wide helices) and a direction along the ring axis. Stage 2
builds each of these linkers as two halves that meet exactly at its waypoint.

The ring is drawn as a torus (PyMOL CGO) around both duplexes: --ring-radius
is the centreline radius, chosen by default as the largest of 60-200 A that
keeps the torus clear of every histone core; replace later with a real
cohesin structure superposed on the same centre/axis.
"""
import argparse, json
import numpy as np

import closure as cl
from build_zigzag import canonical_1kx5, place

HALF_GAP = 17.5     # duplex axis to ring centre (35 A between duplex axes: phosphates >= ~6 A apart)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--restraints", default="data/homology_restraints.json")
    ap.add_argument("--out", default="data/tethers.json")
    ap.add_argument("--ring-radius", type=float, default=None)
    ap.add_argument("--tube", type=float, default=12.0)
    a = ap.parse_args()

    F = np.load(a.frames); O, Rm = F["origins"], F["rotations"]
    dy = F["dyad_pos"].astype(int); nrl = F["nrl"].astype(int)
    le, lx = canonical_1kx5()
    R = json.load(open(a.restraints))
    out = []
    for t in R["tethers"]:
        sites = [("cohesin", 0.5 * (t["cohesin"][0] + t["cohesin"][1])), ("mirror", float(t["mirror"]))]
        L = []
        for kind, g in sites:
            i = int(np.searchsorted(dy, g) - 1)                 # linker i joins core i -> i+1
            ex = place((O[i], Rm[i]), lx); en = place((O[i + 1], Rm[i + 1]), le)
            nb = nrl[i] - 147
            first_bp = dy[i] + 74                               # first linker bp (genomic)
            k = int(round(g - first_bp))                        # bp index of the site within the linker
            k = int(np.clip(k, 8, nb - 9))
            q = ex[0] + (k / nb) * (en[0] - ex[0])
            L.append(dict(kind=kind, linker=i, nb=int(nb), split_bp=k, genomic=int(g),
                          exit=ex[0], entry=en[0], q=q, d=(en[0] - ex[0]) / np.linalg.norm(en[0] - ex[0])))
        d1, d2 = L[0]["d"], L[1]["d"]
        s = 1.0 if np.dot(d1, d2) >= 0 else -1.0
        n_ax = d1 + s * d2; n_ax /= np.linalg.norm(n_ax)
        C = 0.5 * (L[0]["q"] + L[1]["q"])
        u = L[1]["q"] - L[0]["q"]; u -= np.dot(u, n_ax) * n_ax
        u = u / np.linalg.norm(u) if np.linalg.norm(u) > 1e-6 else np.cross(n_ax, [1, 0, 0])
        for sign, lk in ((-1, L[0]), (1, L[1])):
            lk["waypoint"] = C + sign * HALF_GAP * u
            lk["dir"] = n_ax * (1.0 if np.dot(lk["d"], n_ax) >= 0 else -1.0)
            # reach: each half must span its straight distance with slack
            c1 = lk["split_bp"] * cl.B_RISE; c2 = (lk["nb"] - lk["split_bp"]) * cl.B_RISE
            lk["ext_half1"] = float(np.linalg.norm(lk["waypoint"] - lk["exit"]) / c1)
            lk["ext_half2"] = float(np.linalg.norm(lk["entry"] - lk["waypoint"]) / c2)
        # ring radius: largest keeping the torus clear of histone cores (~55 A disc radius)
        radii = [a.ring_radius] if a.ring_radius else list(range(200, 55, -10))
        best = None
        for r in radii:
            # distance from each core centre to the torus centreline circle
            v = O - C; h = v @ n_ax; rad = np.linalg.norm(v - np.outer(h, n_ax), axis=1)
            dcirc = np.sqrt(h ** 2 + (rad - r) ** 2)
            clr = float(dcirc.min() - a.tube - 55.0)
            if best is None or clr > best[1]:
                best = (r, clr)
            if a.ring_radius is None and clr > 0:
                best = (r, clr); break
        ring = dict(centre=C.tolist(), axis=n_ax.tolist(), radius=float(best[0]), tube=a.tube,
                    clearance_to_cores_A=round(best[1], 1))
        out.append(dict(id=t["id"], block=t["block"], ring=ring, duplex_axis_separation_A=2 * HALF_GAP,
                        linkers=[dict(kind=lk["kind"], linker=lk["linker"], nb=lk["nb"], split_bp=lk["split_bp"],
                                      genomic=lk["genomic"], waypoint=lk["waypoint"].tolist(),
                                      dir=lk["dir"].tolist(), ext_half1=round(lk["ext_half1"], 2),
                                      ext_half2=round(lk["ext_half2"], 2)) for lk in L]))
        print(f"{t['id']}: ring centre set; radius {best[0]:.0f} A (torus clearance to nearest core {best[1]:+.0f} A)")
        for lk in L:
            print(f"   {lk['kind']:<8} linker {lk['linker']:5d} ({lk['nb']} bp naked, site at bp {lk['split_bp']}): "
                  f"extension to waypoint {lk['ext_half1']:.2f}, from waypoint {lk['ext_half2']:.2f}"
                  + ("   <-- too tight" if max(lk["ext_half1"], lk["ext_half2"]) > 0.9 else ""))
    json.dump(out, open(a.out, "w"), indent=1)

    lines = ['"""cohesin_rings.py -- PyMOL: draw one torus per cohesin tether.  run cohesin_rings.py"""',
             "from pymol.cgo import BEGIN, END, TRIANGLES, NORMAL, VERTEX, COLOR",
             "from pymol import cmd", "import math", "",
             "def torus(name, c, ax, R, r, col=(1.0, 0.75, 0.1), nu=72, nv=16):",
             "    ax = [x / math.sqrt(sum(y*y for y in ax)) for x in ax]",
             "    t = [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0]",
             "    e1 = [t[1]*ax[2]-t[2]*ax[1], t[2]*ax[0]-t[0]*ax[2], t[0]*ax[1]-t[1]*ax[0]]",
             "    m = math.sqrt(sum(x*x for x in e1)); e1 = [x/m for x in e1]",
             "    e2 = [ax[1]*e1[2]-ax[2]*e1[1], ax[2]*e1[0]-ax[0]*e1[2], ax[0]*e1[1]-ax[1]*e1[0]]",
             "    def pt(i, j):",
             "        u, v = 2*math.pi*i/nu, 2*math.pi*j/nv",
             "        rad = [math.cos(u)*e1[k] + math.sin(u)*e2[k] for k in range(3)]",
             "        nrm = [math.cos(v)*rad[k] + math.sin(v)*ax[k] for k in range(3)]",
             "        return [c[k] + R*rad[k] + r*nrm[k] for k in range(3)], nrm",
             "    obj = [BEGIN, TRIANGLES, COLOR, *col]",
             "    for i in range(nu):",
             "        for j in range(nv):",
             "            for (a, b) in ((i, j), (i+1, j), (i+1, j+1), (i, j), (i+1, j+1), (i, j+1)):",
             "                p, n = pt(a, b); obj += [NORMAL, *n, VERTEX, *p]",
             "    obj.append(END)",
             "    cmd.load_cgo(obj, name)", ""]
    for t in out:
        r = t["ring"]
        lines.append(f"torus('ring_{t['id']}', {[round(x, 2) for x in r['centre']]}, "
                     f"{[round(x, 4) for x in r['axis']]}, {r['radius']}, {r['tube']})")
    lines.append("cmd.group('cohesin_rings', 'ring_*')")
    open("cohesin_rings.py", "w").write("\n".join(lines) + "\n")
    print(f"wrote {a.out} and cohesin_rings.py")


if __name__ == "__main__":
    main()
