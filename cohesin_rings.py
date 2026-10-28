"""cohesin_rings.py -- PyMOL: draw one torus per cohesin tether.  run cohesin_rings.py"""
from pymol.cgo import BEGIN, END, TRIANGLES, NORMAL, VERTEX, COLOR
from pymol import cmd
import math

def torus(name, c, ax, R, r, col=(1.0, 0.75, 0.1), nu=72, nv=16):
    ax = [x / math.sqrt(sum(y*y for y in ax)) for x in ax]
    t = [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0]
    e1 = [t[1]*ax[2]-t[2]*ax[1], t[2]*ax[0]-t[0]*ax[2], t[0]*ax[1]-t[1]*ax[0]]
    m = math.sqrt(sum(x*x for x in e1)); e1 = [x/m for x in e1]
    e2 = [ax[1]*e1[2]-ax[2]*e1[1], ax[2]*e1[0]-ax[0]*e1[2], ax[0]*e1[1]-ax[1]*e1[0]]
    def pt(i, j):
        u, v = 2*math.pi*i/nu, 2*math.pi*j/nv
        rad = [math.cos(u)*e1[k] + math.sin(u)*e2[k] for k in range(3)]
        nrm = [math.cos(v)*rad[k] + math.sin(v)*ax[k] for k in range(3)]
        return [c[k] + R*rad[k] + r*nrm[k] for k in range(3)], nrm
    obj = [BEGIN, TRIANGLES, COLOR, *col]
    for i in range(nu):
        for j in range(nv):
            for (a, b) in ((i, j), (i+1, j), (i+1, j+1), (i, j), (i+1, j+1), (i, j+1)):
                p, n = pt(a, b); obj += [NORMAL, *n, VERTEX, *p]
    obj.append(END)
    cmd.load_cgo(obj, name)

torus('ring_T_I5', [2605.08, -1807.67, -114.34], [0.2242, 0.9741, -0.03], 200.0, 12.0)
torus('ring_T_I7', [-1516.82, -2032.11, -1879.08], [-0.5962, -0.7981, 0.0874], 200.0, 12.0)
torus('ring_T_I2', [2002.56, -3853.72, -532.86], [0.8243, 0.2536, 0.5062], 200.0, 12.0)
cmd.group('cohesin_rings', 'ring_*')
