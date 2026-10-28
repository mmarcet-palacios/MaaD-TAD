# fig7c_render.pml -- Fig 7C from the final model with the open-source PyMOL of the chromatin environment
# (no licence watermark). Run headless:  pymol -cq fig7c_render.pml   (from MaaD-TAD\pymol)
# Only the tether atoms are shown, so ray-tracing is quick and light on memory.
reinitialize
set max_threads, 24
set hash_max, 255
set auto_dss, 0
python
from pymol import cmd
import os
for i in range(68):
    cmd.load(os.path.join(r"../blocks_final", "block_%02d.cif" % i), "block_%02d" % i)
print("model loaded", flush=True)
python end
@fig7c.pml
ring_size 45
ring_info
view fig7C_side
png ../figures/final/fig7c_side_v2.png, width=6559, height=3531, ray=1
view fig7C_through
png ../figures/final/fig7c_through_v2.png, width=6559, height=3531, ray=1
print("FIG 7C DONE")
