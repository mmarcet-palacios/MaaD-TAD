# fig7c_setup.pml -- Fig 7C (tether T_I5) from the FINAL model, loading only the 4 blocks around the ring
# Blocks 14, 15, 51 and 52 hold every nucleosome within 100 nm of the T_I5 ring, so the full model is not needed
# (about 3 M atoms instead of 52 M; loads in a minute or two).
# Fresh PyMOL:   @fig7c_setup.pml
# Then check the ring (ring_info prints which duplexes pass through it), adjust if needed
#   (ring_size R | ring_move du, dv, da | ring_tilt a, b | ring_fit | ring_fit min),
# and ray yourself:   view fig7C_side
#                     png ../figures/final/fig7c_side_ray.png, width=6559, height=3531, ray=1
# (T_I7 and T_I2 are not loaded, so the measurement table reports them as "not loaded"; that is expected.)
reinitialize
set auto_dss, 0
python
from pymol import cmd
import os
for i in (14, 15, 51, 52):
    cmd.load(os.path.join(r"../blocks_final", "block_%02d.cif" % i), "block_%02d" % i)
print("blocks 14, 15, 51, 52 of the final model loaded: %d atoms" % cmd.count_atoms("all"))
python end
@fig7c.pml
ring_size 45
view fig7C_side
print("Fig 7C ready: check ring_info above, then view fig7C_side and ray")
