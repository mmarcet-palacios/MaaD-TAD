# fig6c_setup.pml -- Fig 6C from the FINAL model: loads only blocks 29-40, colours genes, stores view fig6C
# Fresh PyMOL:   cd MaaD-TAD/pymol
#                @fig6c_setup.pml
# Then ray yourself:   view fig6C
#                      png ../figures/final/fig6c_ray.png, width=6559, height=3531, ray=1
reinitialize
set auto_dss, 0
python
from pymol import cmd
import os
for i in range(29, 41):                       # nucleosomes 1455-2012 live in blocks 29-40
    cmd.load(os.path.join(r"../blocks_final", "block_%02d.cif" % i), "block_%02d" % i)
print("blocks 29-40 of the final model loaded")
python end
@fig6c.pml
view fig6C
print("Fig 6C ready: view fig6C, then ray")
