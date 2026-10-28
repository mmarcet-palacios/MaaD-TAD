# s1_video.pml -- S1 Video of the MaaD-TAD paper: one 360 degree turn of the full-atom model in Fig 5A colours
#   DNA blue, histones orange, white background, depth cue on (as Fig 5A); no gene colouring, no labels
#   check_plain.png        : first frame, written before the turn so the colours and framing can be checked
#   s1_video_frames/       : 900 frames (0.4 degrees per frame, 30 s at 30 fps), 1920 x 1080
# Camera: a rotation of the Fig 5A view, centred on the middle of the model,
# far enough back that the whole model stays in frame at every angle. Draw mode only (no ray tracing).
#
# Run in the PyMOL window, one line at a time:
#     cd MaaD-TAD/pymol
#     @load_final_model.pml          (wait until all 68 blocks have loaded)
#     @s1_video.pml
# For a quick test, change N = 900 to N = 12 in the python block below.
# Delete the old frames in s1_video_frames first: a 900-frame run overwrites f0000-f0449 but not the folder itself.

python
from pymol import cmd
import os, math

OUT = r"../figures/movie"
FOLDER = "s1_video_frames"
W, H = 1920, 1080
N = 900                             # frames per turn (0.4 degrees per frame)
FOV = 20.0                          # degrees, orthoscopic (as Fig 5A)
R_FIT = 6200.0                      # Angstrom: radius kept inside the frame at every angle (model radius about its middle: 5626 A)
COL_DNA = [0.247, 0.435, 0.659]     # Fig 5A colours (as in the Fig 5 legend)
COL_HIST = [0.851, 0.451, 0.122]
ROT = [-0.122136, 0.832715, 0.540059, -0.508779, 0.414676, -0.754440, -0.852188, -0.366915, 0.373025]
BLOCKS = ["block_%02d" % i for i in range(68)]

# ---- styling: every block set explicitly
cmd.set_color("m_dna", COL_DNA)
cmd.set_color("m_hist", COL_HIST)
cmd.set("suspend_updates", 1)
for obj in ("labels_movie", "leaders_movie"):          # left over from the gene video, if loaded
    if obj in cmd.get_names("all"):
        cmd.delete(obj)
cmd.show_as("lines", "all")
cmd.color("grey50", "all")
for obj in BLOCKS:
    cmd.color("m_dna", "%s and chain A+B" % obj)
    cmd.color("m_hist", "%s and chain P*" % obj)
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("orthoscopic", 1)
cmd.set("field_of_view", FOV)
cmd.set("depth_cue", 1)
cmd.deselect()
cmd.set("suspend_updates", 0)

# ---- check: each pair of numbers must match
print("DNA atoms %d, coloured %d" % (cmd.count_atoms("chain A+B"), cmd.count_atoms("chain A+B and color m_dna")), flush=True)
print("histone atoms %d, coloured %d" % (cmd.count_atoms("chain P*"), cmd.count_atoms("chain P* and color m_hist")), flush=True)

# ---- camera: Fig 5A orientation, centred on the middle of the model
ext = cmd.get_extent("all")
C = [(ext[0][k] + ext[1][k]) / 2.0 for k in range(3)]
D = R_FIT / math.tan(math.radians(FOV / 2.0))
SLAB = R_FIT + 300.0
VIEW = ROT + [0.0, 0.0, -D] + C + [D - SLAB, D + SLAB, -FOV]

def save_frame(path):
    cmd.do("draw %d, %d" % (W, H))
    cmd.png(path)

os.makedirs(OUT, exist_ok=True)
cmd.set_view(VIEW)
save_frame(os.path.join(OUT, "check_plain.png"))
print("wrote check_plain.png", flush=True)

d = os.path.join(OUT, FOLDER)
os.makedirs(d, exist_ok=True)
cmd.set_view(VIEW)
for i in range(N):
    save_frame(os.path.join(d, "f%04d.png" % i))
    cmd.turn("y", 360.0 / N)
print(FOLDER, "done:", N, "frames in", d, flush=True)
python end
