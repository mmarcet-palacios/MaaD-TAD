"""make_fig2a.py -- Fig 2A: SEQSIM map of the locus as a PNG (3 px per 2 kb section).
Input: data/seqsim_full.npy, the 650 x 650 SEQSIM matrix (sections 1..650 of 2 kb from chr1:12,447,321, GRCh38;
       sections overlapping the assembly gap are NaN). Output: fig2a_seqsim.png (975 x 975)."""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from PIL import Image
S = np.load('data/seqsim_full.npy')
C = S[162:487, 162:487]                     # sections 163-487 = chr1:12,771,321-13,421,320 (325 x 325)
stops = ['#fffaf2', '#fbebd5', '#f8dbb7', '#f2c493', '#e9a466', '#e0833a', '#ca6a2d', '#b34f20', '#953a16', '#6f280e', '#4a1606']
cm = LinearSegmentedColormap.from_list('seqsim', stops)
norm = LogNorm(vmin=0.4, vmax=100, clip=True)  # log scale, 0.4 (floor) to 100 (self)
rgb = cm(norm(np.nan_to_num(C, nan=0.4)))[:, :, :3]
rgb[np.isnan(C)] = np.array([227, 223, 210]) / 255   # assembly gap (no sequence) in grey
Image.fromarray((rgb * 255).round().astype(np.uint8)).resize((975, 975), Image.NEAREST).save('fig2a_seqsim.png')
