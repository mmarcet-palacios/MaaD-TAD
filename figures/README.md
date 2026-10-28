# figures/ — scripts for Figs 1 and 2 and S1 Fig (S1 File)

Run from this folder (`cd figures`) in the `chromatin` environment. Inputs are in `data/`. Each script writes an SVG fragment or PNG that was placed into the figure layout; figure numbers, titles and legends are in the paper. S1 Text gives the colour scales and other drawing details.

| Panel | Script | Inputs | Output |
|---|---|---|---|
| Fig 1A | `fig1a.py` | data/genes.tsv, data/exons.tsv, data/nucleosomes_GRCh38_resolved.bed | SVG fragment |
| Fig 1B | `../control_loci/extract_fig1b.py`, then `make_fig1b.py` | 4DNFI9GMP2J8.mcool (4DN portal, not included) | fig1b_microc_H1.npz (copy it here), fig1b_microc_map.png, f1b_svg.frag |
| Fig 1C | `fig1c.py` | data/homology_restraints.json, data/genes.tsv, data/new_sites.json | SVG fragment |
| Fig 2A | `make_fig2a.py` | data/seqsim_full.npy (SEQSIM matrix of the locus, 2 kb sections) | fig2a_seqsim.png |
| S1 Fig | `s1fig/stage1_s1fig.py`, then `s1fig/fig_s1.py` | the Placer inputs (data/) | s1fig/anch_s1fig.npz, s1fig/free_s1fig.npz, fig_s1_panels.svgfrag |

`fig1b_microc_H1.npz` is also in the Zenodo data record (figure_source_data.zip), so Fig 1B can be redrawn without the 13 GB .mcool file.

S1 Fig: `stage1_s1fig.py` is `stage1_frames.py` with extra outputs (steering rotation and along-path offset of every nucleosome). Run it from the repository root with the Placer settings, once anchored and once without anchoring:

```
python figures/s1fig/stage1_s1fig.py --bed data/nucleosomes_tethered.bed --guide data/guide_tethered.xyz --gaps data/gap_sections.txt --domain 1011.5 --coil-radius 0 --smooth 1 --calib 0 --look 4 --out figures/s1fig/anch.npz
python figures/s1fig/stage1_s1fig.py --bed data/nucleosomes_tethered.bed --guide data/guide_tethered.xyz --gaps data/gap_sections.txt --domain 1011.5 --coil-radius 0 --smooth 1 --calib 0 --look 4 --no-anchor --out figures/s1fig/free.npz
python figures/s1fig/fig_s1.py
```

Coordinates: GRCh38 chr1. Locus chr1:12,771,346–13,423,345. Assembly gap 12,954,385–13,004,384.
