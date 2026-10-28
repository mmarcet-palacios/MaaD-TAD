# pymol/ — full-atom renders (Figs 4A, 5, 6C, 7C, 9E; S1 Video)

Unpack the model (Zenodo, PRAMEF_full_atom_model.zip) so that `blocks_final/` sits in the repository root, then start PyMOL in this folder:

```
cd MaaD-TAD/pymol
```

| Panel | Script | Notes |
|---|---|---|
| Fig 4A stage 5 | `@load_final_model.pml`, then `@fig4_stage5_view.pml` | all 68 blocks (5–20 min to load) |
| Fig 5A–C | `@load_final_model.pml`, then `@fig5_style_views.pml`; `view fig5A` / `fig5B` / `fig5C` | views also in `fig5_views.pml` |
| Fig 6C | `@fig6c_setup.pml` (loads blocks 29–40 and runs `fig6c.pml`); `view fig6C` | |
| Fig 7C | `@fig7c_setup.pml` (blocks 14, 15, 51, 52 and `fig7c.pml`, ring radius 4.5 nm); `view fig7C_side`. Headless, all blocks: `run_fig7c.bat` (runs `fig7c_render.pml`) | `fig7c.pml` fits the cohesin ring and prints the measurements of all three tethers |
| Fig 9E | `@fig9e_setup.pml` (loads blocks 42–64 and keeps chr1:13,184,371–13,395,000); `view d3d4` | |
| S1 Video | `@load_final_model.pml`, then `@s1_video.pml` | 900 frames, 1,920 × 1,080, in ../figures/movie/ |
| all panels in one run | `@render_all.pml` | draw mode (graphics card) |

Each setup script stores the view and leaves the final ray-trace to you; the size used for the paper is given in each file (6,559 × 3,531 pixels for Figs 5, 6C and 7C). `dna_constant.pml` and `fig4_stage5.pml` are the neutral colours of Fig 4.
