# MaaD-TAD

**Molecular all-atom Display of TADs**: a Python pipeline that turns a 3D path of a chromatin domain (from CSynth, folded on a homology map or on Micro-C) into an all-atom model of DNA and histones, one nucleosome at a time.

This repository holds the code, inputs and figure scripts of:

> Marcet-Palacios M, MacIsaac M, Santos JRL, Moon EB, Befus AD. *MaaD-TAD: an all-atom model of the duplicated PRAMEF/HNRNPCL chromatin domain built from sequence homology.* Submitted to PLOS Computational Biology.

It is S1 File of the paper. The full-atom model and the larger data files are in the Zenodo data record (see [Data](#data)).

## The four modules

| Module | Driver | What it does | Scripts it runs |
|---|---|---|---|
| Folder | `folder.py` | Refolds the CSynth path so that it can hold a chromatin fibre; applies the homology pairing restraints and cohesin tether terms | `refold_guide.py`, `make_naked.py`, `refold_tether.py` |
| Placer | `placer.py` | Places the nucleosomes along the refolded path (two-start zigzag, anchored to genomic position) and removes collisions | `stage1_frames.py`, `relax_frames.py`, `tethers_layer3.py`, `pairing_scorecard.py` |
| Builder | `builder.py` | Builds DNA and histones in atoms in blocks of 50 nucleosomes, closes every linker, threads the GRCh38 sequence and repairs the histone tails | `run_stage2.py` (+ `stage2_block.py`), `thread_sequence.py`, `repair_tails.py`, `compare_base_contacts.py` |
| Checker | `checker.py` | Checks the whole model (continuity across seams, backbone geometry, long-range clashes) and rebuilds only the blocks that fail | `stage4_qc.py`, `relax_frames.py --problems`, `check_standin.py`, `run_stage2.py --adopt`, `repair_tails.py --from-qc` |

The other modules in the root folder are libraries used by these scripts (`bpstep.py`, `closure.py`, `backbone_solver.py`, `fibre_*.py`, `build_zigzag.py`, `analyse_b2.py`, `build_array9b.py`, `rebuild_test.py`, `dna_proxy.py`, `mini_qc.py`, `cohesin_rings.py`, `resolve_alternatives.py`). They import each other by name, so keep them together in this folder. Every parameter value is listed in S2 Table of the paper.

## Install

Tested on Windows 11 (Intel Core i9-14900K, 24 cores, 128 GB RAM) with Python 3.11. The scripts also run on Linux and macOS.

```
conda env create -f environment.yml
conda activate chromatin
```

`pymol-open-source` is needed by `analyse_b2.py` (it reads the 1KX5 zigzag template in `data/b2_onelayer.pse`) and for the figure renders in `pymol/`. The two C++ programs (`SEQSIM/seqsim.cpp`, `control_loci/seqsim_windows_fast.cpp`) need a C++ compiler with OpenMP (Visual Studio 2022 or g++).

## Run

From this folder, in the `chromatin` environment:

```
RUN_PRAMEF.bat            (Windows)
./run_pramef.sh           (Linux, macOS)
```

which runs, in order:

```
python folder.py
python placer.py
python builder.py --steps build
python checker.py --steps qc repair
python builder.py --steps seq tails
python checker.py --steps final
```

Each driver takes `--steps` (run part of a module), `--dry-run` (print the commands only) and `--workers`. `builder.py --size N` changes the block size (default 50 nucleosomes). The full build takes a few hours on 24 cores and about 15 GB of disk; the model ends up in `blocks_final/` (`block_00.cif` ... `block_67.cif`, mmCIF, about 55 million atoms). Folder's outputs (`data/guide_refolded.xyz`, `data/guide_tethered.xyz`, `data/nucleosomes_tethered.bed`) are included, so the build can start at Placer.

Residue numbering: chain A residue r = genomic position − 12,771,345; chain B residue r = 13,423,346 − genomic position (GRCh38 chr1). Histone chains are named P*.

## Folders

| Folder | Contents |
|---|---|
| (root) | the pipeline: drivers, scripts and libraries; `locus_GRCh38.fa` (chr1:12,771,346–13,423,345) |
| `data/` | inputs: 1KX5 nucleosome (`1KX5.cif`), zigzag template (`b2_onelayer.pse`), CSynth path (`gap-filled.xyz`), assembly gap (`gap_sections.txt`), nucleosome positions (`nucleosomes_GRCh38_resolved.bed`), homology blocks (`homology_restraints.json`), tether sites (`ctcf_cohesin_mirror_sites.csv`), Folder outputs, and `verify_nucleosome_source.py` |
| `control_loci/` | homology maps and CSynth-like folds of the control loci, the contact enrichment test and the calibration of the homology weight (Fig 8, S3 Table), the 145 kb line tests (Fig 9) |
| `hic_test/` | the builds from Micro-C alone, IGHC and PWRN (Fig 10) |
| `figures/` | scripts for Figs 1 and 2 and S1 Fig, with their inputs in `figures/data/` |
| `pymol/` | PyMOL scripts for the full-atom renders (Figs 4A, 5, 6C, 7C, 9E; S1 Video) |
| `SEQSIM/` | source of the SEQSIM sequence-similarity program used for the homology map (Fig 2A) |

## Figures and analyses

| Figure / table | Script(s) | Folder |
|---|---|---|
| Fig 1A | `fig1a.py` | figures |
| Fig 1B | `extract_fig1b.py` (needs 4DNFI9GMP2J8.mcool), `make_fig1b.py` | control_loci, figures |
| Fig 1C | `fig1c.py` | figures |
| Fig 2A | `make_fig2a.py` (SEQSIM map `figures/data/seqsim_full.npy`) | figures |
| Table 1 | `check_table1b.py` (re-exports the ChIP evidence from the UCSC API) | root |
| Fig 3, Fig 4 | `refold_guide.py`, `refold_tether.py`, `stage1_frames.py`; stage 5: `fig4_stage5_view.pml` | root, pymol |
| Fig 5 | views in `fig5_views.pml`; colours and views together: `fig5_style_views.pml` | pymol |
| Fig 6C | `fig6c_setup.pml` (loads blocks 29–40, runs `fig6c.pml`) | pymol |
| Fig 7C | `fig7c_setup.pml` or `run_fig7c.bat` (runs `fig7c.pml`, which fits the cohesin ring) | pymol |
| Fig 8 | `csynth_like.py`, `compare_hic.py --fold`, `sweep_calibrate.py`, `pramef_calibration.py`; identity of homologous copies (S3 Table): `control_identity.py` | control_loci |
| Fig 9 | `misassign_pramef.py` (mapping test), `check_stripe.py` (HFFc6), `extract_stripe_maps.py` (panels B, C); panel E: `fig9e_setup.pml` | control_loci, pymol |
| Fig 10 | `extract_hic_test.py`, `to_contacts.py`, `build_locus.py`, `thread_locus.py`, `score_models.py`, `homology_test_loci.py`, `fig_seqsim_maps.py` (A), `fig_models_maps.py` (B), `make_figS_models_pml.py` (C, D) | hic_test |
| S1 Fig | `stage1_s1fig.py` (stage 1 with extra outputs), `fig_s1.py` | figures/s1fig |
| S2 Fig | `scale_feasibility.py` | root |
| S3 Fig | `backfold_axis.py` | root |
| S1 Video | `load_final_model.pml`, then `s1_video.pml` | pymol |
| nucleosome source check | `verify_nucleosome_source.py` (needs pyliftover) | data |

PyMOL scripts load the blocks from `../blocks_final`, so start PyMOL in `pymol/` (`cd pymol`) after unpacking the model there (see below). Renders are written to `figures/final/`.

## Data

Large files are not in this repository. They are in the Zenodo data record [doi:10.5281/zenodo.23718429](https://doi.org/10.5281/zenodo.23718429) (CC BY 4.0):

- `PRAMEF_full_atom_model.zip`: the 68 blocks of the final model (`blocks_final/`, mmCIF), seam files and the tail-repair table. Unpack into this folder so that `blocks_final/` sits next to `folder.py`.
- `PRAMEF_paths_frames_maps.zip`: CSynth, refolded and tethered paths, nucleosome frames, homology and tether restraints, the SEQSIM map and the H1 Micro-C extract, QC reports.
- `test_loci_IGHC_PWRN.zip`: the IGHC and PWRN builds of Fig 10 (frames, blocks, Micro-C and SEQSIM maps, CSynth paths).
- `figure_source_data.zip`: the numbers behind the figures and tables.

Micro-C maps are from the 4D Nucleome Data Portal: H1, 4DNFI9GMP2J8; HFFc6, 4DNFI643OYP9 (`.mcool`, 13–17 GB each). Scripts that read them expect the files in `control_loci/`. The GRCh38 genome (NCBI GCF_000001405.40, `GCF_000001405.40_GRCh38.p14_genomic.fna.gz`) goes in `data/` for `hic_test/extract_hic_test.py`.

## Cite

If you use MaaD-TAD, please cite the paper above and this software (`CITATION.cff`; Zenodo software DOI [10.5281/zenodo.23718187](https://doi.org/10.5281/zenodo.23718187)).

## Licence

Code: MIT (`LICENSE`). Data files in this repository and in the Zenodo record: CC BY 4.0. `data/1KX5.cif` is from the Protein Data Bank (PDB 1KX5; Davey et al., 2002, J Mol Biol 319:1097–1113).

Contact: Marcelo Marcet-Palacios, marcelo@ualberta.ca
