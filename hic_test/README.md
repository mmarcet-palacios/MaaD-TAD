# hic_test/ — builds from Micro-C alone: IGHC and PWRN (Fig 10)

Run from this folder (`cd hic_test`) in the `chromatin` environment. The PRAMEF pipeline scripts in the root folder are used unchanged; `build_locus.py` sets the locus start for its own process only.

| Step | Script | Output |
|---|---|---|
| Micro-C (H1, HFFc6; 2 kb) and SEQSIM maps of both windows | `extract_hic_test.py` (needs the two .mcool files in ../control_loci, the GRCh38 FASTA in ../data and ../control_loci/seqsim_windows_fast.exe) | `<name>_microc.npz/.csv`, `<name>_SEQSIM.csv`, `<name>.fa` |
| CSynth contact files | `to_contacts.py` | `<name>_<map>.contacts` (folded in CSynth -> `<name>_<map>.xyz`) |
| Full-atom build (Folder, Placer, Builder, Checker; no homology, no tethers) | `build_locus.py PWRN_chr15 IGHC_chr14` | `<name>_H1/` (frames.npz, blocks_final/, qc) |
| GRCh38 bases on the blocks | `thread_locus.py` (called by `build_locus.py`) | `<name>_H1/blocks_seq/` |
| Agreement with Micro-C | `score_models.py` | score_models_results.txt |
| Homology enrichment on the test loci | `homology_test_loci.py` | homology_test_loci_results.txt |
| Fig 10A, 10B | `fig_seqsim_maps.py`, `fig_models_maps.py` | figS_seqsim_A_*.png/svg, figS_models_A_*.png/svg |
| Fig 10C, 10D | `make_figS_models_pml.py` writes `figS_IGHC.pml`, `figS_PWRN.pml` (run in PyMOL from this folder); `label_overlay.py`, `add_scalebar.py` add labels and scale bars | figS_*_B_labelled.png, figS_*_C_bar.png |

Windows: IGHC chr14:105,540,001–105,750,000 (105 bins); PWRN chr15:24,130,001–24,440,000 (155 bins); GRCh38.
The maps, CSynth paths and the two builds are in the Zenodo data record (test_loci_IGHC_PWRN.zip). Unpack it here.
