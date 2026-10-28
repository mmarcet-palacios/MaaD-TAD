# control_loci/ — homology maps, contact enrichment and calibration (Figs 1B, 8 and 9; S3 Table)

Run from this folder (`cd control_loci`). Scripts that read Micro-C expect the 4D Nucleome files here:
4DNFI9GMP2J8.mcool (H1) and 4DNFI643OYP9.mcool (HFFc6).

Windows (GRCh38, 2 kb bins): PRAMEF chr1:12,771,321 (325 bins) · PCDH chr5:140,700,000 (450) · ZNF_chr19 chr19:57,000,000 (500).

| Analysis | Script | Output |
|---|---|---|
| Fig 1B extract (2 kb, chr1:11.8–14.4 Mb) | `extract_fig1b.py` | fig1b_microc_H1.npz |
| Homology maps and identity of homologous copies (S3 Table) | `control_identity.py` (parasail; fetches sequence from UCSC, cached in identity_cache/) | identity_pairs.csv, identity_summary.txt |
| CSynth-like fold of a homology map | `csynth_like.py` | `<locus>_csynth_like.xyz` |
| Contact enrichment test (Fig 8B) and fold from Micro-C | `compare_hic.py --map 4DNFI9GMP2J8.mcool --locus PCDH --res 2000 --fold` | `<locus>_hic_compare.json`, `<locus>_hic_fold.xyz` |
| Calibration of the homology weight (Fig 8C) | `sweep_calibrate.py`, then `pramef_calibration.py` | sweep_calibrate.json, pramef_cal.json |
| Mapping test of the 145 kb line (Fig 9) | `misassign_pramef.py ../locus_GRCh38.fa fig1b_microc_H1.npz` (mappy) | misassign_report.txt, misassign_table.csv |
| The line in HFFc6 (Fig 9D) | `check_stripe.py` | printed table |
| Fig 9B, C maps (5 kb) | `extract_stripe_maps.py` | stripe_maps.npz |
| Fast SEQSIM for 2 kb windows | `build_seqsim_fast.bat` builds `seqsim_windows_fast.exe` from `seqsim_windows_fast.cpp` (used by `../hic_test/extract_hic_test.py`) | |

Homology map scale: alignment scores (repeats removed, both strands) are rescaled to the SEQSIM scale on PRAMEF by quantile matching; homologous pairs score > 1.1 and controls ≤ 0.6 (Materials and methods). The input maps (`<locus>_seqsim_like.csv`) and all outputs above are in the Zenodo data record (figure_source_data.zip).
