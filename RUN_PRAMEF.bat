@echo off
REM RUN_PRAMEF.bat -- the whole MaaD-TAD build of the PRAMEF/HNRNPCL locus (chr1:12,771,346-13,423,345, GRCh38),
REM in the order used for the paper. Run from the MaaD-TAD folder in the "chromatin" conda environment.
REM Takes a few hours on a 24-core workstation with 128 GB RAM; the model needs about 15 GB of disk.
REM Folder's outputs are supplied in data\, so the first line can be skipped.
cd /d %~dp0

python folder.py                                || exit /b 1
python placer.py                                || exit /b 1
python builder.py --steps build                 || exit /b 1
python checker.py --steps qc repair             || exit /b 1
python builder.py --steps seq tails             || exit /b 1
python checker.py --steps final                 || exit /b 1
echo MaaD-TAD build finished: the model is in blocks_final\ (block_00.cif ... block_67.cif)
