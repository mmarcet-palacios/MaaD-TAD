#!/usr/bin/env bash
# run_pramef.sh -- the same build as RUN_PRAMEF.bat, for Linux or macOS. Run from the MaaD-TAD folder
# in the "chromatin" conda environment (conda env create -f environment.yml).
set -e
cd "$(dirname "$0")"
python folder.py
python placer.py
python builder.py --steps build
python checker.py --steps qc repair
python builder.py --steps seq tails
python checker.py --steps final
echo "MaaD-TAD build finished: the model is in blocks_final/ (block_00.cif ... block_67.cif)"
