#!/bin/sh
# run_ensemble.sh FIRST LAST WORKERS [K OUTDIR] -- build members FIRST..LAST (Folder + stage 1), skipping finished ones
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
K=${4:-0.007}; OUT=${5:-ensemble}; mkdir -p $OUT
seq "$1" "$2" | xargs -P "$3" -I{} sh -c "m=\$(printf '$OUT/m%03d' {}); [ -f \$m/member.json ] || python ensemble_member.py --seed {} --k $K --out $OUT --no-relax > \$m.out 2>&1"
