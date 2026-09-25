#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
RB=${RB:-$ROOT/build/src/rabbitbin}
GOLD_ROOT=${GOLD_ROOT:-/home/bigssd/zt/runs/cami2_benchmark}
THREADS=${THREADS:-64}
OUTROOT=${OUTROOT:-$ROOT/experiments/candidate_stage/cluster_results}

for dataset in marine plant_associated strain_madness; do
  cache="$OUTROOT/$dataset/baseline.cache"
  prefix="$OUTROOT/$dataset/no_split"
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.rebin.time" \
    "$RB" bin --load-cache "$cache" --no-split \
    --output "$prefix" --threads "$THREADS" --seed 42 \
    --sketch-m 500 --max-edges 200 --dual-depth 5 \
    --min-bin-size 200000 --no-bin-fasta > "$prefix.log" 2>&1
  "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
    --members "$prefix.members.tsv" --min-bin-size 200000 \
    --threads "$THREADS" --quiet > "$prefix.amber" 2>&1
done
