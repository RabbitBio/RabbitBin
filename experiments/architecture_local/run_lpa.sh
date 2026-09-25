#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
RB=${RB:-$ROOT/build/src/rabbitbin}
OUTROOT=${OUTROOT:-$ROOT/experiments/architecture_local/results}
THREADS=${THREADS:-64}
CACHE_ROOT=${CACHE_ROOT:-$ROOT/experiments/candidate_stage/cluster_results}
CACHE_NAME=${CACHE_NAME:-baseline}
LABEL_PREFIX=${LABEL_PREFIX:-lpa}
read -r -a scores <<< "${SCORES:-fisher tail logsum}"
if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(marine plant_associated strain_madness); fi
for dataset in "${datasets[@]}"; do
  mkdir -p "$OUTROOT/$dataset"
  cache="$CACHE_ROOT/$dataset/$CACHE_NAME.cache"
  for score in "${scores[@]}"; do
    prefix="$OUTROOT/$dataset/${LABEL_PREFIX}_$score"
    /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
      env RB_LPA_PROF=1 RB_LPA_TRACE=1 RABBIT_LPA_SCORE="$score" \
      "$RB" bin --load-cache "$cache" --output "$prefix" \
      --threads "$THREADS" --seed 42 --min-bin-size 200000 --no-bin-fasta \
      > "$prefix.log" 2>&1
    "$RB" amber --gold "/home/bigssd/zt/runs/cami2_benchmark/$dataset/prep/gold_len.binning" \
      --members "$prefix.members.tsv" --min-bin-size 200000 --threads "$THREADS" \
      --quiet > "$prefix.amber" 2>&1
    if [[ "$score" == fisher ]]; then
      cmp "$prefix.members.tsv" "$CACHE_ROOT/$dataset/$CACHE_NAME.members.tsv"
    fi
    printf '%s %s done\n' "$dataset" "$score"
  done
done
