#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
RB=${RB:-$ROOT/build/src/rabbitbin}
OUTROOT=${OUTROOT:-$ROOT/experiments/architecture_local/results}
THREADS=${THREADS:-64}
REPEATS=${REPEATS:-1}
read -r -a modes <<< "${MODES:-baseline coverage_first}"
if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(marine plant_associated strain_madness); fi
for dataset in "${datasets[@]}"; do
  mkdir -p "$OUTROOT/$dataset"
  for ((rep=1; rep<=REPEATS; ++rep)); do
    for mode in "${modes[@]}"; do
      prefix="$OUTROOT/$dataset/$mode.r$rep"
      extra=()
      if [[ "$mode" == coverage_first ]]; then extra=(RABBIT_CANDIDATE_COVERAGE=1); fi
      /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
        env -u RABBIT_CANDIDATE_COVERAGE -u RABBIT_LPA_SCORE -u RB_LPA_AUDIT \
        RB_TIMING=1 "${extra[@]}" "$RB" bin \
        --assembly "/home/bigssd/zt/cami2/$dataset/CAMI2_${dataset}_GoldStandardAssembly.fasta" \
        --bam-list "/home/bigssd/zt/CAMI2_remap_bams/$dataset/bam.list" \
        --output "$prefix" --threads "$THREADS" --seed 42 \
        --sketch-m 500 --max-edges 200 --dual-depth 5 \
        --min-bin-size 200000 --no-bin-fasta \
        --export-retained-graph "$prefix.rbedge" --save-cache "$prefix.cache" \
        > "$prefix.log" 2>&1
      "$RB" amber --gold "/home/bigssd/zt/runs/cami2_benchmark/$dataset/prep/gold_len.binning" \
        --members "$prefix.members.tsv" --min-bin-size 200000 --threads "$THREADS" \
        --quiet > "$prefix.amber" 2>&1
      if [[ "$mode" == baseline ]]; then
        cmp "$prefix.members.tsv" "$ROOT/experiments/candidate_stage/cluster_results/$dataset/baseline.members.tsv"
      fi
      printf '%s %s repeat %s done\n' "$dataset" "$mode" "$rep"
    done
  done
done
