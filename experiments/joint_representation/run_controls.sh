#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
HERE="$ROOT/experiments/joint_representation"
PY=${PY:-$ROOT/experiments/candidate_stage/venv/bin/python}
RB=${RB:-$ROOT/build/src/rabbitbin}
THREADS=${THREADS:-64}
BASE_ROOT=${BASE_ROOT:-$ROOT/experiments/candidate_stage/cluster_results}
OUTROOT=${OUTROOT:-$HERE/results}
GOLD_ROOT=${GOLD_ROOT:-/home/bigssd/zt/runs/cami2_benchmark}
read -r -a methods <<< "${METHODS:-concat whiten cca}"
read -r -a specs <<< "${CONTROL_SPECS:-sum:joint:sum guarded:guarded:fisher guarded_sum:guarded:sum neutral:neutral:fisher weights_sum:retained:sum}"
if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(strain_madness marine plant_associated); fi
for dataset in "${datasets[@]}"; do
  out="$OUTROOT/$dataset"
  cache="$BASE_ROOT/$dataset/baseline.cache"
  "$PY" "$HERE/controls.py" --cache "$cache" \
    --retained "$BASE_ROOT/$dataset/retained.rbedge" --output "$out" --methods "${methods[@]}"
  for spec in identity:identity:fisher baseline_sum:identity:sum; do
    IFS=: read -r label graph score <<< "$spec"
    prefix="$out/$label"
    /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
      env RABBIT_LPA_SCORE="$score" "$RB" bin --load-cache "$cache" \
      --external-graph "$out/$graph.rbedge" --output "$prefix" \
      --threads "$THREADS" --seed 42 --min-bin-size 200000 --no-bin-fasta \
      --export-retained-graph "$prefix.scored.rbedge" > "$prefix.log" 2>&1
    "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
      --members "$prefix.members.tsv" --min-bin-size 200000 \
      --threads "$THREADS" --quiet > "$prefix.amber" 2>&1
  done
  cmp "$out/identity.members.tsv" "$BASE_ROOT/$dataset/baseline.members.tsv"
  for method in "${methods[@]}"; do
    for spec in "${specs[@]}"; do
      IFS=: read -r label graph score <<< "$spec"
      prefix="$out/$method.$label"
      /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
        env RABBIT_LPA_SCORE="$score" "$RB" bin --load-cache "$cache" \
        --external-graph "$out/$method.$graph.rbedge" --output "$prefix" \
        --threads "$THREADS" --seed 42 --min-bin-size 200000 --no-bin-fasta \
        --export-retained-graph "$prefix.scored.rbedge" > "$prefix.log" 2>&1
      "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
        --members "$prefix.members.tsv" --min-bin-size 200000 \
        --threads "$THREADS" --quiet > "$prefix.amber" 2>&1
    done
  done
done
