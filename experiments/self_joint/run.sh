#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
PY=${PY:-$ROOT/experiments/candidate_stage/venv/bin/python}
RB=${RB:-$ROOT/build/src/rabbitbin}
OUTROOT=${OUTROOT:-$ROOT/experiments/self_joint/results}
BASE_ROOT=${BASE_ROOT:-$ROOT/experiments/candidate_stage/cluster_results}
FEATURE_ROOT=${FEATURE_ROOT:-$ROOT/experiments/joint_representation/results}
GOLD_ROOT=${GOLD_ROOT:-/home/bigssd/zt/runs/cami2_benchmark}
THREADS=${THREADS:-64}
RUN_LABEL=${RUN_LABEL:-learning}
TRAINING_LABEL=${TRAINING_LABEL:-}
export OPENBLAS_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS
read -r -a methods <<< "${METHODS:-joint block coverage}"
if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(marine plant_associated strain_madness); fi
for dataset in "${datasets[@]}"; do
  out="$OUTROOT/$dataset"
  base="$BASE_ROOT/$dataset"
  cache="$base/baseline.cache"
  if [[ "${SKIP_LEARN:-0}" != 1 ]]; then
    /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$out/$RUN_LABEL.time" \
      "$PY" "$ROOT/experiments/self_joint/learn.py" --cache "$cache" \
      --features "$FEATURE_ROOT/$dataset/features" \
      --halves "$out/halves" --retained "$base/retained.rbedge" --output "$out" \
      --threads "$THREADS" --methods "${methods[@]}" --training-label "$TRAINING_LABEL" > "$out/$RUN_LABEL.log" 2>&1
  fi
  for method in "${methods[@]}"; do
    for stage in candidates weights joint; do
      prefix="$out/sj_$method.$stage"
      graph="$out/sj_$method.joint.rbedge"
      extra=()
      if [[ "$stage" == weights ]]; then graph="$out/sj_$method.retained.rbedge"; fi
      if [[ "$stage" == candidates ]]; then extra=(--external-graph-coverage); fi
      /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
        env RABBIT_LPA_SCORE=fisher "$RB" bin --load-cache "$cache" \
        --external-graph "$graph" "${extra[@]}" --output "$prefix" --threads "$THREADS" \
        --seed 42 --min-bin-size 200000 --no-bin-fasta \
        --export-retained-graph "$prefix.scored.rbedge" > "$prefix.log" 2>&1
      "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
        --members "$prefix.members.tsv" --min-bin-size 200000 --threads "$THREADS" \
        --quiet > "$prefix.amber" 2>&1
    done
  done
  # Same kernel/topology/aggregation controls as the previous representation experiments.
  prefixed=()
  for method in "${methods[@]}"; do prefixed+=("sj_$method"); done
  env METHODS="${prefixed[*]}" OUTROOT="$OUTROOT" THREADS="$THREADS" RB="$RB" PY="$PY" \
    BASE_ROOT="$BASE_ROOT" GOLD_ROOT="$GOLD_ROOT" \
    bash "$ROOT/experiments/joint_representation/run_controls.sh" "$dataset"
  printf '%s self-supervised metric experiments done\n' "$dataset"
done
