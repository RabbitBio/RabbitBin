#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
HERE="$ROOT/experiments/joint_representation"
PY=${PY:-$ROOT/experiments/candidate_stage/venv/bin/python}
RB=${RB:-$ROOT/build/src/rabbitbin}
THREADS=${THREADS:-64}
OUTROOT=${OUTROOT:-$HERE/results}
BASE_ROOT=${BASE_ROOT:-$ROOT/experiments/candidate_stage/cluster_results}
ASSEMBLY_ROOT=${ASSEMBLY_ROOT:-/home/bigssd/zt/cami2}
GOLD_ROOT=${GOLD_ROOT:-/home/bigssd/zt/runs/cami2_benchmark}
export OPENBLAS_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS
read -r -a methods <<< "${METHODS:-concat whiten cca}"
run_label=${RUN_LABEL:-representation}

if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(strain_madness marine plant_associated); fi
for dataset in "${datasets[@]}"; do
  out="$OUTROOT/$dataset"
  mkdir -p "$out"
  cache="$BASE_ROOT/$dataset/baseline.cache"
  features="$out/features"
  if [[ "${SKIP_FEATURES:-0}" != 1 ]]; then
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$out/features.time" \
    "$HERE/build/k4_features" \
    "$ASSEMBLY_ROOT/$dataset/CAMI2_${dataset}_GoldStandardAssembly.fasta" \
    "$features" 1000 "$THREADS" > "$out/features.log" 2>&1
  fi
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$out/$run_label.time" \
    "$PY" "$HERE/joint_space.py" --cache "$cache" --features "$features" \
    --retained "$BASE_ROOT/$dataset/retained.rbedge" --output "$out" \
    --threads "$THREADS" --neighbours 200 --methods "${methods[@]}" > "$out/$run_label.log" 2>&1
  for method in "${methods[@]}"; do
    for stage in candidates weights joint; do
      prefix="$out/$method.$stage"
      if [[ "$stage" == weights ]]; then
        graph="$out/$method.retained.rbedge"
      else
        graph="$out/$method.joint.rbedge"
      fi
      extra=()
      if [[ "$stage" == candidates ]]; then extra=(--external-graph-coverage); fi
      /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
        "$RB" bin --load-cache "$cache" --external-graph "$graph" "${extra[@]}" \
        --output "$prefix" --threads "$THREADS" --seed 42 \
        --sketch-m 500 --max-edges 200 --dual-depth 5 \
        --min-bin-size 200000 --no-bin-fasta \
        --export-retained-graph "$prefix.scored.rbedge" > "$prefix.log" 2>&1
      "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
        --members "$prefix.members.tsv" --min-bin-size 200000 \
        --threads "$THREADS" --quiet > "$prefix.amber" 2>&1
    done
  done
done
