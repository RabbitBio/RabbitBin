#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
RB=${RB:-$ROOT/build/src/rabbitbin}
PY=${PY:-$ROOT/experiments/candidate_stage/venv/bin/python}
ASSEMBLY_ROOT=${ASSEMBLY_ROOT:-/home/bigssd/zt/cami2}
BAM_ROOT=${BAM_ROOT:-/home/bigssd/zt/CAMI2_remap_bams}
GOLD_ROOT=${GOLD_ROOT:-/home/bigssd/zt/runs/cami2_benchmark}
THREADS=${THREADS:-64}
OUTROOT=${OUTROOT:-$ROOT/experiments/candidate_stage/cluster_results}

if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(marine plant_associated strain_madness); fi
for dataset in "${datasets[@]}"; do
  mkdir -p "$OUTROOT/$dataset"
  base="$OUTROOT/$dataset/baseline"
  graph="$OUTROOT/$dataset/retained.rbedge"
  cache="$OUTROOT/$dataset/baseline.cache"
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$base.time" \
    "$RB" bin \
    --assembly "$ASSEMBLY_ROOT/$dataset/CAMI2_${dataset}_GoldStandardAssembly.fasta" \
    --bam-list "$BAM_ROOT/$dataset/bam.list" \
    --output "$base" --threads "$THREADS" --seed 42 \
    --sketch-m 500 --max-edges 200 --dual-depth 5 \
    --min-bin-size 200000 --no-bin-fasta \
    --save-cache "$cache" --export-retained-graph "$graph" \
    > "$base.log" 2>&1
  "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
    --members "$base.members.tsv" --min-bin-size 200000 \
    --threads "$THREADS" --quiet > "$base.amber" 2>&1

  for method in infomap leiden; do
    labels="$OUTROOT/$dataset/$method.labels.tsv"
    prefix="$OUTROOT/$dataset/$method"
    /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.cluster.time" \
      "$PY" "$ROOT/experiments/candidate_stage/cluster_graph.py" \
      "$graph" --method "$method" --output "$labels" \
      > "$prefix.cluster.log" 2>&1
    /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.rebin.time" \
      "$RB" bin --load-cache "$cache" --external-labels "$labels" \
      --output "$prefix" --threads "$THREADS" --seed 42 \
      --sketch-m 500 --max-edges 200 --dual-depth 5 \
      --min-bin-size 200000 --no-bin-fasta > "$prefix.log" 2>&1
    "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
      --members "$prefix.members.tsv" --min-bin-size 200000 \
      --threads "$THREADS" --quiet > "$prefix.amber" 2>&1
  done
done
