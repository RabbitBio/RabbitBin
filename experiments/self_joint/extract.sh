#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
RB=${RB:-$ROOT/build/src/rabbitbin}
OUTROOT=${OUTROOT:-$ROOT/experiments/self_joint/results}
BASE_ROOT=${BASE_ROOT:-$ROOT/experiments/candidate_stage/cluster_results}
ASSEMBLY_ROOT=${ASSEMBLY_ROOT:-/home/bigssd/zt/cami2}
BAM_ROOT=${BAM_ROOT:-/home/bigssd/zt/CAMI2_remap_bams}
K4=${K4:-$ROOT/build/k4_features}
THREADS=${THREADS:-64}
if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(marine plant_associated strain_madness); fi
for dataset in "${datasets[@]}"; do
  out="$OUTROOT/$dataset"
  mkdir -p "$out"
  assembly="$ASSEMBLY_ROOT/$dataset/CAMI2_${dataset}_GoldStandardAssembly.fasta"
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$out/extract.time" \
    env RB_TIMING=1 RB_DEPTH_PROF=1 "$RB" bin --assembly "$assembly" \
    --bam-list "$BAM_ROOT/$dataset/bam.list" \
    --output "$out/baseline_export" --threads "$THREADS" --seed 42 \
    --sketch-m 500 --max-edges 200 --dual-depth 5 --min-bin-size 200000 --no-bin-fasta \
    --export-fragment-depth "$out/halves" > "$out/extract.log" 2>&1
  cmp "$out/baseline_export.members.tsv" \
      "$BASE_ROOT/$dataset/baseline.members.tsv"
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$out/composition.time" \
    "$K4" "$assembly" "$out/halves" 2500 "$THREADS" --halves \
    > "$out/composition.log" 2>&1
  printf '%s fragment extraction done; baseline members identical\n' "$dataset"
done
