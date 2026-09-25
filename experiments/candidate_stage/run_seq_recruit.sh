#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
RB=${RB:-$ROOT/build/src/rabbitbin}
ASSEMBLY_ROOT=${ASSEMBLY_ROOT:-/home/bigssd/zt/cami2}
BAM_ROOT=${BAM_ROOT:-/home/bigssd/zt/CAMI2_remap_bams}
GOLD_ROOT=${GOLD_ROOT:-/home/bigssd/zt/runs/cami2_benchmark}
THREADS=${THREADS:-64}
OUTROOT=${OUTROOT:-$ROOT/experiments/candidate_stage/seq_results}
MODE=${MODE:-plain}
if [[ "$MODE" == agree ]]; then
  seq_args=(--seq-recruit-agree)
  label=seq_recruit_agree
elif [[ "$MODE" == plain ]]; then
  seq_args=(--seq-recruit)
  label=seq_recruit
else
  echo "MODE must be plain or agree" >&2
  exit 2
fi

if [[ $# -gt 0 ]]; then datasets=("$@"); else datasets=(marine plant_associated strain_madness); fi
for dataset in "${datasets[@]}"; do
  mkdir -p "$OUTROOT/$dataset"
  prefix="$OUTROOT/$dataset/$label"
  /usr/bin/time -f 'wall_s=%e peak_kb=%M' -o "$prefix.time" \
    "$RB" bin \
    --assembly "$ASSEMBLY_ROOT/$dataset/CAMI2_${dataset}_GoldStandardAssembly.fasta" \
    --bam-list "$BAM_ROOT/$dataset/bam.list" \
    --output "$prefix" --threads "$THREADS" --seed 42 \
    --sketch-m 500 --max-edges 200 --dual-depth 5 \
    --min-bin-size 200000 --no-bin-fasta "${seq_args[@]}" \
    > "$prefix.log" 2>&1
  "$RB" amber --gold "$GOLD_ROOT/$dataset/prep/gold_len.binning" \
    --members "$prefix.members.tsv" --min-bin-size 200000 \
    --threads "$THREADS" --quiet > "$prefix.amber" 2>&1
done
