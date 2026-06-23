#!/usr/bin/env bash
# rabbitbin_markers.sh — build a contig->marker map for `rabbitbin qc`.
#
# This is the EXPENSIVE, run-once-per-assembly step: predict genes (Prodigal)
# and search a single-copy marker-gene HMM set (hmmsearch). The resulting map is
# reused to score any number of binnings cheaply with `rabbitbin qc`.
#
# Usage:
#   rabbitbin_markers.sh <contigs.fa> <out.markers.tsv> [threads]
#
# Environment overrides:
#   RABBITBIN_MARKER_HMM   path to marker HMM   (default: ../data/markers.hmm)
#   PRODIGAL               prodigal binary      (default: auto-detect)
#   HMMSEARCH              hmmsearch binary      (default: auto-detect)
#   PRODIGAL_MODE          prodigal -p mode      (default: meta)
#
# Output: tab-separated `contig <TAB> marker`, one gene hit per line, prefixed by
#   `# marker_set_size=<N>` so `rabbitbin qc` knows the marker-set size.
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <contigs.fa> <out.markers.tsv> [threads]" >&2
  exit 1
fi

contigs=$1
out=$2
threads=${3:-8}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HMM=${RABBITBIN_MARKER_HMM:-"$SCRIPT_DIR/../data/markers.hmm"}

# Resolve prodigal / hmmsearch: PATH first, then common conda envs on this host.
find_tool() {
  local name=$1 envvar=$2 val
  val=$(eval "echo \${$envvar:-}")
  if [[ -n "$val" ]]; then echo "$val"; return; fi
  if command -v "$name" >/dev/null 2>&1; then command -v "$name"; return; fi
  local c
  for c in "$HOME"/anaconda3/envs/*/bin/"$name" \
           "$HOME"/miniconda3/envs/*/bin/"$name"; do
    [[ -x "$c" ]] && { echo "$c"; return; }
  done
  echo ""
}

PRODIGAL=$(find_tool prodigal PRODIGAL)
HMMSEARCH=$(find_tool hmmsearch HMMSEARCH)
PRODIGAL_MODE=${PRODIGAL_MODE:-meta}

[[ -z "$PRODIGAL"  ]] && { echo "[Error] prodigal not found (set PRODIGAL=/path/to/prodigal)" >&2; exit 1; }
[[ -z "$HMMSEARCH" ]] && { echo "[Error] hmmsearch not found (set HMMSEARCH=/path/to/hmmsearch)" >&2; exit 1; }
[[ -f "$HMM" ]] || { echo "[Error] marker HMM not found: $HMM" >&2; exit 1; }

echo "[markers] prodigal : $PRODIGAL" >&2
echo "[markers] hmmsearch: $HMMSEARCH" >&2
echo "[markers] hmm      : $HMM" >&2

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
faa="$tmp/proteins.faa"
tbl="$tmp/hits.tbl"

# Prodigal's meta mode scores each contig independently, so splitting the
# assembly into N chunks and running Prodigal in parallel is exact and gives a
# near-linear speedup (Prodigal itself is single-threaded). N = threads.
echo "[markers] (1/3) predicting genes with Prodigal (-p $PRODIGAL_MODE, $threads-way parallel)..." >&2
chunkdir="$tmp/chunks"
mkdir -p "$chunkdir"
awk -v n="$threads" -v dir="$chunkdir" '
  /^>/ { f = dir "/chunk." (++c % n) ".fa" }
  { print > f }
' "$contigs"

ls "$chunkdir"/chunk.*.fa 2>/dev/null | xargs -P "$threads" -I{} \
  "$PRODIGAL" -p "$PRODIGAL_MODE" -q -a {}.faa -i {} -o /dev/null
cat "$chunkdir"/chunk.*.fa.faa > "$faa"

echo "[markers] (2/3) hmmsearch against marker set (--cut_tc, $threads threads)..." >&2
"$HMMSEARCH" --cpu "$threads" --cut_tc --noali --tblout "$tbl" "$HMM" "$faa" >/dev/null

echo "[markers] (3/3) writing contig->marker map: $out" >&2
G=$(grep -c '^NAME' "$HMM")
{
  echo "# marker_set_size=$G"
  # tblout: $1=protein (contig_geneNum), $3=marker (query name).
  # Strip the trailing _<geneNum> Prodigal appends to recover the contig id.
  awk 'BEGIN{OFS="\t"} !/^#/ { prot=$1; marker=$3; sub(/_[0-9]+$/,"",prot); print prot, marker }' "$tbl"
} > "$out"

nhits=$(grep -vc '^#' "$out" || true)
echo "[markers] done: $nhits marker hits, marker_set_size=$G -> $out" >&2
