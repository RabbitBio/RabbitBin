#!/usr/bin/env bash
# Optional installed-binary check: Bash/coreutils/awk/gzip only, no Python/samtools.
set -euo pipefail
if [[ $# != 1 || ! -x "$1" ]]; then
  printf 'Usage: bash %s /path/to/rabbitbin\n' "$0" >&2
  exit 2
fi
binary=$(cd -- "$(dirname -- "$1")" && pwd)/${1##*/}
fixtures=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
while IFS='=' read -r key value; do
  case "$key" in
    RB_*|RABBIT_*|OMP_*|BADMAP|PCTID|MINDEPTH|LD_LIBRARY_PATH|LD_PRELOAD|DYLD_LIBRARY_PATH) unset "$key" ;;
  esac
done < <(env)
work=$(mktemp -d "${TMPDIR:-/tmp}/rabbitbin-core-smoke.XXXXXXXX")
trap 'printf "Smoke-test files retained at: %s\n" "$work"' EXIT
cd -- "$work"
"$binary" --version
"$binary" --help >/dev/null
cp "$fixtures/contigs.fa" 'assembly space.fa'
cp "$fixtures/contigs-1000.fastq.bam" 'sample one.bam'
cp "$fixtures/contigs-1000.fastq.bam" 'sample two.bam'
gzip -c 'assembly space.fa' > 'assembly space.fa.gz'
printf '%s\n' "$work/sample one.bam" "$work/sample two.bam" > 'bams list.txt'

# The fixture contains one ~5 kb contig: lower only the output-size cutoff.
# Reusing the BAM as a second sample exercises list parsing, not biology.
common=(--threads 2 --min-bin-size 0 --seed 42 --bin-fasta)
"$binary" bin --assembly 'assembly space.fa' --bam 'sample one.bam' \
  --output bam "${common[@]}"
"$binary" bin --assembly 'assembly space.fa' --bam-list 'bams list.txt' \
  --output bam-list "${common[@]}"
"$binary" bin --assembly 'assembly space.fa.gz' --depth "$fixtures/contigs_depth.txt" \
  --output gz-depth "${common[@]}"
"$binary" depth --bam-list 'bams list.txt' --out generated.depth.tsv --threads 2
"$binary" bin --assembly 'assembly space.fa' --depth generated.depth.tsv \
  --output generated-depth "${common[@]}"

sequence=$(awk 'NR==1 {sub(/^>/, "", $1); print $1; exit}' 'assembly space.fa')
for prefix in bam bam-list gz-depth generated-depth; do
  # BAM and ordinary depth TSV inputs have different coverage features; test
  # membership and positive coverage, not byte identity of TotalDepth.
  awk -F '\t' -v sequence="$sequence" '
    NR==1 {
      for (i=1; i<=NF; ++i) column[$i]=i
      if (!column["BinNum"] || !column["SequenceName"] || !column["TotalDepth"]) exit 1
    }
    NR==2 {
      if ($(column["BinNum"]) != 1 || $(column["SequenceName"]) != sequence ||
          $(column["TotalDepth"]) <= 0) exit 1
    }
    END {if (NR != 2) exit 1}
  ' "$prefix.members.tsv"
  test -s "$prefix.bins.tsv"
  test -s "${prefix}_bin_001.fa"
done
printf 'Core smoke tests passed (single BAM, BAM list, gzip FASTA, and depth TSV).\n'
