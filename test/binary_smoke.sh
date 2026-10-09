#!/usr/bin/env bash
# Requires only Bash/coreutils/gzip and the extracted binary package.
set -euo pipefail
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
while IFS='=' read -r key value; do
  case "$key" in
    RB_*|RABBIT_*|OMP_*|BADMAP|PCTID|MINDEPTH|LD_LIBRARY_PATH|LD_PRELOAD) unset "$key" ;;
  esac
done < <(env)
work=$(mktemp -d "${TMPDIR:-/tmp}/rabbitbin-binary-test.XXXXXXXX")
cd -- "$work"
version=$(sed -n 's/^version=//p' "$package/BUILDINFO.txt")
revision=$(sed -n 's/^commit=//p' "$package/BUILDINFO.txt")
test "$("$package/bin/rabbitbin" --version)" = "RabbitBin $version (commit ${revision:0:12})"
"$package/bin/rabbitbin" --help >/dev/null
cp "$package/test/contigs.fa" 'assembly space.fa'
cp "$package/test/contigs-1000.fastq.bam" 'sample one.bam'
gzip -c 'assembly space.fa' > 'assembly space.fa.gz'
"$package/bin/rabbitbin" bin --assembly 'assembly space.fa' --bam 'sample one.bam' \
  --output bam --min-bin-size 0 --threads 2 --seed 42
"$package/bin/rabbitbin" bin --assembly 'assembly space.fa.gz' \
  --depth "$package/test/contigs_depth.txt" --output gz --min-bin-size 0 --threads 2 --seed 42
RB_LABEL='binary test' "$package/bin/run_rabbitbin.sh" \
  --min-bin-size 0 --threads 2 --seed 42 'assembly space.fa' 'sample one.bam'
for members in bam.members.tsv gz.members.tsv ./*.rabbitbin-*/out.members.tsv; do
  test "$(wc -l < "$members")" -eq 2
  awk -F '\t' 'NR==2 {if ($1 != 1 || $NF <= 0) exit 1}' "$members"
done
printf 'Binary smoke tests passed: %s\n' "$work"
