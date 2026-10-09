#!/bin/bash
# RabbitBin convenience wrapper: BAM -> depth TSV -> binning

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
PATH=$SCRIPTPATH:$PATH
RB=rabbitbin
SUM=rabbit_depth
BADMAP=${BADMAP:=0}
PCTID=${PCTID:=97}
MINDEPTH=${MINDEPTH:=1.0}
RB_LABEL=${RB_LABEL:="bins"}

if ! "$RB" --help >/dev/null 2>&1; then
  echo "Please ensure RabbitBin is in PATH: could not find $RB" 1>&2
  exit 1
fi

if ! "$SUM" 2>/dev/null; then
  echo "Please ensure depth utility is in PATH: could not find $SUM" 1>&2
  exit 1
fi

USAGE="$0 <rabbitbin options> assembly.fa sample1.bam [ sample2.bam ...]
Pass any rabbitbin options except:
  -a/--assembly  -o/--output  -d/--depth

Depth-stage environment variables:
  PCTID=${PCTID}       discard reads below this %% identity
  BADMAP=${BADMAP}     write discarded reads to a subdirectory
  MINDEPTH=${MINDEPTH}  minimum contig depth to emit
  RB_LABEL=${RB_LABEL} label in output directory name

Full options: $RB --help
"

rbopts=()
for arg in "$@"; do
  if [ -f "$arg" ]; then
    break
  fi
  rbopts+=("$arg")
  shift
done

if [ $# -lt 2 ]; then
  echo "$USAGE" 1>&2
  exit 1
fi

assembly=$1
shift
if [ ! -f "$assembly" ]; then
  echo "Assembly not found: $assembly" 1>&2
  exit 1
fi

for bam in "$@"; do
  if [ ! -f "$bam" ]; then
    echo "BAM not found: $bam" 1>&2
    exit 1
  fi
done

set -e

depth=${assembly##*/}.depth.txt
lock=$depth.BUILDING
owns_lock=0

waitforlock() {
  while [ -L "$lock" ]; do
    echo "Waiting for $lock ($(date))"
    sleep 60
  done
}

cleanup() {
  if [ "$owns_lock" = 1 ]; then
    rm -f -- "$lock"
  fi
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

badmap=${assembly##*/}.d
badmapopts=()
if [ "$BADMAP" != '0' ]; then
  mkdir -p "${badmap}"
  badmapopts=(--unmappedFastq "${badmap}/badmap")
fi

while :; do
  waitforlock
  if [ -f "$depth" ]; then
    echo "Using existing depth file: $depth"
    break
  fi
  if ln -s -- "$(uname -n) $$" "$lock"; then
    owns_lock=1
    # Another process may have finished between our existence check and lock.
    if [ ! -f "$depth" ]; then
      sumopts=(--outputDepth "${depth}.tmp" --percentIdentity "$PCTID"
               --minContigLength 1000 --minContigDepth "$MINDEPTH"
               "${badmapopts[@]}" --referenceFasta "$assembly")
      printf 'Running depth:'
      printf ' %q' "$SUM" "${sumopts[@]}" "$@"
      printf '\n'
      if ! "$SUM" "${sumopts[@]}" "$@"; then
        echo "Depth generation failed; binning was not started." >&2
        exit 1
      fi
      mv -- "${depth}.tmp" "$depth"
    fi
    cleanup
    owns_lock=0
    break
  fi
  # A concurrent writer owns a symlink lock; wait and retry. Other failures
  # (permissions, a non-lock file at this path) must not masquerade as a cache.
  if [ ! -L "$lock" ]; then
    echo "Cannot acquire depth lock: $lock" >&2
    exit 1
  fi
done

outname=${assembly##*/}.rabbitbin-${RB_LABEL}-$(date '+%Y%m%d_%H%M%S')/out
printf 'Running RabbitBin:'
printf ' %q' "$RB" "${rbopts[@]}" --assembly "$assembly" --output "$outname" --depth "$depth"
printf '\n'
"$RB" "${rbopts[@]}" --assembly "$assembly" --output "$outname" --depth "$depth"
echo "Finished RabbitBin ($(date))"
