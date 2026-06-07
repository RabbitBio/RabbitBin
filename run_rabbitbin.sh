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

if ! $RB --help 2>/dev/null; then
  echo "Please ensure RabbitBin is in PATH: could not find $RB" 1>&2
  exit 1
fi

if ! $SUM 2>/dev/null; then
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

rbopts=""
for arg in "$@"; do
  if [ -f "$arg" ]; then
    break
  fi
  rbopts="$rbopts $arg"
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

waitforlock() {
  while [ -L "$lock" ]; do
    echo "Waiting for $lock ($(date))"
    sleep 60
  done
}

cleanup() { rm -f $lock; }
trap cleanup 0 1 2 3 15

waitforlock
badmap=${assembly##*/}.d
badmapopts=
if [ "$BADMAP" != '0' ]; then
  mkdir -p "${badmap}"
  badmapopts="--unmappedFastq ${badmap}/badmap"
fi

if [ ! -f "${depth}" ] && ln -s "$(uname -n) $$" $lock; then
  if [ ! -f "${depth}" ]; then
    sumopts="--outputDepth ${depth}.tmp --percentIdentity ${PCTID} --minContigLength 1000 --minContigDepth ${MINDEPTH} ${badmapopts} --referenceFasta ${assembly}"
    echo "Running depth: $SUM $sumopts $* ($(date))"
    $SUM $sumopts "$@" && mv ${depth}.tmp ${depth}
  fi
  rm -f $lock
else
  waitforlock
  echo "Using existing depth file: $depth"
fi

outname=${assembly##*/}.rabbitbin-${RB_LABEL}-$(date '+%Y%m%d_%H%M%S')/out
echo "Running RabbitBin: $RB $rbopts -a $assembly -o $outname -d ${depth} ($(date))"
$RB $rbopts --assembly "$assembly" --output "$outname" --depth "${depth}"
echo "Finished RabbitBin ($(date))"
