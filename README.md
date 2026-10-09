# RabbitBin

Fast, sketch-based metagenome binning. The default RabbitBin pipeline uses
canonical 4-mer ProbMinHash (PMH) sketches weighted by enrichment relative to
each contig's own base composition to construct a
bounded mutual-nearest-neighbour candidate graph, uses abundance profiles as
the edge evidence when coverage is available, clusters the retained
graph with Fisher label propagation, and re-splits multi-modal bins. One final
coverage-based pass recruits all remaining long and short contigs against the
split, frozen bin cores. In other words, PMH proposes where to look; it does not
by itself determine the final biological grouping.

RabbitBin uses [RabbitBAM](https://github.com/RabbitBio/RabbitBAM/tree/sortedbam)
for parallel BAM I/O. The `sortedbam` branch provides the BAM reading, sorting,
and indexing modules.

| Command | What it does |
|---------|--------------|
| `rabbitbin bin`    | Bin contigs into genomes (the main pipeline) |
| `rabbitbin depth`  | Turn sorted BAM(s)/CRAM into a MetaBAT/JGI depth TSV |
| `rabbitbin qc`     | Score a binning by SCG completeness/contamination (no gold standard) |
| `rabbitbin refine` | DAS Tool-style SCG consensus over several independent binnings |
| `rabbitbin amber`  | Fast, multithreaded AMBER-compatible binning evaluation |

`rabbitbin bin` is the default: a bare `rabbitbin -a contigs.fa -o out` still works.
Read-mapping subcommands (`map`, `bwa`, `sortbam`, `bai`) exist only in builds
configured with `-DRABBITBIN_ENABLE_MAP=ON`.

## Requirements

- C++17 compiler with **OpenMP** (GCC ≥ 7 or a recent Clang)
- **CMake** ≥ 3.16
- **Boost** ≥ 1.66 (`program_options filesystem system graph serialization iostreams regex`)
- `git`, `make`, and Autotools when dependencies must be downloaded
- **zlib** ≥ 1.2.11, **HTSlib** ≥ 1.13, and **libdeflate** — pinned copies
  are downloaded automatically when development packages are unavailable
- Python 3.6+ and `samtools` for the complete test suite (not core binning)

## Build

On Ubuntu/Debian, the complete system-dependency route is:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git autoconf automake libtool \
  pkg-config libboost-all-dev zlib1g-dev libhts-dev libdeflate-dev
```

On Rocky/RHEL, zlib, HTSlib and libdeflate may instead be left to the automatic
dependency build:

```bash
sudo dnf install -y gcc gcc-c++ cmake git make autoconf automake libtool \
  pkgconf-pkg-config boost-devel
```

Then clone and build. Conda and Docker are not required.

```bash
git clone https://github.com/RabbitBio/RabbitBin.git
cd RabbitBin
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
(cd build && ctest --output-on-failure)
# binary: build/src/rabbitbin
```

An optional per-user installation is:

```bash
cmake --install build --prefix "$HOME/.local"
"$HOME/.local/bin/rabbitbin" --version
```

The default binary uses a portable CPU baseline. For a binary that will only
run on the machine where it is compiled, enable local CPU tuning explicitly:

```bash
cmake -S . -B build-native -DRABBITBIN_NATIVE_ARCH=ON
cmake --build build-native --parallel
```

Conda can be used as an optional dependency manager:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate rabbitbin
cmake -S . -B build-conda -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
cmake --build build-conda --parallel
(cd build-conda && ctest --output-on-failure)
cmake --install build-conda
rabbitbin --version
```

The environment file specifies both `conda-forge` and `bioconda` (HTSlib and
samtools), plus the full test dependencies. The explicit CMake prefix keeps
the build pointed at this environment. Use a fresh build directory when
switching environments so cached system-library paths are not reused.
See [INSTALL.md](INSTALL.md) for source archives and installation smoke tests.

After pulling new changes, configure a fresh build directory:

```bash
cmake -S . -B build-new -DCMAKE_BUILD_TYPE=Release
cmake --build build-new --parallel
```

## Input modes

`rabbitbin bin` accepts coverage in one of two forms, or none at all. The mode
is inferred from the flags and echoed in the log line beginning `Edge weighting:`.

| Mode | Flags | Edge weight | Abundance stages |
|------|-------|-------------|------------------|
| BAM/CRAM | `--fasta` + `--bam`/`--bam-list` | coverage only (depth computed in-process) | enabled |
| Precomputed depth | `--assembly` + `--depth` | coverage only | enabled |
| Sequence-only | `--assembly` alone | composition only | **disabled** |

PMH selects candidate neighbours and coverage determines their edge weights.
With one or two input samples, RabbitBin averages min/max ratios over
informative coverage features. With at least three input samples, it combines
Spearman correlation and weighted Jaccard; see the
[pipeline documentation](docs/pipeline.md). Coverage also drives splitting
and recruitment.

Here `S` counts input samples: one per BAM, or one per mean-depth column in a
precomputed table. BAM input additionally computes MAPQ-filtered depth features
by default (`--dual-depth 5`): each sample contributes ordinary coverage and
coverage from reads with MAPQ >= 5. Both features are used for graph weighting
and recruitment, giving two coverage dimensions for one BAM and four for two
BAMs. With `--dual-depth 0`, only ordinary coverage is used. Precomputed depth
tables supply one mean-depth feature per sample; variance columns are not
coverage dimensions. The additional BAM features do not increase `S`. The log
reports both the input sample count and the total number of coverage columns.

**Sequence-only mode runs, but it is not the configuration used for the reported
benchmarks.** Without coverage there is no abundance signal, so contig
recruitment and abundance-guided bin splitting are both skipped and edges are
weighted by PMH composition similarity alone. Use it only for assemblies with no
reads available; expect materially lower bin quality on multi-sample datasets.

## Usage

### 1. Bin from BAMs in one shot (depth computed internally)

The pipeline takes coordinate-sorted BAMs directly; a `.bai` index is optional.
The BAM-list file contains one sorted BAM or CRAM path per line; alternatively,
repeat `--bam` with one or more paths.

```bash
rabbitbin bin \
  --fasta contigs.fa \
  --bam-list bams.txt \
  --percent-identity 97 \
  --threads 64 \
  --output results/out
```

### 2. Bin from a precomputed depth file

```bash
rabbitbin bin \
  --assembly contigs.fa \
  --depth depth.tsv \
  --output results/out \
  --threads 32
```

### 3. Bin from the assembly alone (sequence-only)

```bash
rabbitbin bin --assembly contigs.fa --output results/out --threads 32
```

### 4. BAM(s) → depth TSV

```bash
rabbitbin depth --bam-list bams.txt --out depth.tsv --threads 64
rabbitbin depth --fasta contigs.fa --bam s1.bam s2.bam -o depth.tsv
```

### 5. Evaluate against a gold standard (AMBER-compatible)

```bash
# CAMI bioboxes gold mapping; the gold file must provide _LENGTH.
# Use --binning preds.binning instead of --members for a two-column prediction.
rabbitbin amber \
  --gold gsa_mapping.binning \
  --members results/out.members.tsv \
  --output metrics_per_bin.tsv \
  --threads 64
```

## Outputs (`bin`)

| File | Description |
|------|-------------|
| `prefix.members.tsv` | Contig-to-bin membership |
| `prefix.bins.tsv` | Per-bin stats |
| `prefix.unbinned.fa` | Unbinned contigs (with `--unbinned`) |
| `prefix_bin_001.fa` | Per-bin FASTA (only with `--bin-fasta`) |

## Key options (`bin`)

| Option | Default | Meaning |
|--------|---------|---------|
| `-a, --assembly` / `--fasta` | — | Input contig FASTA (gzip ok) |
| `-o, --output` | — | Output path prefix |
| `-d, --depth` | — | Coverage depth TSV (MetaBAT/JGI format) |
| `--bam` / `--bam-list` | — | Sorted BAM input (compute depth in-process) |
| `-t, --threads` | 0 (all) | Worker threads |
| `-m, --min-contig` | 2500 | Minimum length of a clustered contig (must be ≥1500) |
| `--min-small-contig` | 1000 | Minimum length of a recruitable short contig (must be ≥500); shorter contigs are discarded |
| `-s, --min-bin-size` | 200000 | Minimum output bin size (bp) |
| `--min-edge-score` | 71.53318629591614 | Minimum coverage edge weight, percent (>1 and <100, decimals accepted) |
| `--max-edges` | 200 | Maximum PMH neighbours per contig among production-feasible pairs, before mutual filtering |
| `--sketch-m` | 500 | Number of ProbMinHash registers |
| `--validate-pmh-gold` | — | Evaluate sequence-only PMH top-N neighbourhoods against CAMI gold, write TSV, and exit |
| `--validate-pmh-queries` | 1000 | Seed-controlled labelled queries used by PMH validation |
| `--validate-pmh-top` | 400 | Largest neighbourhood retained by PMH validation |
| `--audit-graph-gold` | — | Audit production candidate and abundance-retained edges against CAMI gold without changing binning |
| `--audit-graph-out` | `<output>.graph_audit.tsv` | Output TSV for `--audit-graph-gold` |
| `--no-recruit` | off | Disable the post-split long/short-contig coverage recruitment |
| `--recruit-max-fpr` | 0.05 | Maximum leave-one-out ROC false-positive rate used to select the recruitment threshold; use 1 for unconstrained Youden |
| `--no-singleton-rescue` | off | Disable promotion of output-sized unassigned long contigs |
| `--no-split` | off | Disable abundance-guided bin splitting |
| `--split-silhouette` | 0.70 | Minimum mean silhouette to accept a split |
| `--split-max-k` | 6 | Maximum sub-clusters per split bin |
| `--split-kmeans-restarts` | 10 | K-means initializations tested for each K |
| `--percent-identity` | 97 | Min read identity when reading BAMs |
| `--markers` | — | Contig→marker map, required by `--qc`/`--purify`/`--auto`/`--autotune` |
| `--qc` | off | Annotate `bins.tsv` with SCG completeness/contamination + MIMAG tier |
| `--bioboxes` | off | Also write a CAMI bioboxes `<prefix>.binning` |
| `--bin-fasta` | off | Also write per-bin FASTA files |
| `--unbinned` | off | Write unbinned contigs to FASTA |
| `--save-cache` / `--load-cache` | — | Cache the graph for fast re-binning |

Run `rabbitbin <command> --help` for the full option list, including the `qc`
and `refine` subcommands.

Marker-based options additionally require Prodigal, HMMER's `hmmsearch`, and a
single-copy marker HMM collection. Set `RABBITBIN_MARKER_HMM` to that HMM file,
then generate the reusable map with
`rabbitbin_markers.sh contigs.fa contigs.markers.tsv 32`.

## Key options (`depth`)

| Option | Default | Meaning |
|--------|---------|---------|
| `--bam` / `--bam-list` | — | Sorted BAM input (required) |
| `-o, --out` | stdout | Output depth TSV |
| `--percent-identity` | 97 | Min mapped-read percent identity |
| `--min-contig-length` | 1 | Min contig length emitted |
| `--max-edge-bases` | 75 | Bases trimmed per contig end |
| `--no-variance` | off | Omit per-sample variance columns |
| `--long-read` | off | Long-read preset (percent-identity default 80) |
| `--reference` | — | Reference FASTA (required for CRAM input) |

## Key options (`amber`)

| Option | Default | Meaning |
|--------|---------|---------|
| `-g, --gold` | — | Gold-standard binning (CAMI bioboxes, needs `_LENGTH`) |
| `-i, --binning` | — | Predicted binning (bioboxes or 2-col `SEQ<TAB>BIN`) |
| `--members` | — | Predicted binning as rabbitbin `members.tsv` |
| `-o, --output` | — | Per-bin metrics TSV (optional) |
| `--min-length` | 0 | Ignore GS contigs shorter than this |
| `-q, --quiet` | off | Print only the summary |

## Reproducibility

`--seed` defaults to `0`, which seeds the RNG from the wall clock. The k-means
restarts in the bin-splitting stage consume that RNG, so two runs on identical
input can differ by a small number of contigs. Pass an explicit `--seed` for any
run you intend to report:

```bash
rabbitbin bin --assembly contigs.fa --depth depth.tsv \
              --output results/out --threads 64 --seed 1
```

The other stages described in the [default pipeline](docs/pipeline.md) are
deterministic given the input files and thread count. No other flags were used for the published
benchmarks.

## Pipeline wrapper

`run_rabbitbin.sh` runs BAM depth summarization then RabbitBin in one call:

```bash
run_rabbitbin.sh assembly.fa sample1.bam sample2.bam
```

## Optional Docker build

Docker is not required. It is provided as an alternative reproducible build:

```bash
docker build \
  --build-arg RABBITBIN_SOURCE_REVISION="$(git rev-parse --short=12 HEAD)" \
  -t rabbitbin:local .
docker run --rm rabbitbin:local rabbitbin --version
```

Mount input and output directories when running analyses, for example
`-v "$PWD:/work" -w /work`.

## Further documentation

- [Installation and troubleshooting](INSTALL.md)
- [Default pipeline and PMH validation](docs/pipeline.md)
- [Performance diagnostics and optimization checks](docs/performance.md)

## License

RabbitBin is distributed under the BSD 3-Clause license in `LICENSE`. Portions
derive from earlier BSD-licensed metagenome-binning work; attribution and
third-party notices are recorded in `license.txt`.
