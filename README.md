# RabbitBin

Fast, sketch-based metagenome binning. The default RabbitBin pipeline uses
4-mer weighted ProbMinHash (PMH) sketches to construct a mutual-nearest-neighbour
candidate graph, weights candidate edges by Spearman abundance correlation, and
clusters the graph with Fisher label propagation.

| Command | What it does |
|---------|--------------|
| `rabbitbin bin`   | Bin contigs into genomes (the main pipeline) |
| `rabbitbin depth` | Turn sorted BAM(s) into a MetaBAT/JGI depth TSV |
| `rabbitbin amber` | Fast, multithreaded AMBER-compatible binning evaluation |

`rabbitbin bin` is the default: a bare `rabbitbin -a contigs.fa -o out` still works.

## Requirements

- C++17 compiler with **OpenMP** (GCC ≥ 7 recommended)
- **CMake** ≥ 3.5
- **Boost** ≥ 1.66 (`program_options filesystem system graph serialization iostreams`)
- **zlib** ≥ 1.2.11 and **HTSlib** ≥ 1.13 — auto-downloaded if not found on the system

## Build

```bash
mkdir build && cd build
cmake ..
make rabbitbin -j
# binary: build/src/rabbitbin
```

After pulling new changes, rebuild from a clean tree:

```bash
rm -rf build && mkdir build && cd build && cmake .. && make -j
```

## Usage

### 1. Bin from BAMs in one shot (depth computed internally)

The pipeline takes coordinate-sorted BAMs directly; a `.bai` index is optional.

```bash
rabbitbin bin \
  --fasta contigs.fa \
  --bam-list bams.txt \      # one sorted-BAM path per line (or repeat --bam s1.bam s2.bam)
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

### 3. BAM(s) → depth TSV

```bash
rabbitbin depth --bam-list bams.txt --out depth.tsv --threads 64
rabbitbin depth --fasta contigs.fa --bam s1.bam s2.bam -o depth.tsv
```

### 4. Evaluate against a gold standard (AMBER-compatible)

```bash
rabbitbin amber \
  --gold gsa_mapping.binning \           # CAMI bioboxes, needs _LENGTH
  --members results/out.members.tsv \    # or --binning preds.binning (bioboxes / 2-col)
  --output metrics_per_bin.tsv \
  --threads 64
```

## Default algorithm

RabbitBin's default binning path contains five stages:

1. **4-mer PMH sketching.** Each contig is represented by a weighted
   ProbMinHash sketch over canonical 4-mer counts.
2. **500-register signatures.** Each PMH sketch contains 500 registers.
3. **Mutual top-N candidate graph.** PMH similarity retrieves at most 200
   neighbours per contig. An undirected candidate edge is retained only when
   the neighbour relation is mutual.
4. **Spearman edge weighting.** For multi-sample data, PMH is used only for
   candidate generation. Surviving edges are weighted by the Spearman
   correlation between contig abundance profiles.
5. **Fisher label propagation.** Label propagation combines the evidence from
   edges incident to each candidate label using Fisher's method and iterates
   until the partition converges.

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
| `-m, --min-contig` | 2500 | Minimum contig length (≥1500) |
| `-s, --min-bin-size` | 200000 | Minimum output bin size (bp) |
| `--min-edge-score` | 70 | Minimum Spearman edge weight (1–99) |
| `--max-edges` | 200 | Maximum PMH neighbours per contig before mutual filtering |
| `--sketch-m` | 500 | Number of ProbMinHash registers |
| `--percent-identity` | 97 | Min read identity when reading BAMs |
| `--bioboxes` | off | Also write a CAMI bioboxes `<prefix>.binning` |
| `--bin-fasta` | off | Also write per-bin FASTA files |
| `--unbinned` | off | Write unbinned contigs to FASTA |

Run `rabbitbin <command> --help` for the full option list.

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

## Pipeline wrapper

`run_rabbitbin.sh` runs BAM depth summarization then RabbitBin in one call:

```bash
run_rabbitbin.sh assembly.fa sample1.bam sample2.bam
```

## License

RabbitBin is released under the LBNL BSD license. Portions of the graph-clustering
pipeline derive from earlier open-source metagenome binning work; see `license.txt`.
