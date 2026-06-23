# RabbitBin

Fast, sketch-based metagenome binning. RabbitBin clusters assembly contigs into
genome bins (MAGs) from a composition–abundance similarity graph, with optional
abundance-guided bin splitting. It ships as a single binary with three commands:

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

### 4. Fast re-binning with a cache

Build the expensive feature graph once, then sweep cheap parameters in seconds.
`--load-cache` skips FASTA parse / sketch / depth / graph build; only lightweight
parameters (`--min-bin-size`, `--split-silhouette`, `--min-edge-score`, recruit…)
take effect.

```bash
# Build once
rabbitbin bin --fasta contigs.fa --depth depth.tsv --save-cache run.cache -o tmp

# Re-bin instantly under different criteria
rabbitbin bin --load-cache run.cache --min-bin-size 200000 -o out_200k
rabbitbin bin --load-cache run.cache --split-silhouette 0.7 -o out_sil07
```

### 5. Evaluate against a gold standard (AMBER-compatible)

```bash
rabbitbin amber \
  --gold gsa_mapping.binning \           # CAMI bioboxes, needs _LENGTH
  --members results/out.members.tsv \    # or --binning preds.binning (bioboxes / 2-col)
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
| `-m, --min-contig` | 2500 | Minimum contig length (≥1500) |
| `-s, --min-bin-size` | 50000 | Minimum output bin size (bp) — *cheap, cache-friendly* |
| `--min-edge-score` | 60 | Minimum graph edge weight (1–99) — *cheap, cache-friendly* |
| `--split-silhouette` | env | Abundance-split silhouette threshold — *cheap, cache-friendly* |
| `--no-split` | off | Disable abundance bin splitting |
| `--no-recruit` | off | Disable small-contig recruiting |
| `--sim-cutoff` | 0 (auto) | Composition similarity cutoff (×100) |
| `--sketch-k` / `--sketch-m` | 8 / 500 | Sketch k-mer size / ProbMinHash registers |
| `--percent-identity` | 97 | Min read identity when reading BAMs |
| `--save-cache` / `--load-cache` | — | Write / re-bin from a feature-graph cache |
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

## Key options (`amber`)

| Option | Default | Meaning |
|--------|---------|---------|
| `-g, --gold` | — | Gold-standard binning (CAMI bioboxes, needs `_LENGTH`) |
| `-i, --binning` | — | Predicted binning (bioboxes or 2-col `SEQ<TAB>BIN`) |
| `--members` | — | Predicted binning as rabbitbin `members.tsv` |
| `-o, --output` | — | Per-bin metrics TSV (optional) |
| `--min-length` | 0 | Ignore GS contigs shorter than this |
| `-q, --quiet` | off | Print only the summary |

## Environment variables (advanced tuning)

- `RABBIT_PMH` — weighted ProbMinHash graph (default on)
- `RABBIT_MUTUAL_KNN` — mutual k-NN graph (default on)
- `RABBIT_GC_NORM` — GC normalization mode
- `RABBIT_W_COMP` — composition vs. abundance edge weight
- `RABBIT_NEG_DEPTH` — negative depth-correlation filter (default −0.3)
- `RABBIT_SPLIT_SIL` — default silhouette split threshold

## Pipeline wrapper

`run_rabbitbin.sh` runs BAM depth summarization then RabbitBin in one call:

```bash
run_rabbitbin.sh assembly.fa sample1.bam sample2.bam
```

## License

RabbitBin is released under the LBNL BSD license. Portions of the graph-clustering
pipeline derive from earlier open-source metagenome binning work; see `license.txt`.
