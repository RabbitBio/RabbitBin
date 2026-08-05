# RabbitBin

Fast, sketch-based metagenome binning. The default RabbitBin pipeline uses
4-mer weighted ProbMinHash (PMH) sketches to construct a mutual-nearest-neighbour
candidate graph, weights candidate edges by contig composition and coverage
agreement, clusters the graph with Fisher label propagation, and then recruits
short contigs and re-splits multi-modal bins.

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

## Input modes

`rabbitbin bin` accepts coverage in one of two forms, or none at all. The mode
is inferred from the flags and echoed in the log line beginning `Fusion:`.

| Mode | Flags | Edge weight | Abundance stages |
|------|-------|-------------|------------------|
| BAM/CRAM | `--fasta` + `--bam`/`--bam-list` | composition + coverage (depth computed in-process) | enabled |
| Precomputed depth | `--assembly` + `--depth` | composition + coverage | enabled |
| Sequence-only | `--assembly` alone | composition only | **disabled** |

With one or two coverage samples, edge weights fall back to composition alone,
but the coverage values still drive recruitment and bin splitting. From three
samples on, coverage drives the edge weights (see *Edge weighting* below).

**Sequence-only mode runs, but it is not the configuration used for the reported
benchmarks.** Without coverage there is no abundance signal, so contig
recruitment and abundance-guided bin splitting are both skipped and edges are
weighted by PMH composition similarity alone. Use it only for assemblies with no
reads available; expect materially lower bin quality on multi-sample datasets.

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
rabbitbin amber \
  --gold gsa_mapping.binning \           # CAMI bioboxes, needs _LENGTH
  --members results/out.members.tsv \    # or --binning preds.binning (bioboxes / 2-col)
  --output metrics_per_bin.tsv \
  --threads 64
```

## Default pipeline

Every stage below runs by default. Stages marked *(needs coverage)* are skipped
in sequence-only mode.

**Contig partitioning**

Contigs are split by length into three groups:

| Group | Length | Role |
|-------|--------|------|
| large | ≥ `--min-contig` (default 2500) | sketched, graphed, clustered |
| small | ≥ `--min-small-contig` and < `--min-contig` (default 1000–2499) | held back, recruited into finished bins |
| discarded | < `--min-small-contig` | dropped at parse time |

Both bounds are user-settable: `--min-contig` accepts any value ≥ 1500 and
`--min-small-contig` any value ≥ 500, so the recruitment window is
`[--min-small-contig, --min-contig)` and is 1000–2500 bp at the defaults.

**Graph construction**

1. **4-mer PMH sketching.** Each large contig is represented by a weighted
   ProbMinHash sketch over canonical 4-mer counts, with `--sketch-m`
   (default 500) registers.
2. **Mutual top-N candidate graph.** PMH similarity retrieves at most
   `--max-edges` (default 200) neighbours per contig; an undirected candidate
   edge is retained only when the neighbour relation is mutual.
3. **Edge weighting.** Surviving edges get
   `w = α · s_comp + (1 − α) · d`, and edges below `--min-edge-score`
   (default 0.70) are dropped. `s_comp` is the PMH composition similarity.
   `d` *(needs coverage)* is the conjunctive `min` of the Spearman rank
   correlation and the weighted Jaccard similarity of the two coverage
   profiles, with a hard cut on strongly negative correlation. `α` is chosen
   from the number of coverage samples `S`: `α = 1` for `S ≤ 2`, because a
   correlation over ≤ 2 paired observations has ≤ 0 degrees of freedom and
   carries no information, and `α = 0` for `S ≥ 3`, where PMH acts purely as a
   candidate filter. Override with `RABBIT_W_COMP`.

**Clustering**

4. **Fisher label propagation.** Each contig moves to the label whose incident
   edges carry the most aggregated support, combined by a Fisher-style
   nonlinear transform of the edge weights. Contigs are visited in a fixed
   order. A contig retains its current label when it is tied for the highest
   score; a contig that revisits an earlier label is frozen to break strict-score
   cycles. Propagation stops after a complete round with no label changes.

**Post-processing**

5. **Contig recruitment** *(needs coverage)*. Unbinned large contigs, then
   small contigs, are tested against the initial bins with a length-aware
   abundance model and assigned only when exactly one bin passes. Disable with
   `--no-recruit`.
6. **Singleton rescue.** Large contigs left unbinned are promoted to their own
   single-contig bins, subject to the output size filter.
7. **Abundance-guided bin splitting** *(needs coverage)*. Bins that are
   multi-modal in per-sample log-abundance are re-split by k-means, with `k`
   chosen by mean silhouette over `k = 2 … --split-max-k` (default 6) and the
   split accepted only when the best silhouette ≥ `--split-silhouette`
   (default 0.70). Disable with `--no-split`. Supplying `--marker-seed`
   replaces this with marker-guided splitting.
8. **Consolidate / secondary recruit** *(needs coverage)*. Same-genome
   fragments are merged and unbinned tails are compared with frozen bin cores
   using co-abundance shape and TNF. Confidence is the ratio of runner-up to
   winning cosine residuals, so absolute fit and separation form one statistic.
   Its boundary is learned per run from actual correct and incorrect leave-one-out
   core predictions (ROC/Youden); no fixed cosine or best-minus-second cutoff is
   used. Paired-end support is reported as corroborating evidence but cannot
   redirect the feature winner. Disable these stages with `RABBIT_BIN_MERGE=0`
   and `RABBIT_BIN_RECRUIT=0`.
9. **Output size filter.** Bins smaller than `--min-bin-size` (default
   200 000 bp) are not emitted.

The default workflow performs no marker-free subtraction/decontamination pass.
Heterogeneous bins are handled by splitting; marker-backed purification remains
available explicitly through `--markers ... --purify`.

### Experimental split-first refinement

`--simple-refinement` replaces post-processing stages 5–8 with a smaller,
order-independent workflow:

1. Split the initial large-contig bins using `log1p(coverage)`, seeded k-means
   with multiple initializations, and the best silhouette over `K=2…6`.
2. Reject a proposed partition in full if any child is smaller than
   `--split-min-sub-contigs` or `--split-min-sub-bp`.
3. Freeze the accepted, output-eligible bin cores that have at least
   `--min-recruit-cluster` valid coverage profiles. Compute each core's mean
   pairwise Spearman correlation once.
4. Assign unbinned long contigs in one batch, then short contigs in a second
   batch. Both batches compare only with the frozen cores, and a contig is
   assigned only when exactly one core's threshold is met.

This mode intentionally skips the default consolidate and secondary composition
recruitment passes. It is available for controlled ablation; it is not the
default quality preset.

Additional single-factor ablation controls are `--fixed-core-recruitment`
(fixed cores followed by the established refinement tail),
`--no-singleton-rescue`, `--stable-split-kmeans`, and
`--split-reject-small-children`. They are off by default.

Off by default, all requiring an explicit flag: SCG quality annotation
(`--qc`), purification (`--purify`), HQ-only output (`--keep-hq-only`),
composition-based recruitment (`--recruit-cutoff`), parameter search
(`--auto`, `--autotune`), consensus (`--ensemble`), and multi-resolution
output (`--resolutions`).

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
| `--min-edge-score` | 70 | Minimum edge weight, percent (1–99) |
| `--max-edges` | 200 | Maximum PMH neighbours per contig before mutual filtering |
| `--sketch-m` | 500 | Number of ProbMinHash registers |
| `--no-recruit` | off | Disable leftover/short-contig recruitment |
| `--simple-refinement` | off | Use experimental split-first, frozen-core batch recruitment |
| `--fixed-core-recruitment` | off | Test fixed-core recruitment while retaining the established refinement tail |
| `--no-singleton-rescue` | off | Disable promotion of output-sized unassigned long contigs |
| `--no-split` | off | Disable abundance-guided bin splitting |
| `--split-silhouette` | 0.70 | Minimum mean silhouette to accept a split |
| `--split-max-k` | 6 | Maximum sub-clusters per split bin |
| `--split-min-sub-contigs` | 3 | In simple refinement, reject a split with a smaller child |
| `--split-min-sub-bp` | 0 (off) | In simple refinement, reject a split with a shorter child |
| `--split-kmeans-restarts` | 10 | K-means initializations tested for each K |
| `--stable-split-kmeans` | off | Seed splitting from canonical bin membership rather than transient bin order |
| `--split-reject-small-children` | off | Apply the minimum-child guards to the established splitter |
| `--percent-identity` | 97 | Min read identity when reading BAMs |
| `--markers` | — | Contig→marker map, required by `--qc`/`--purify`/`--auto`/`--autotune` |
| `--qc` | off | Annotate `bins.tsv` with SCG completeness/contamination + MIMAG tier |
| `--bioboxes` | off | Also write a CAMI bioboxes `<prefix>.binning` |
| `--bin-fasta` | off | Also write per-bin FASTA files |
| `--unbinned` | off | Write unbinned contigs to FASTA |
| `--save-cache` / `--load-cache` | — | Cache the graph for fast re-binning |

Run `rabbitbin <command> --help` for the full option list, including the `qc`
and `refine` subcommands.

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

Everything else in the *Default pipeline* section above is deterministic given
the input files and thread count. No other flags were used for the published
benchmarks.

## Pipeline wrapper

`run_rabbitbin.sh` runs BAM depth summarization then RabbitBin in one call:

```bash
run_rabbitbin.sh assembly.fa sample1.bam sample2.bam
```

## License

RabbitBin is released under the LBNL BSD license. Portions of the graph-clustering
pipeline derive from earlier open-source metagenome binning work; see `license.txt`.
