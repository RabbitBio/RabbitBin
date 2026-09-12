# RabbitBin

Fast, sketch-based metagenome binning. The default RabbitBin pipeline uses
canonical, count-weighted 4-mer ProbMinHash (PMH) sketches to construct a
bounded mutual-nearest-neighbour candidate graph, uses abundance profiles as
the edge evidence in the standard multi-sample setting, clusters the retained
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
is inferred from the flags and echoed in the log line beginning `Edge weighting:`.

| Mode | Flags | Edge weight | Abundance stages |
|------|-------|-------------|------------------|
| BAM/CRAM | `--fasta` + `--bam`/`--bam-list` | coverage only for `S >= 3`; composition for `S <= 2` (depth computed in-process) | enabled |
| Precomputed depth | `--assembly` + `--depth` | coverage only for `S >= 3`; composition for `S <= 2` | enabled |
| Sequence-only | `--assembly` alone | composition only | **disabled** |

The manuscript's multi-sample method uses PMH to select candidate neighbours
and coverage alone to weight their edges when `S >= 3`.

For one or two coverage samples, the original composition-weighted fallback is
retained; coverage still drives recruitment and splitting. In the two-sample
path, the existing negative-correlation gate and minimum-edge cutoff still
filter candidates before composition weights are used. Spearman is undefined with one
observation and, without ties, is restricted to +1 or -1 with two. Three
observations are the first sample count with non-binary rank-correlation
resolution. This fallback is separate from the manuscript's `S >= 3` method.

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
2. **Bounded mutual candidate graph.** In the standard multi-sample path, an
   exact abundance-feasibility bound first removes pairs that cannot pass the
   final edge threshold. PMH similarity then retains at most `--max-edges`
   (default 200) neighbours per contig from the feasible pairs. The production
   candidate graph keeps a pair only when the neighbour relation is mutual.
   This order prevents abundance-incompatible high-PMH pairs from consuming
   the bounded neighbourhood.
3. **Coverage-only edge weighting.** In the default method with `S >= 3`,
   each candidate edge receives

   $$w_{ij}=\min\{\max(\rho_{ij},0),J^{\mathrm{cov}}_{ij}\}.$$

   Here `rho` is the Spearman correlation of the contigs' coverage profiles
   across samples. `Jcov = sum_s min(x_is, x_js) / sum_s max(x_is, x_js)` is
   weighted Jaccard on per-sample mean depths. By default, `x_is` is divided
   by the mean depth of sample `s` across the retained large contigs, so library
   size differences do not dominate this magnitude term. Constant rank profiles
   have zero correlation support; an all-zero Jaccard denominator gives zero
   support. PMH composition similarity does not enter `w`.

   Edges with `w < 0.70` are dropped by default (`--min-edge-score 70`, expressed
   as a percentage). This is a survival threshold on the coverage weight, not a
   mixing coefficient. No edge-power transform is applied by default. The former
   `RABBIT_W_COMP` override is ignored, with a message when it is set.

   Optional coverage-metric ablations (`RABBIT_DEPTH_SIM`, `RABBIT_DEPTH_FUSE`),
   edge-power/SNN transforms, and auxiliary GFA/SNV evidence are separate from
   this default formula. Parameter search and certification also use coverage
   weights for `S >= 3`; they never reintroduce composition mixing.

### PMH representation validation

`--validate-pmh-gold` is a terminal, gold-aware diagnostic for testing the PMH
candidate representation independently of abundance and clustering. It samples
large contigs under `--seed`, computes their exact sequence-only PMH top-N lists
from the same packed winner representation used by graph construction, writes a
TSV, and exits before constructing the production graph. Gold labels never enter
the normal binning path.

```bash
rabbitbin bin \
  --assembly contigs.fa \
  --output validation/run \
  --sketch-m 500 \
  --seed 42 \
  --validate-pmh-gold gold.binning \
  --validate-pmh-queries 1000 \
  --validate-pmh-top 400
```

The report includes candidate precision, recall against all recoverable
same-genome contigs, the fractions of queries with at least 1, 5, and 10
same-genome neighbours, and MRR at N = 50, 100, 200, and 400. This mode is for
controlled evaluation only and requires CAMI/bioboxes-style gold assignments.

`--audit-graph-gold` is the complementary production-path diagnostic. It loads
gold labels only after candidate generation and abundance edge scoring, reports
true/false candidate and retained edges plus per-contig true-neighbour coverage,
and then lets the unchanged binning path continue. The pure-PMH and production
reports must not be conflated: the latter also reflects the abundance-feasibility
gate and mutual-neighbour requirement.

**Clustering**

4. **Fisher label propagation.** Each contig moves to the label whose incident
   edges carry the most aggregated support, combined by a Fisher-style
   nonlinear transform of the edge weights. Contigs are visited in a fixed
   order. A contig retains its current label when it is tied for the highest
   score; a contig that revisits an earlier label is frozen to break strict-score
   cycles. Propagation stops after a complete round with no label changes.

**Post-processing**

The default refinement path is: singleton output safeguard → abundance-guided
splitting → one selective coverage recruitment → output-size filtering.

5. **Singleton output safeguard.** Before splitting, an unassigned large contig
   that is itself at least `--min-bin-size` long is retained as a single-contig
   bin. This is an output safeguard, not a recruitment pass.
6. **Abundance-guided bin splitting** *(needs coverage)*. Bins that are
   multi-modal in per-sample log-abundance are re-split by k-means, with `k`
   chosen by mean silhouette over `k = 2 … --split-max-k` (default 6) and the
   split accepted only when the best silhouette ≥ `--split-silhouette`
   (default 0.70). Disable with `--no-split`. Supplying `--marker-seed`
   replaces this with marker-guided splitting.
   Silhouette is averaged over contigs, using a random subset of up to 600
   contigs for large bins. Singleton observations contribute zero and remain
   in the average. For threshold sweeps, `--resolutions` reuses the full
   in-memory state and starts each refinement from the same pre-split bins;
   `RB_SPLIT_AUDIT=1` writes per-parent decisions to `.split_audit.tsv`.
7. **Selective post-split recruitment** *(needs coverage)*. The split bins at
   least `--min-bin-size` long are frozen as cores. All remaining unbinned long
   and short contigs are compared with every core using their coverage
   trajectories. With at least three samples, RabbitBin rank-transforms and
   normalizes each coverage profile, represents each core by the normalized sum
   of its member profiles, and uses cosine similarity to score contig-core
   matches (equivalent to correlation on ranks for individual profiles). With
   one or two samples, for which rank correlation is undefined or nearly binary,
   it instead uses weighted Jaccard between the contig's normalized raw coverage
   and the core's mean coverage, preserving abundance magnitude. If the best and
   second-best core scores are `s_best` and `s_second`, the recruitment
   confidence is

   $$C=\log\frac{1-s_{\mathrm{second}}}{1-s_{\mathrm{best}}}.$$

   RabbitBin learns one acceptance boundary per run without reference labels.
   Each existing core member is temporarily left out of its own centroid and
   classified against the frozen cores. Correct returns to the source core and
   incorrect returns to another core supply the positive and negative confidence
   distributions. Only the best other core is required, so calibration also
   works with exactly two cores. If one outcome class is absent, RabbitBin uses
   the paired source-core and strongest-wrong-core counterfactual scores from
   the same leave-one-out members to supply that class instead of applying a
   fixed cutoff. The resulting values form an internal ROC curve, and RabbitBin
   selects the boundary that maximizes Youden's
   `J = TPR - FPR`. An unassigned contig is recruited into its best-scoring core
   only when its confidence reaches that boundary. The same learned boundary
   and scoring rule are used for long and short contigs. This pass
   never moves an already binned contig and does not merge bins. Disable it with
   `--no-recruit` or set `RABBIT_BIN_RECRUIT=0` for an environment-controlled
   ablation.
8. **Output size filter.** Bins smaller than `--min-bin-size` (default
   200 000 bp) are not emitted.

The default workflow performs no marker-free subtraction/decontamination pass.
Heterogeneous bins are handled by splitting; marker-backed purification remains
available explicitly through `--markers ... --purify`.

Off by default, all requiring an explicit flag: SCG quality annotation
(`--qc`), purification (`--purify`), HQ-only output (`--keep-hq-only`),
parameter search (`--auto`, `--autotune`), consensus (`--ensemble`), and
multi-resolution output (`--resolutions`).

The optional graph-reuse search sweeps edge powers (`--no_gold`,
`--auto`, `--ensemble`, `--autotune`); `--autotune` also searches split
silhouettes. Its default powers are 1, 1.25, 1.5, 2 and 3. A custom
`RABBIT_REUSE_SWEEP="1.0;1.5;2.0"` lists powers only; the former `alpha:power`
syntax is rejected. Label-free modularity is measured on the fixed baseline
coverage graph for `S >= 3` (the composition fallback for `S <= 2`). Rebuild
candidate caches when changing graph-construction settings: cached topology is
reused, while edge weights are recomputed.

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
| `--min-edge-score` | 70 | Minimum edge weight, percent (2–99); coverage weight for `S >= 3` |
| `--max-edges` | 200 | Maximum PMH neighbours per contig among production-feasible pairs, before mutual filtering |
| `--sketch-m` | 500 | Number of ProbMinHash registers |
| `--validate-pmh-gold` | — | Evaluate sequence-only PMH top-N neighbourhoods against CAMI gold, write TSV, and exit |
| `--validate-pmh-queries` | 1000 | Seed-controlled labelled queries used by PMH validation |
| `--validate-pmh-top` | 400 | Largest neighbourhood retained by PMH validation |
| `--audit-graph-gold` | — | Audit production candidate and abundance-retained edges against CAMI gold without changing binning |
| `--audit-graph-out` | `<output>.graph_audit.tsv` | Output TSV for `--audit-graph-gold` |
| `--no-recruit` | off | Disable the post-split long/short-contig coverage recruitment |
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
