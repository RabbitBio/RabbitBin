# Additional commands and options

For installation and the standard BAM/depth workflow, start with the
[README](../README.md). These details are optional; ordinary binning needs
neither marker tools nor reference genome labels.

## Commands

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

## Coverage modes

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
[pipeline documentation](pipeline.md). Coverage also drives splitting
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

## Generate a depth table

```bash
rabbitbin depth --bam-list bams.txt --out depth.tsv --threads 8
```

## Evaluate against a gold standard (AMBER-compatible)

```bash
# CAMI bioboxes gold mapping; the gold file must provide _LENGTH.
# Use --binning preds.binning instead of --members for a two-column prediction.
rabbitbin amber \
  --gold gsa_mapping.binning \
  --members results/out.members.tsv \
  --output metrics_per_bin.tsv \
  --threads 64
```

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

To generate a marker map, install Prodigal, HMMER's `hmmsearch`, and a
single-copy marker HMM collection. Set `RABBITBIN_MARKER_HMM` to that HMM file,
then run `rabbitbin_markers.sh contigs.fa contigs.markers.tsv 32`.
These external tools are not needed if you already have a compatible marker
map, and neither tools nor markers are required for default binning.

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

The default wrapper does not require external `samtools` or Python. For new
analyses, direct `rabbitbin bin --bam ...` is the simpler entry point.
