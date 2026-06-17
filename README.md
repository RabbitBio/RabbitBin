# RabbitBin

Sketch-based metagenome binning: fast composition–abundance graph clustering with optional abundance-guided bin splitting.

## Quick start

```bash
cd build && cmake .. && make rabbitbin -j

# Bin with a precomputed depth file
./build/src/rabbitbin \
  --assembly contigs.fa \
  --depth depth.tsv \
  --output results/out \
  --threads 32
```

## Main outputs (default naming)

| File | Description |
|------|-------------|
| `prefix.members.tsv` | Contig-to-bin membership |
| `prefix.bins.tsv` | Per-bin stats |
| `prefix_bin_001.fa` | Bin FASTA files |
| `prefix.unbinned.fa` | Unbinned contigs (with `--unbinned`) |

After pulling these changes, rebuild from a clean tree:

```bash
rm -rf build && mkdir build && cd build && cmake .. && make -j
```

## Pipeline wrapper

`run_rabbitbin.sh` runs depth summarization from BAMs, then RabbitBin:

```bash
run_rabbitbin.sh assembly.fa sample1.bam sample2.bam
```

## Key options

| Option | Meaning |
|--------|---------|
| `-a, --assembly` | Input FASTA |
| `-o, --output` | Output prefix |
| `-d, --depth` | Coverage depth TSV |
| `-t, --threads` | Threads (0 = all CPUs) |
| `--sim-cutoff` | Composition similarity cutoff (0 = auto) |
| `--sketch-k/m` | Sketch k-mer size / sketch size |
| `--no-split` | Disable abundance bin splitting |

## Environment variables

- `RABBIT_PMH` — weighted ProbMinHash graph (default on)
- `RABBIT_MUTUAL_KNN` — mutual k-NN graph (default on)
- `RABBIT_GC_NORM` — GC normalization mode
- `RABBIT_W_COMP` — composition vs abundance edge weight
- `RABBIT_NEG_DEPTH` — negative depth-correlation filter (default -0.3)
- `RABBIT_SPLIT_SIL` — silhouette split threshold

## Source layout

```text
src/
  rabbitbin.cpp          # CLI + main pipeline
  rabbit_depth.cpp       # BAM → depth TSV utility
  rabbit_overlap.cpp     # cross-assembly overlap utility
  rabbitbin.h            # shared types and state
  probmh.cpp/h           # ProbMinHash composition metric
  rabbit_sketch.h        # k-mer sketch embedding
  rabbit_invidx.h        # inverted index for sketch keys
  impl/
    rb_cluster.cpp       # label propagation
    rb_split.cpp         # abundance / marker bin splitting
    rb_graph.cpp         # similarity graph build + calibration
    rb_abundance.cpp     # abundance distance helpers
    rb_output.cpp        # bin output writers
```

## License

RabbitBin is released under the LBNL BSD license. Portions of the graph-clustering
pipeline derive from earlier open-source metagenome binning work; see `license.txt`.
