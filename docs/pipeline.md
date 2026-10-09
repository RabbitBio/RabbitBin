# Default pipeline and PMH validation

[Back to the README](../README.md)

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
   ProbMinHash sketch over canonical 4-mers, with `--sketch-m`
   (default 500) registers. The default weight measures enrichment relative
   to the contig's own base composition:

   $$w_c(x)=\frac{f_c(x)/\sum_y f_c(y)}{e_c(x)},
   \qquad
   e_c(x)=P_c(x)+\mathbf{1}_{x\neq\operatorname{rc}(x)}P_c(\operatorname{rc}(x)),
   \qquad
   P_c(x)=\prod_{\ell=1}^{4}\pi_c(x_\ell).$$

   Here `f_c(x)` counts canonical 4-mer occurrences, `pi_c(b)` is the frequency
   of base `b` among the contig's A/C/G/T bases, and `rc(x)` is the reverse
   complement of `x`. Palindromic 4-mers (`x = rc(x)`) are counted once and
   therefore use `e_c(x)=P_c(x)` rather than `2P_c(x)`. A denominator at or
   below `1e-14` gives zero weight.
2. **Bounded mutual candidate graph.** In the standard multi-sample path, an
   exact abundance-feasibility bound first removes pairs that cannot pass the
   final edge threshold. PMH similarity then retains at most `--max-edges`
   (default 200) neighbours per contig from the feasible pairs. The production
   candidate graph keeps a pair only when the neighbour relation is mutual.
   This order prevents abundance-incompatible high-PMH pairs from consuming
   the bounded neighbourhood.
3. **Coverage-only edge weighting.** With one or two input samples, each
   candidate edge receives the mean ratio over informative coverage dimensions:

   $$w_{ij}=A_{ij}=\frac{1}{|\mathcal I_{ij}|}\sum_{\ell\in\mathcal I_{ij}}
   \frac{\min(x_{i\ell},x_{j\ell})}
        {\max(x_{i\ell},x_{j\ell})},\qquad 1\le S\le2.$$

   Here `x_i` contains the available coverage features, including ordinary and
   MAPQ-filtered depth for default BAM inputs. `I_ij` contains dimensions with
   nonzero coverage in at least one contig; its size is the number of valid
   coverage dimensions, not the sample count `S`. A shared zero is omitted; a
   zero on only one side contributes zero. An empty `I_ij` gives zero support.
   Per-feature normalization cancels within each ratio, so every informative
   dimension contributes equally. The same magnitude-based definition applies
   to single- and two-sample inputs, without correlation calculation or gating.

   With `S >= 3`, each candidate edge receives

   $$w_{ij}=\min\{\max(\rho_{ij},0),J^{\mathrm{cov}}_{ij}\}.$$

   Here `rho` is the Spearman correlation of the contigs' coverage feature
   vectors. `Jcov = sum_l min(x_il, x_jl) / sum_l max(x_il, x_jl)` is weighted
   Jaccard over the same dimensions. By default, each feature is divided by its
   mean across the retained large contigs, so library size differences do not
   dominate this magnitude term. Constant rank profiles
   have zero correlation support; an all-zero Jaccard denominator gives zero
   support.

   Edges with `w < tau` are dropped, with the default
   `tau = 0.7153318629591614` (`--min-edge-score 71.53318629591614`, expressed
   as a percentage). The cutoff is the neutral point at which two equal-weight
   Fisher supports have the same score as one support. It is the unique solution
   in `(0, 1)` of

   $$(1-\tau)[1-2\ln(1-\tau)]=1.$$

   This adopts the two-term Fisher soft truncation boundary on the weight scale
   ([Zaykin et al., 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC2569904/))
   as a local scoring consistency rule. No edge-power transform is applied by
   default. The former
   `RABBIT_W_COMP` override is ignored, with a message when it is set.

   Optional coverage-metric ablations (`RABBIT_DEPTH_SIM`, `RABBIT_DEPTH_FUSE`),
   edge-power/SNN transforms, and auxiliary GFA/SNV evidence are separate from
   this default formula. Parameter search and certification also use coverage
   weights whenever coverage is available.

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
   it uses the same mean coverage-feature min/max ratio as graph weighting,
   comparing the contig's coverage with the arithmetic mean of the core's
   nonzero long-contig profiles. Default BAM profiles include ordinary and
   MAPQ-filtered coverage. Jointly zero dimensions are omitted as above. The
   core mean excludes the current member during leave-one-out calibration.
   If the best and second-best core scores are `s_best` and `s_second`, the recruitment
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
   selects the boundary that maximizes Youden's `J = TPR - FPR` among operating
   points with `FPR <= 0.05`. Tied optima use the largest (most conservative)
   threshold. An unassigned contig is recruited into its best-scoring core
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
coverage graph. Rebuild candidate caches when changing graph-construction
settings: cached topology is reused, while edge weights are recomputed.
Cache format v4 stores the input sample count separately from coverage columns;
older caches must be rebuilt to preserve the low-sample graph semantics.
