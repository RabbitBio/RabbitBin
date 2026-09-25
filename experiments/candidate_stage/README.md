# CAMI2 candidate-stage and clustering audit

Base: `RabbitBio/RabbitBin` `origin/main` at
`bf9e1013941bf7864501e1ac771c77db73525d6e` (verified remotely on
2026-09-25). This was an isolated, uncommitted worktree at stage completion;
it was subsequently archived to an [opt-in research branch](../README.md) at
the user's request. The production default path is unchanged.
Inputs are the same existing CAMI2 gold-standard assemblies
and remapped BAM lists as the MetaCAT comparison. Contigs below 1,000 bp are
excluded; RabbitBin builds its graph on contigs >=2,500 bp and recruits shorter
ones later. RabbitBin `amber` evaluates only bins >=200,000 bp. Seed 42,
64 threads, PMH sketch `m=500`, mutual top-200, dual-depth 5 throughout.

## 1. Which stage loses correct long-contig links?

`--audit-candidate-stages` loads gold labels solely for diagnostics. For each
labelled same-genome pair it evaluates the *actual final* coverage gate, exact
PMH score, both saved directed top-200 heap cutoffs, and final retained edge.
The audit also computes connected components of true retained edges. The
`mutual_top_same` and `retained_same` counts match for every eligible node in
all three datasets; RabbitBin membership files are byte-identical to the
previous non-audit baseline. Gold never changes candidate construction,
clustering, splitting, recruitment, or bins.

The following denominator is gold-labelled contigs from genomes with at least
two >=2,500 bp contigs that are **outside that genome's largest RabbitBin bin**.
This is a diagnostic fragmentation measure, not AMBER HQ loss. Categories are
exclusive and ordered by earliest definitive missing stage. In particular,
"PMH top" means a coverage-feasible, positive-PMH same-genome peer exists but
none survives the contig's own or reciprocal top-200 heap.

| Dataset | Missing large bp | Coverage-infeasible | PMH own/reciprocal top | True-edge connected to best bin but assigned elsewhere | True-edge retained, disconnected from best bin |
| --- | ---: | ---: | ---: | ---: | ---: |
| Marine | 143.76 Mb | 8.65% | 41.61% | 30.78% | 8.08% |
| Plant-associated | 351.38 Mb | 56.18% | 1.07% | 15.64% | 24.67% |
| Strain-madness | 897.69 Mb | 5.60% | 63.34% | 10.39% | 19.81% |

Other categories (genome without a RabbitBin bin, zero positive PMH) account
for the remainder. In Strain-madness, contigs with no correct neighbour in
their own top-200 have a median **200 known different-genome contigs** in that
top-200; reciprocal failures have median 198. This directly verifies
top-200 displacement **at the candidate-edge level**. It does not establish
that enlarging top-200 would improve HQ: a previous N-sweep and a joint
PMH×coverage ranking both added many incorrect edges and degraded purity.

The default final coverage weight is `min(max(Spearman,0), Jcov)`, with cutoff
0.715331862959. Among missing large bp with no coverage-feasible same-genome
peer, the subgate decomposition is:

| Dataset | Jcov-only failure | Both fail | Spearman-only failure | Each passes, but on different peers |
| --- | ---: | ---: | ---: | ---: |
| Marine | 91.64% | 6.48% | 0.31% | 1.58% |
| Plant-associated | 39.10% | 36.55% | 13.09% | 11.26% |
| Strain-madness | 99.68% | 0.31% | 0% | 0.01% |

These are contig-level exclusive categories within the coverage-infeasible
subset, weighted by contig bp. The final two-stage coverage score, not just
Spearman, is needed to diagnose the Plant failure.

## 2. Focus on the genomes behind MetaCAT's HQ advantage

The same local gold file and bin-size/purity/completeness rules reproduce the
`amber` HQ bin counts exactly. MetaCAT-only HQ genomes number 68/10/16 for
Marine/Plant/Strain; RabbitBin-only HQ genomes number 4/5/7, giving net
differences +64/+5/+9. Consider the true contigs *inside those MetaCAT HQ
bins* but outside the corresponding genome's largest RabbitBin bin:

| Cause / fraction of those missing true bp | Marine (35.11 Mb) | Plant (10.31 Mb) | Strain (36.00 Mb) |
| --- | ---: | ---: | ---: |
| 1,000–2,499 bp, excluded from RabbitBin's initial graph | **55.5%** | 3.4% | 0.2% |
| True-edge connected to best RabbitBin bin but assigned elsewhere | 29.3% | **72.4%** | **44.0%** |
| True-edge retained, disconnected from best bin | 1.5% | 0.7% | 22.7% |
| PMH own/reciprocal top-200 excluded | 11.7% | 12.5% | **31.3%** |
| No coverage-feasible same-genome peer | 1.8% | 10.7% | 1.6% |

Marine's short-contig category is 19.50 Mb, of which 19.21 Mb is unbinned by
RabbitBin. In Plant, 7.23 of 7.47 Mb in the connected-but-elsewhere category
is in another RabbitBin bin; in Strain the corresponding number is 15.51 of
15.86 Mb. Thus global missing-edge statistics cannot stand in for the HQ gap.
These fractions describe already recovered true bp, **not the number of HQ
bins each intervention would gain**.

On the fixed retained graph, 93%/92%/85% of the bp assigned to another bin in
the connected-but-elsewhere category have at least as much *direct total edge
weight* to their current bin as to the best RabbitBin bin (Marine/Plant/Strain).
Simple local reassignment by graph-edge sum is therefore unlikely to solve
the fragmentation; this is an evidence check, not a normalized-score proof.

## 3. Same-graph clustering and simple ablations

`--export-retained-graph` writes the post-coverage weighted graph (the recorded
runs export before LPA's tiny `<1` numerical clamp);
`--external-labels` substitutes a graph-clustering membership at the LPA
step. Every variant then goes through RabbitBin's unchanged singleton rescue,
abundance splitting, and short-contig recruitment. Infomap and Leiden use
`igraph 1.0.0`; Infomap uses its default 10 trials and Leiden uses weighted
modularity to convergence. No resolution sweep or gold-driven selection.

| Dataset | Method | HQ / MQ | Weighted purity | Extra clustering time |
| --- | --- | ---: | ---: | ---: |
| Marine | RabbitBin LPA | 286 / 353 | 0.9264 | — |
| Marine | Infomap | 292 / 355 | 0.9233 | 46.94 s |
| Marine | Leiden | 271 / 333 | 0.8586 | 5.41 s |
| Plant | RabbitBin LPA | 84 / 94 | 0.9797 | — |
| Plant | Infomap | 78 / 94 | 0.9452 | 10.64 s |
| Plant | Leiden | 73 / 91 | 0.9299 | 1.21 s |
| Strain | RabbitBin LPA | 33 / 38 | 0.4038 | — |
| Strain | Infomap | 36 / 43 | 0.3706 | 4.94 s |
| Strain | Leiden | 37 / 44 | 0.3315 | 0.41 s |

Infomap is not a robust default improvement (Plant regression, all three
purities lower, Marine clustering alone slower than the entire ~14 s baseline).
Leiden degrades purity sharply. On the identical cached candidate graph,
changing the coverage edge score to correlation-only lowers HQ from
286/84/33 to 180/64/30 and weighted purity from 0.9264/0.9797/0.4038 to
0.8009/0.8930/0.1582. Disabling abundance splitting lowers HQ to
283/83/28. Marine with recruitment disabled yields 288 HQ and 0.9333
weighted purity, but loses MQ and completeness; the present recruiter adds
both correct and incorrect short contigs.

## 4. Length-aware self-supervised short-contig rescue

An optional `--seq-recruit` prototype borrows the *principle* of fragment
views from SemiBin2/COMEBin without adding a neural-network dependency or
dataset-specific tuning. It learns the sequence-confidence boundary from
held-out central fragments of existing long contigs, with leave-parent-out
bin profiles and the **existing** recruitment FPR limit (0.05). Queries use
canonical 4-mer counts; pseudo-fragment length is the observed median length
of 1,000–2,499 bp contigs. Gold is never used for calibration or binning.
The stricter `--seq-recruit-agree` additionally requires that the sequence and
existing coverage recruiter's best core be identical; no new numerical weight
or threshold is introduced.

| Dataset | Baseline HQ/MQ, purity | Sequence-only HQ/MQ, purity | Agreement HQ/MQ, purity | Agreement recruited shorts |
| --- | ---: | ---: | ---: | ---: |
| Marine | 286/353, 0.9264 | 271/343, 0.9177 | 289/356, 0.9269 | 9,273 |
| Plant-associated | 84/94, 0.9797 | 77/93, 0.9449 | 84/95, 0.9795 | 2,193 |
| Strain-madness | 33/38, 0.4038 | 33/38, 0.4032 | 33/38, 0.4038 | 71 |

Pure sequence rescue harms Marine and Plant. Evidence agreement gives a small
Marine gain but no HQ gain in the other two datasets. The original agreement
run took 25.05/19.41/74.47 s (Marine/Plant/Strain) and peaked at
3.25/4.01/1.83 GB; baseline clustering runs took 14.01/24.84/66.67 s and
peaked at 1.49/2.74/0.57 GB. Wall-time comparisons are sensitive to BAM and
FASTA cache warmth, especially Marine/Plant. The prototype remains opt-in:
three-dataset quality is not sufficiently improved, and the memory overhead
is material. `run_seq_recruit.sh` reproduces it.

An extra feature-cache experiment kept only the 4-mer profiles and released
the FASTA bytes before graph construction. It reproduced byte-identical bin
memberships on all three datasets. Its observed peaks were 3.58/2.91/1.75
GB: an improvement on Plant/Strain but a regression on Marine because the
original bytes and profile arrays overlap during feature construction. It
therefore does not remove the cost objection.

## 5. What to borrow, and what the experiments rule out

- [SemiBin2](https://pmc.ncbi.nlm.nih.gov/articles/PMC10311329/) creates
  self-supervised must-links by cutting long contigs and trains a contrastive
  embedding; its randomly sampled cannot-links are a useful training signal
  but can accidentally pair same-genome contigs. The transferable idea
  is *length-matched, held-out same-contig evidence*, not the assumption that
  every random contig pair differs. Our cheap raw-4-mer version was not robust
  enough to replace the learned representation.
- [COMEBin](https://www.nature.com/articles/s41467-023-44290-z) creates
  multiple fragment views and learns a joint coverage/composition embedding,
  then clusters with Leiden and marker-aware selection. Our fixed-graph Leiden
  ablation establishes that its gain cannot be attributed to Leiden alone;
  the representation and graph construction are the important coupled parts.
- [MetaCAT](https://www.nature.com/articles/s41564-026-02472-7) uses single-copy
  gene seeds, multiple coverage/composition affinity models, marker-based
  selection, and either semi-supervised label propagation or a sparse weighted
  mixture model. Its coverage model also uses variance; RabbitBin's BAM route
  currently retains means only. The prior direct marker transplant into
  RabbitBin's late split did not improve HQ, and marker generation alone cost
  86–277 s across these datasets (see sibling `metacat_borrow` report).

**Priority under RabbitBin's runtime constraint:** (1) improve the *evidence*
used to assign or combine surviving same-genome fragments, with a
self-supervised confidence estimate and explicit negative controls; (2)
explore an inexpensive, length-aware second candidate channel in Strain,
but require measured false-positive control before admitting extra edges;
(3) treat short-contig recruitment as a selective Marine opportunity, not a
global fix. Avoid defaulting a larger top-N, raw 4-mer recruitment, a different
community detector, or late marker splitting: the completed three-dataset
ablations already show their quality/performance problems. Any future model
must be selected without gold labels or per-dataset parameter tuning and must
beat all three frozen baseline AMBER reports in both quality and resource cost.

## Implications

1. **Marine:** prioritize an uncertainty-calibrated 1–2.5 kb recruiter that
   combines sequence and coverage evidence. Do not merely recruit more by
   weakening the existing FPR gate; its current net purity cost is measurable.
2. **Plant and Strain:** investigate why long contigs with surviving true paths
   still form other bins. Existing graph edges often support their current bin
   more strongly, and both Leiden/Infomap and no-split ablations fail as
   universal fixes. Better evidence or a calibrated representation is needed,
   not just a different community-detection algorithm.
3. **Strain:** PMH top-200 crowding is now directly verified and relevant to
   MetaCAT-only HQ genomes, but globally widening candidates or dropping Jcov
   loses purity. A new candidate channel needs independent discrimination
   (for example reliable marker constraints when precomputed, or a
   self-supervised length-aware score) and must be evaluated before defaulting.

No dataset-specific thresholds or gold-derived parameters were added. The
new command-line options are strictly opt-in diagnostics/experimental hooks.

Reproduce with `run.sh`, `summarize.py`, `compare_hq.py`, `edge_support.py`,
`run_clustering.sh`, `run_corr_ablation.sh`, `run_no_split.sh`, and
`run_seq_recruit.sh` in this
directory. Full generated results are in the ignored `results_subgates/` and
`cluster_results/` and `seq_results/` directories.
