# Opt-in RabbitBin research experiments

Snapshot branch: `experiments/joint-representation-20260925`.
Base: `bf9e1013941bf7864501e1ac771c77db73525d6e` (origin/main,
verified again on 2026-09-25). This branch is not a new production default.

## Status: paused at the user's request

The base self-supervised joint/block/coverage CAMI2 round is complete.
The signal_joint/signal_block/signal_coverage variants have code and synthetic
algebra checks only; their CAMI2 experiments have **not** been run.
No bin-level merging experiment has started. Do not launch experiments just
because this branch has been checked out; wait for explicit authorization.

The native experimental paths require explicit CLI flags or environment
variables. Without those switches, the original binning defaults are retained:

| Entry point | Explicit opt-in |
| --- | --- |
| Same-graph alternative clustering | `--external-labels` |
| Learned/external graph | `--external-graph`, optionally `--external-graph-coverage` |
| Actual half-contig BAM depth export | `--export-fragment-depth` |
| Sequence-based short recruitment | `--seq-recruit` or `--seq-recruit-agree` |
| Full coverage gate before PMH top-N | `RABBIT_CANDIDATE_COVERAGE=1` |
| Stable Fisher log-tail / other aggregation | `RABBIT_LPA_SCORE=logtail` / explicit alternative |
| Read-only diagnostics | `--audit-candidate-stages`, `--export-retained-graph`, `RB_LPA_AUDIT` |

Some environment flags are checked by **presence**, not value: to disable them,
unset them rather than assigning `0`. No experimental variable is installed in
shell startup files by this repository. Python packages are research-script
dependencies only, not new dependencies for ordinary native RabbitBin binning.

## Reading order

1. [Current stage summary](self_joint/README.md)
2. [Current protocol and pending signal metric](self_joint/DESIGN.md)
3. [Within-architecture ablations](architecture_local/README.md)
4. [Earlier joint-space controls](joint_representation/README.md)
5. [Candidate-stage, MetaCAT gap and recruitment audit](candidate_stage/README.md)
6. [Migration and resumption instructions](RESUME.md)

The same CAMI2 Marine, Plant-associated and Strain-madness inputs are used.
Contigs >=1,000 bp, graph contigs >=2,500 bp, bins >=200,000 bp, RabbitBin AMBER,
seed 42, PMH m=500, mutual top-200, dual-depth 5. No dataset-specific algorithm
selection or tuning; gold is evaluation/diagnostic input, never training input.

Generated results, caches, graph files, BAMs, assemblies, private conversation
exports, credentials and virtual environments are not committed to Git.
