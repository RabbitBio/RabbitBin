# Joint representation experiments (specified before evaluation)

Base: RabbitBio/RabbitBin origin/main bf9e1013941bf7864501e1ac771c77db73525d6e.
Use the same BAM-derived cached coverage, assemblies, seed 42, 64 threads,
1,000 bp input floor, 2,500 bp graph floor, mutual top-200, 200,000 bp bin floor,
and RabbitBin AMBER as the prior experiments. Gold is evaluation-only.

The hypothesis is that switching similarity measures between candidate
selection and edge weighting loses a coherent geometry. Test actual joint
feature spaces; earlier multiplication of PMH and coverage scores is not a
test of a learned joint representation.

## Representations

Use canonical 4-mer frequencies (136 dimensions, square-root transform) and
the existing all-read plus MAPQ>=5 coverage dimensions (column library
normalization, log1p). Center each feature block and divide it by its global
root-mean-square row norm: both views have equal average energy, independently
of their dimension. No per-dataset weight or dimension sweep.

1. **concat**: Euclidean distance in the normalized concatenation.
2. **whiten**: OAS shrinkage covariance of the concatenation, then its inverse
   square-root transform. Cross-view covariance enters the full metric.
3. **cca**: OAS-regularized CCA between composition and coverage for the same
   contigs. Average the aligned canonical coordinates, weighting each by the
   square-root of its canonical correlation. Keep all numerically nonzero
   directions; no tuned latent dimension. This is a linear same-contig
   cross-view representation, not a neural contrastive model.

Exact Euclidean top-200 neighbours are searched in each representation;
mutual neighbours define edges. For each node, the squared distance to its
200th neighbour is the local scale. Edge affinity is
exp(-squared_distance / sqrt(scale_i * scale_j)). Thus one joint distance
controls both topology and weights, with no new scalar similarity cutoff.
The local scale uses RabbitBin's existing neighbour count; no labels or
dataset-specific constants enter it. Zero distances use a numerical floor
only to avoid division by zero.

## Factorial controls

For each representation compare:

- new candidate topology + original coverage weight/filter (candidate-only);
- original **retained** topology + joint affinity (weight-only, identical edges);
- new candidate topology + joint affinity (coherent joint graph).

Keep RabbitBin LPA, singleton rescue, abundance splitting, and recruitment
unchanged initially. This isolates the graph stages; a successful graph can
subsequently justify extending the same geometry to recruitment/refinement.
Report HQ/MQ, purity, graph true/false edges, same-genome neighbour coverage,
representation/search time and peak memory. Python/FAISS prototyping time is
not claimed as integrated RabbitBin end-to-end runtime. A promising variant
must subsequently be measured from BAM and pass a performance check before
becoming a default.

No CAMI labels are read by feature extraction, representation training,
neighbour search, or RabbitBin binning. Freeze all three methods before
looking at the three-dataset outcomes; report failed variants as well.

## Aggregation controls (added after the first three-dataset run)

The first run showed a strong Strain gain but Marine/Plant regressions.
Inspection identified a semantic mismatch: a Gaussian kernel affinity is fed
to RabbitBin's Fisher combination, whose existing 0.715331862959 cutoff is
specifically the two-support neutral point. Therefore add these controls,
again identically on all three datasets, without fitting to gold:

- Standard weighted-sum LPA on the baseline, each joint graph, and each
  weight-only graph (the latter retains exactly the original surviving edges).
- Joint affinities on the subset surviving the unchanged coverage gate,
  with both Fisher and sum LPA; this isolates loss of the original rejection
  constraint from the choice of edge weights.
- A monotonic mapping w'=tau+(1-tau)*w, using the existing analytically defined
  Fisher neutral point tau, keeping the entire joint topology. This checks the
  numerical domain of the aggregator and adds no fitted constant.

Reimport the original graph with unchanged weights as an identity control;
its final member TSV must be byte-identical to the baseline for each dataset.

## Representation and modality controls

The raw-frequency feature map also differs from RabbitBin's existing per-base
background-corrected PMH input. To avoid confusing that normalization change
with joint representation, repeat concat/whiten/CCA using exactly that
background correction (then probability normalization and square roots).
Record these as gc_concat/gc_whiten/gc_cca. Additionally run comp, gc_comp,
and coverage alone with the identical neighbour/weight/aggregation pipeline.
This determines whether a gain actually requires both input modalities.
These are common algorithmic ablations on all three datasets, never per-dataset
switches. All results, including regressions, remain in the report.
