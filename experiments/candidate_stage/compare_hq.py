#!/usr/bin/env python3
"""Inspect contigs that MetaCAT recovers in HQ bins but RabbitBin misses.

Gold is used only for this offline evaluation, never for a clustering call.
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
STAGE_ROOT = ROOT / "RabbitBin_stage_audit/experiments/candidate_stage/results_subgates"
RABBIT_ROOT = ROOT / "RabbitBin_stage_audit/experiments/candidate_stage/cluster_results"
METACAT_ROOT = ROOT / "MetaCAT/benchmark_cami2/runs_from_bam"
GOLD_ROOT = Path("/home/bigssd/zt/runs/cami2_benchmark")


def rows(path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def load_gold(path):
    gold = {}
    genome_bp = Counter()
    with path.open() as handle:
        for line in handle:
            if not line or line[0] in "@#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            contig, genome, _, length = fields[:4]
            length = int(length)
            gold[contig] = (genome, length)
            genome_bp[genome] += length
    return gold, genome_bp


def evaluate(assignments, gold, genome_bp):
    by_bin_genome = defaultdict(Counter)
    by_bin_contigs = defaultdict(list)
    unknown = 0
    for contig, bin_id in assignments.items():
        if contig not in gold:
            unknown += 1
            continue
        genome, length = gold[contig]
        by_bin_genome[bin_id][genome] += length
        by_bin_contigs[bin_id].append(contig)
    hq = defaultdict(set)
    for bin_id, counts in by_bin_genome.items():
        size = sum(counts.values())
        if size < 200_000:
            continue
        genome, true_bp = max(counts.items(), key=lambda x: x[1])
        if true_bp / size > 0.95 and true_bp / genome_bp[genome] > 0.90:
            hq[genome].add(bin_id)
    return hq, by_bin_genome, by_bin_contigs, unknown


def stage_cause(row, core_components):
    if row is None:
        return "short_not_in_graph"
    if int(row["coverage_feasible_same"]) == 0:
        return "no_coverage_feasible_same"
    if int(row["pmh_positive_same"]) == 0:
        return "no_positive_pmh_same"
    if int(row["own_top_same"]) == 0:
        return "own_top_excluded"
    if int(row["mutual_top_same"]) == 0:
        return "reciprocal_top_excluded"
    if int(row["retained_same"]) == 0:
        return "mutual_but_no_retained"
    if row["true_component_id"] not in core_components:
        return "retained_but_disconnected_from_core"
    return "connected_to_core_but_other_bin"


def main(dataset):
    gold, genome_bp = load_gold(GOLD_ROOT / dataset / "prep/gold_len.binning")
    rabbit = {
        row["SequenceName"]: row["BinNum"]
        for row in rows(RABBIT_ROOT / dataset / "baseline.members.tsv")
    }
    meta = {
        row["Sequence ID"]: row["Cluster ID"]
        for row in rows(METACAT_ROOT / dataset / "metacat.mapping")
    }
    rb_hq, rb_counts, _, rb_unknown = evaluate(rabbit, gold, genome_bp)
    mc_hq, _, mc_contigs, mc_unknown = evaluate(meta, gold, genome_bp)
    stage = {
        row["contig"]: row
        for row in rows(STAGE_ROOT / dataset / "baseline.candidate_stages.tsv")
    }
    if rb_unknown or mc_unknown:
        raise ValueError(f"unlabelled predicted contigs: Rabbit={rb_unknown}, MetaCAT={mc_unknown}")
    rb_hq_n = sum(map(len, rb_hq.values()))
    mc_hq_n = sum(map(len, mc_hq.values()))
    print(f"{dataset}: HQ bins Rabbit={rb_hq_n}, MetaCAT={mc_hq_n}; "
          f"HQ genomes Rabbit={len(rb_hq)}, MetaCAT={len(mc_hq)}")
    only_meta = set(mc_hq) - set(rb_hq)
    only_rabbit = set(rb_hq) - set(mc_hq)
    print(f"  MetaCAT-only HQ genomes={len(only_meta)}, "
          f"Rabbit-only HQ genomes={len(only_rabbit)}")

    best_bin = {}
    for bin_id, counts in rb_counts.items():
        for genome, value in counts.items():
            if genome not in best_bin or value > best_bin[genome][1]:
                best_bin[genome] = (bin_id, value)
    core_components = defaultdict(set)
    for contig, row in stage.items():
        genome = gold[contig][0]
        if genome in best_bin and rabbit.get(contig) == best_bin[genome][0]:
            core_components[genome].add(row["true_component_id"])

    bp = Counter()
    contigs = Counter()
    assignment_bp = Counter()
    seen = set()
    for genome in only_meta:
        for bin_id in mc_hq[genome]:
            for contig in mc_contigs[bin_id]:
                if contig in seen or gold[contig][0] != genome:
                    continue
                seen.add(contig)
                if genome in best_bin and rabbit.get(contig) == best_bin[genome][0]:
                    continue
                cause = stage_cause(stage.get(contig), core_components[genome])
                bp[cause] += gold[contig][1]
                contigs[cause] += 1
                assignment = "other_bin" if contig in rabbit else "unbinned"
                assignment_bp[(cause, assignment)] += gold[contig][1]
    total = sum(bp.values())
    print(f"  MetaCAT-HQ true contigs outside RabbitBin best bin: {total:,} bp")
    for cause, value in bp.most_common():
        print(f"    {cause}: {contigs[cause]} contigs, {value:,} bp, "
              f"{value / total:.1%}; "
              f"other_bin={assignment_bp[(cause, 'other_bin')]:,}, "
              f"unbinned={assignment_bp[(cause, 'unbinned')]:,}" if total else "")


if __name__ == "__main__":
    for dataset in sys.argv[1:] or ("marine", "plant_associated", "strain_madness"):
        main(dataset)
