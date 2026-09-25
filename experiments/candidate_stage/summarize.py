#!/usr/bin/env python3
"""Summarize diagnostic stages without using gold to change any binning call."""

import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def ratio(numerator, denominator):
    return numerator / denominator if denominator else 0.0


def summarize(directory):
    prefix = Path(directory) / "baseline"
    stage = list(read_tsv(Path(str(prefix) + ".candidate_stages.tsv")))
    assignments = {}
    for row in read_tsv(Path(str(prefix) + ".members.tsv")):
        assignments[row["SequenceName"]] = row["BinNum"]

    by_genome_bin = Counter()
    for row in stage:
        row["length_bp"] = int(row["length_bp"])
        for field in ("genome_large_contigs", "nonzero_same",
                      "coverage_feasible_same", "pmh_positive_same",
                      "own_top_same", "mutual_top_same", "candidate_same",
                      "retained_same", "foreign_in_top",
                      "true_component_id", "true_component_bp"):
            row[field] = int(row[field])
        for field in ("rho_ge_cut_same", "jcov_ge_cut_same", "both_ge_cut_same"):
            if field in row:
                row[field] = int(row[field])
        row["bin"] = assignments.get(row["contig"], "unbinned")
        if row["bin"] != "unbinned":
            by_genome_bin[(row["genome_id"], row["bin"])] += row["length_bp"]

    core_bin = {}
    for (genome, bin_id), bp in by_genome_bin.items():
        if genome not in core_bin or bp > core_bin[genome][1]:
            core_bin[genome] = (bin_id, bp)
    core_components = defaultdict(set)
    for row in stage:
        core = core_bin.get(row["genome_id"])
        if core and row["bin"] == core[0]:
            core_components[row["genome_id"]].add(row["true_component_id"])

    stages = Counter()
    bp = Counter()
    foreign_by_stage = defaultdict(list)
    mismatches = Counter()
    coverage_causes_bp = Counter()
    coverage_causes_n = Counter()
    eligible = [row for row in stage if row["genome_large_contigs"] >= 2]
    for row in eligible:
        genome = row["genome_id"]
        core = core_bin.get(genome)
        missed = core is not None and row["bin"] != core[0]
        if not core:
            label = "no_core_bin"
        elif not missed:
            label = "in_core"
        elif row["coverage_feasible_same"] == 0:
            label = "no_coverage_feasible_same"
        elif row["pmh_positive_same"] == 0:
            label = "no_positive_pmh_same"
        elif row["own_top_same"] == 0:
            label = "own_top_excluded"
        elif row["mutual_top_same"] == 0:
            label = "reciprocal_top_excluded"
        elif row["retained_same"] == 0:
            label = "mutual_but_no_retained"
        elif row["true_component_id"] not in core_components[genome]:
            label = "retained_but_disconnected_from_core"
        else:
            label = "connected_to_core_but_other_bin"

        stages[label] += 1
        bp[label] += row["length_bp"]
        if label == "no_coverage_feasible_same" and "rho_ge_cut_same" in row:
            if row["nonzero_same"] == 0:
                cause = "no_nonzero_pair"
            elif row["rho_ge_cut_same"] == 0 and row["jcov_ge_cut_same"] == 0:
                cause = "both_gates_fail"
            elif row["rho_ge_cut_same"] == 0:
                cause = "rho_gate_fail"
            elif row["jcov_ge_cut_same"] == 0:
                cause = "jcov_gate_fail"
            elif row["both_ge_cut_same"] == 0:
                cause = "gates_pass_on_different_peers"
            else:
                cause = "other_gate"
            coverage_causes_bp[cause] += row["length_bp"]
            coverage_causes_n[cause] += 1
        foreign_by_stage[label].append(row["foreign_in_top"])
        if row["mutual_top_same"] != row["retained_same"]:
            mismatches["feasible_mutual_vs_retained_nodes"] += 1
        if row["retained_same"] > row["candidate_same"]:
            mismatches["retained_gt_candidate_nodes"] += 1

    missed_bp = sum(v for k, v in bp.items() if k != "in_core")
    all_bp = sum(bp.values())
    print(f"\n{Path(directory).name}: eligible={len(eligible)} large_bp={all_bp:,} "
          f"missed_from_dominant_bin_bp={missed_bp:,} ({ratio(missed_bp, all_bp):.3f})")
    print("stage\tcontigs\tbp\tshare_of_missed_bp\tmedian_foreign_in_top")
    for label, count in stages.most_common():
        median = statistics.median(foreign_by_stage[label])
        share = "-" if label == "in_core" else f"{ratio(bp[label], missed_bp):.4f}"
        print(f"{label}\t{count}\t{bp[label]}\t{share}\t{median:g}")
    print("consistency:", dict(mismatches))
    if coverage_causes_bp:
        print("coverage subgate\tcontigs\tbp\tshare_of_no_coverage_bp")
        for cause, count in coverage_causes_n.most_common():
            print(f"{cause}\t{count}\t{coverage_causes_bp[cause]}\t"
                  f"{ratio(coverage_causes_bp[cause], bp['no_coverage_feasible_same']):.4f}")


if __name__ == "__main__":
    for item in sys.argv[1:]:
        summarize(item)
