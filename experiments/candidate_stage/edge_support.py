#!/usr/bin/env python3
"""Offline support check for MetaCAT-only HQ contigs on RabbitBin's fixed graph."""

import csv
import struct
import sys
from collections import Counter, defaultdict

from compare_hq import GOLD_ROOT, METACAT_ROOT, RABBIT_ROOT, STAGE_ROOT
from compare_hq import evaluate, load_gold, rows, stage_cause


EDGE = struct.Struct("<IIf")
HEADER = struct.Struct("<8sQQ")


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
    rb_hq, rb_counts, _, _ = evaluate(rabbit, gold, genome_bp)
    mc_hq, _, mc_contigs, _ = evaluate(meta, gold, genome_bp)
    stage = {
        row["contig"]: row
        for row in rows(STAGE_ROOT / dataset / "baseline.candidate_stages.tsv")
    }
    only_meta = set(mc_hq) - set(rb_hq)
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

    targets = {}
    for genome in only_meta:
        for bin_id in mc_hq[genome]:
            for contig in mc_contigs[bin_id]:
                if gold[contig][0] != genome or contig not in stage:
                    continue
                if rabbit.get(contig) == best_bin.get(genome, (None, 0))[0]:
                    continue
                if stage_cause(stage[contig], core_components[genome]) == \
                        "connected_to_core_but_other_bin":
                    targets[contig] = genome

    graph_path = RABBIT_ROOT / dataset / "retained.rbedge"
    names = []
    with (graph_path.with_name(graph_path.name + ".nodes.tsv")).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            names.append(row["contig"])
    target_indices = {i: name for i, name in enumerate(names) if name in targets}
    support = {name: [0.0, 0.0, 0, 0] for name in targets}
    with graph_path.open("rb") as handle:
        magic, n, e = HEADER.unpack(handle.read(HEADER.size))
        if magic != b"RBEDGE1\0" or n != len(names):
            raise ValueError("graph header mismatch")
        for _ in range(e):
            i, j, weight = EDGE.unpack(handle.read(EDGE.size))
            for target_index, neighbor_index in ((i, j), (j, i)):
                name = target_indices.get(target_index)
                if name is None:
                    continue
                genome = targets[name]
                other_bin = rabbit.get(names[neighbor_index])
                own_bin = rabbit.get(name)
                if other_bin == best_bin[genome][0]:
                    support[name][0] += weight
                    support[name][2] += 1
                if own_bin is not None and other_bin == own_bin:
                    support[name][1] += weight
                    support[name][3] += 1

    totals = Counter()
    for name, genome in targets.items():
        length = gold[name][1]
        core_w, own_w, core_edges, own_edges = support[name]
        totals["all_bp"] += length
        if core_edges:
            totals["direct_core_bp"] += length
        if rabbit.get(name) is None:
            totals["unbinned_bp"] += length
        elif core_w > own_w:
            totals["core_support_stronger_bp"] += length
        else:
            totals["own_support_at_least_core_bp"] += length
        if core_edges == 0 and own_edges == 0:
            totals["neither_direct_bp"] += length
    print(f"{dataset}: {len(targets)} contigs, {dict(totals)}")


if __name__ == "__main__":
    for dataset in sys.argv[1:] or ("marine", "plant_associated", "strain_madness"):
        main(dataset)
