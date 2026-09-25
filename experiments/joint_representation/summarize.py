#!/usr/bin/env python3
"""Offline gold-labelled evaluation; never called by representation training."""
import argparse
import csv
import json
import re
import struct
from pathlib import Path

import numpy as np

from joint_space import EDGE_DTYPE, read_cache


def amber(path):
    text = path.read_text()
    result = {}
    for key, pattern in (
        ("HQ", r"^HQ .*: (\d+)$"),
        ("MQ", r"^MQ .*: (\d+)$"),
        ("LQ", r"^LQ .*: (\d+)$"),
        ("purity", r"^weighted purity .*: ([.\d]+)$"),
        ("completeness", r"^avg completeness/genome: ([.\d]+)$"),
    ):
        match = re.search(pattern, text, re.M)
        if not match:
            raise ValueError(f"missing {key}: {path}")
        result[key] = match[1]
    return result


def graph_stats(path, labels, lengths):
    with path.open("rb") as f:
        magic, n, e = struct.unpack("<8sQQ", f.read(24))
        if magic != b"RBEDGE1\0" or n != len(labels):
            raise ValueError("graph mismatch")
        edges = np.fromfile(f, dtype=EDGE_DTYPE, count=e)
    i, j = edges["i"], edges["j"]
    known = (labels[i] >= 0) & (labels[j] >= 0)
    same = known & (labels[i] == labels[j])
    degree = np.bincount(np.concatenate((i, j)), minlength=n)
    true_degree = np.bincount(np.concatenate((i[same], j[same])), minlength=n)
    represented = np.bincount(labels[labels >= 0])
    eligible = labels >= 0
    eligible[eligible] &= represented[labels[eligible]] >= 2
    denominator = lengths[eligible].sum()
    return dict(edges=e, true_edges=int(same.sum()), false_edges=int((known & ~same).sum()),
                same_neighbour_bp_frac=float(lengths[eligible & (true_degree > 0)].sum() / denominator),
                connected_bp_frac=float(lengths[eligible & (degree > 0)].sum() / denominator))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--datasets", nargs="+", default=("marine", "plant_associated", "strain_madness"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--methods", nargs="+", help="optional exact method names to include")
    args = parser.parse_args()
    output = []
    for dataset in args.datasets:
        base = args.root / "experiments/candidate_stage/cluster_results" / dataset
        result_dir = args.root / "experiments/joint_representation/results" / dataset
        cache = read_cache(base / "baseline.cache")
        gold = {}
        genome_ids = {}
        path = Path("/home/bigssd/zt/runs/cami2_benchmark") / dataset / "prep/gold_len.binning"
        with path.open() as f:
            for line in f:
                if line.startswith(("@", "#")):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 4:
                    continue
                genome_ids.setdefault(parts[1], len(genome_ids))
                gold[parts[0]] = genome_ids[parts[1]]
        labels = np.array([gold.get(name, -1) for name in cache["names"]], dtype=np.int32)
        entries = [("baseline", base / "baseline.amber", base / "retained.rbedge")]
        for f in sorted(result_dir.glob("*.amber")):
            name = f.name[:-6]
            entries.append((name, f, result_dir / f"{name}.scored.rbedge"))
        for name, report, graph in entries:
            if args.methods and name not in args.methods:
                continue
            if not report.exists():
                # Some original baseline summaries were printed but not saved;
                # use the identical diagnostic baseline's AMBER report.
                report = args.root / "experiments/candidate_stage/results_subgates" / dataset / "baseline.amber"
            row = dict(dataset=dataset, method=name, **amber(report),
                       **graph_stats(graph, labels, cache["lengths"]))
            output.append(row)
            print("\t".join(str(row[k]) for k in row), flush=True)
    if args.output and output:
        with args.output.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(output[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(output)


if __name__ == "__main__":
    main()
