#!/usr/bin/env python3
"""Cluster an exported RabbitBin graph without gold labels or parameter search."""

import argparse
import csv
import random
import struct
import time
from pathlib import Path

import igraph


EDGE = struct.Struct("<IIf")
HEADER = struct.Struct("<8sQQ")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=Path)
    parser.add_argument("--method", choices=("infomap", "leiden"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    start = time.monotonic()
    with args.graph.open("rb") as handle:
        magic, n, e = HEADER.unpack(handle.read(HEADER.size))
        if magic != b"RBEDGE1\0":
            raise ValueError("unexpected graph format")
        raw = handle.read()
    if len(raw) != e * EDGE.size:
        raise ValueError("truncated graph")
    edges = []
    weights = []
    for i, j, weight in EDGE.iter_unpack(raw):
        edges.append((i, j))
        weights.append(weight)
    graph = igraph.Graph(n=n, edges=edges, directed=False)
    igraph.set_random_number_generator(random.Random(42))
    loaded = time.monotonic()
    if args.method == "infomap":
        partition = graph.community_infomap(edge_weights=weights)
    else:
        # Standard weighted modularity, converged without a tuned resolution.
        partition = graph.community_leiden(
            objective_function="modularity", weights=weights, n_iterations=-1
        )
    clustered = time.monotonic()
    membership = partition.membership
    nodes = args.graph.with_name(args.graph.name + ".nodes.tsv")
    with nodes.open(newline="") as handle, args.output.open("w") as out:
        rows = csv.DictReader(handle, delimiter="\t")
        for index, row in enumerate(rows):
            if int(row["index"]) != index:
                raise ValueError("node order mismatch")
            out.write(f"{index}\t{row['contig']}\t{membership[index]}\n")
    print(
        f"method={args.method} nodes={n} edges={e} communities={len(partition)} "
        f"load_s={loaded-start:.2f} cluster_s={clustered-loaded:.2f} "
        f"total_s={time.monotonic()-start:.2f}"
    )


if __name__ == "__main__":
    main()
