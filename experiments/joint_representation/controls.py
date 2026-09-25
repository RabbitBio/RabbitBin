#!/usr/bin/env python3
"""Matched topology/aggregation controls; no gold labels are read."""
import argparse
import struct
from pathlib import Path

import numpy as np

from joint_space import EDGE_DTYPE, read_cache, write_graph


def weighted(path):
    with path.open("rb") as f:
        magic, n, e = struct.unpack("<8sQQ", f.read(24))
        if magic != b"RBEDGE1\0":
            raise ValueError("bad graph header")
        values = np.fromfile(f, dtype=EDGE_DTYPE, count=e)
        if len(values) != e or f.read(1):
            raise ValueError("graph size mismatch")
    return n, values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--retained", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--methods", nargs="+", default=("concat", "whiten", "cca"))
    args = parser.parse_args()
    cache = read_cache(args.cache)
    n, baseline = weighted(args.retained)
    if n != len(cache["names"]):
        raise ValueError("node count mismatch")
    write_graph(args.output / "identity.rbedge", cache["names"], cache["lengths"],
                baseline["i"], baseline["j"], baseline["w"])
    # Analytic neutral point already used by RabbitBin: A_2(tau)=A_1(tau).
    # No fit/sweep: this only maps kernel affinities into Fisher's reinforcing
    # interval while preserving their ordering and the complete topology.
    tau = 0.7153318629591614
    for method in args.methods:
        _, joint = weighted(args.output / f"{method}.joint.rbedge")
        _, passed = weighted(args.output / f"{method}.candidates.scored.rbedge")
        keys = joint["i"].astype(np.uint64) * n + joint["j"]
        requested = passed["i"].astype(np.uint64) * n + passed["j"]
        order = np.argsort(keys)
        pos = np.searchsorted(keys[order], requested)
        if np.any(pos == len(keys)) or np.any(keys[order[pos]] != requested):
            raise ValueError("coverage-surviving edge is absent from imported topology")
        weights = joint["w"][order[pos]]
        write_graph(args.output / f"{method}.guarded.rbedge", cache["names"],
                    cache["lengths"], passed["i"], passed["j"], weights)
        weights = (tau + (1-tau) * joint["w"].astype(np.float64)).astype(np.float32)
        write_graph(args.output / f"{method}.neutral.rbedge", cache["names"],
                    cache["lengths"], joint["i"], joint["j"], weights)


if __name__ == "__main__":
    main()
