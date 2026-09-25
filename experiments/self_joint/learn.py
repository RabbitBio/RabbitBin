#!/usr/bin/env python3
"""Positive-pair Mahalanobis learning; no gold, taxonomy or negative labels."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from sklearn.covariance import oas
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "joint_representation"))
from joint_space import (read_cache, load_composition, invsqrt, neighbours,
                         mutual_edges, affinity, read_edges, write_graph)


def rows(path):
    with path.open() as f: return list(csv.DictReader(f, delimiter="\t"))


def training_features(cache, whole_prefix, split_prefix):
    composition = load_composition(whole_prefix, cache["names"], cache["lengths"], gc=True)
    coverage = cache["depth"].astype(np.float64)
    colmean = coverage.mean(axis=0)
    colmean = np.maximum(colmean, np.finfo(np.float64).tiny)
    coverage = np.log1p(coverage / colmean)
    cc, dc = composition.mean(axis=0), coverage.mean(axis=0)
    cs = np.sqrt(np.mean(np.sum((composition-cc)**2, axis=1)))
    ds = np.sqrt(np.mean(np.sum((coverage-dc)**2, axis=1)))
    if not cs > 0 or not ds > 0: raise ValueError("constant feature block")
    whole = np.concatenate(((composition-cc)/cs, (coverage-dc)/ds), axis=1)
    del composition, coverage
    nodes = rows(Path(str(split_prefix) + ".depth.nodes.tsv"))
    samples = rows(Path(str(split_prefix) + ".depth.samples.tsv"))
    dims = 2 * len(samples)
    depth = np.fromfile(str(split_prefix) + ".depth.f32", dtype="<f4").reshape(len(nodes), 2, dims)
    if dims != cache["depth"].shape[1]: raise ValueError("depth dimensions differ")
    lookup = {name: i for i, name in enumerate(cache["names"])}
    if len(lookup) != len(cache["names"]): raise ValueError("duplicate parent name")
    kept, names, lengths, parent_names = [], [], [], []
    for r, node in enumerate(nodes):
        parent = node["contig"]
        if parent not in lookup: continue
        i = lookup[parent]
        if int(node["length_bp"]) != cache["lengths"][i]: raise ValueError("BAM/FASTA length mismatch")
        parent_names.append(parent)
        kept.append(r)
        for half in range(2):
            names.append(parent + "/fragment" + str(half))
            lengths.append(int(node["right_bp" if half else "left_bp"]))
    if len(kept) < 2: raise ValueError("too few split parents")
    comp = load_composition(split_prefix, names, lengths, gc=True)
    cov = np.log1p(depth[kept].reshape(-1, dims).astype(np.float64) / colmean)
    split = np.concatenate(((comp-cc)/cs, (cov-dc)/ds), axis=1).reshape(len(kept), 2, -1)
    if not np.isfinite(split).all() or not np.isfinite(whole).all(): raise ValueError("nonfinite feature")
    metadata = dict(parents=len(kept), bam_parents=len(nodes), dimensions=whole.shape[1],
                    composition_dims=136, coverage_dims=dims, gc_correction=True,
                    composition_scale=float(cs), coverage_scale=float(ds),
                    fragment_length_min=int(min(lengths)), fragment_length_median=float(np.median(lengths)),
                    all_zero_coverage_halves=int(np.sum(np.all(depth[kept] == 0, axis=2))))
    return whole, split, parent_names, metadata


def signal_factor(within, midpoint):
    # For left/right = latent + independent noise:
    # W = Cov(left-right), C = Cov((left+right)/2).
    # Random-parent difference covariance U = 2*C + W/2.
    # Gaussian same-vs-independent log likelihood orders by
    # delta.T [W^-1 - U^-1] delta. Project estimated signal onto the PSD
    # cone; lambda/(1+lambda) is derived, not a tuned dimension weight.
    white = invsqrt(within)
    signal = 2 * white @ midpoint @ white - 0.5*np.eye(len(within))
    value, vector = np.linalg.eigh((signal + signal.T)*0.5)
    value = np.maximum(value, 0.0)
    return white @ (vector * np.sqrt(value/(1.0+value)))


def transform(within, method, midpoint=None):
    if method.startswith("signal_"):
        if midpoint is None: raise ValueError("signal metric requires midpoint covariance")
        kind = method.removeprefix("signal_")
        if kind == "joint": return signal_factor(within, midpoint)
        if kind == "block":
            result = np.zeros_like(within)
            result[:136, :136] = signal_factor(within[:136, :136], midpoint[:136, :136])
            result[136:, 136:] = signal_factor(within[136:, 136:], midpoint[136:, 136:])
            return result
        if kind == "coverage":
            result = np.zeros((len(within), len(within)-136))
            result[136:, :] = signal_factor(within[136:, 136:], midpoint[136:, 136:])
            return result
        raise ValueError(method)
    if method == "joint": return invsqrt(within)
    if method == "block":
        result = np.zeros_like(within)
        result[:136, :136] = invsqrt(within[:136, :136])
        result[136:, 136:] = invsqrt(within[136:, 136:])
        return result
    if method == "coverage":
        result = np.zeros((len(within), len(within)-136))
        result[136:, :] = invsqrt(within[136:, 136:])
        return result
    raise ValueError(method)


def cross_validate(split, parents, methods):
    # Split by parent, never by half. Hash is independent of names' ordering.
    fold = np.array([hashlib.blake2b(name.encode(), digest_size=1).digest()[0] & 1
                     for name in parents])
    results = []
    for held in (0, 1):
        train, test = np.flatnonzero(fold != held), np.flatnonzero(fold == held)
        if len(train) < 2 or len(test) < 2: raise ValueError("insufficient fold size")
        delta_train = split[train, 0] - split[train, 1]
        within, shrinkage = oas(delta_train, assume_centered=True)
        midpoint, _ = oas(split[train].mean(axis=1))
        permutation = np.random.default_rng(42).permutation(test)
        other = np.empty_like(test)
        pairing = dict(zip(permutation.tolist(), np.roll(permutation, 1).tolist()))
        for j, i in enumerate(test): other[j] = pairing[int(i)]
        positive = split[test, 0] - split[test, 1]
        unlabeled = split[test, 0] - split[other, 1]
        for method in ("identity", *methods):
            matrix = np.eye(split.shape[2]) if method == "identity" else transform(within, method, midpoint)
            pd = np.sum((positive @ matrix)**2, axis=1)
            ud = np.sum((unlabeled @ matrix)**2, axis=1)
            results.append(dict(method=method, fold=held, train_parents=len(train),
                                test_parents=len(test), shrinkage=float(shrinkage),
                                median_same_parent_distance=float(np.median(pd)),
                                median_other_parent_distance=float(np.median(ud)),
                                same_closer_fraction=float(np.mean((pd < ud) + 0.5*(pd == ud)))))
    return results


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--features", type=Path, required=True)
    p.add_argument("--halves", type=Path, required=True)
    p.add_argument("--retained", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--threads", type=int, default=64)
    p.add_argument("--neighbours", type=int, default=200)
    p.add_argument("--methods", nargs="+", choices=("joint", "block", "coverage", "signal_joint", "signal_block", "signal_coverage"), default=("joint", "block", "coverage"))
    p.add_argument("--training-label", default="")
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    cache = read_cache(args.cache)
    whole, split, parents, details = training_features(cache, args.features, args.halves)
    prepared = time.monotonic()
    with threadpool_limits(limits=args.threads):
        crossval = cross_validate(split, parents, args.methods)
        validated = time.monotonic()
        delta = split[:, 0] - split[:, 1]
        within, shrinkage = oas(delta, assume_centered=True)
        midpoint, midpoint_shrinkage = oas(split.mean(axis=1))
    trained = time.monotonic()
    details.update(oas_shrinkage=float(shrinkage), feature_s=prepared-start,
                   fit_and_validation_s=trained-prepared, validation_s=validated-prepared,
                   fit_s=trained-validated, midpoint_shrinkage=float(midpoint_shrinkage), cv=crossval,
                   cross_block_frobenius=float(np.linalg.norm(within[:136, 136:])),
                   within_frobenius=float(np.linalg.norm(within)),
                   max_abs_half_bias=float(np.max(np.abs(delta.mean(axis=0)))))
    np.save(args.output / "within.npy", within)
    np.save(args.output / "midpoint.npy", midpoint)
    training_name = "training" + ("_" + args.training_label if args.training_label else "")
    with (args.output / f"{training_name}.json").open("w") as f: json.dump(details, f, indent=2)
    print(json.dumps(details), flush=True)
    n, oldi, oldj = read_edges(args.retained)
    if n != len(whole): raise ValueError("graph/cache node mismatch")
    for method in args.methods:
        begin = time.monotonic()
        with threadpool_limits(limits=args.threads):
            matrix = transform(within, method, midpoint)
            embedding = np.ascontiguousarray(whole @ matrix, dtype=np.float32)
        np.save(args.output / f"sj_{method}.transform.npy", matrix)
        projected = time.monotonic()
        print(f"{method}: projected {embedding.shape}", flush=True)
        with threadpool_limits(limits=args.threads):
            distance, neighbour = neighbours(embedding, args.neighbours, args.threads)
        searched = time.monotonic()
        scale = distance[:, -1].astype(np.float64)
        source, target = mutual_edges(neighbour)
        del neighbour, distance
        weights = affinity(embedding, scale, source, target)
        write_graph(args.output / f"sj_{method}.joint.rbedge", cache["names"], cache["lengths"], source, target, weights)
        old_weights = affinity(embedding, scale, oldi, oldj)
        write_graph(args.output / f"sj_{method}.retained.rbedge", cache["names"], cache["lengths"], oldi, oldj, old_weights)
        record = dict(method=method, nodes=n, dimensions=embedding.shape[1], mutual_edges=len(source),
                      project_s=projected-begin, search_s=searched-projected,
                      graph_s=time.monotonic()-searched, total_s=time.monotonic()-begin,
                      oas_shrinkage=float(shrinkage), neighbours=args.neighbours)
        with (args.output / f"sj_{method}.json").open("w") as f: json.dump(record, f, indent=2)
        print(json.dumps(record), flush=True)


if __name__ == "__main__": main()
