#!/usr/bin/env python3
"""Construct label-free joint feature spaces and their matched graph controls."""

import argparse
import csv
import json
import struct
import time
from pathlib import Path

import faiss
import numpy as np
from scipy.linalg import eigh, svd
from sklearn.covariance import oas
from threadpoolctl import threadpool_limits

EDGE_DTYPE = np.dtype([("i", "<u4"), ("j", "<u4"), ("w", "<f4")])


def read_cache(path):
    with path.open("rb", buffering=0) as f:
        def pod(fmt):
            fmt = "<" + fmt
            return struct.unpack(fmt, f.read(struct.calcsize(fmt)))

        def vector(dtype):
            n, = pod("Q")
            result = np.fromfile(f, dtype=dtype, count=n)
            if len(result) != n:
                raise ValueError("truncated cache array")
            return result

        def strings():
            n, = pod("Q")
            return [f.read(pod("I")[0]).decode() for _ in range(n)]

        def matrix():
            rows, cols = pod("QQ")
            return np.fromfile(f, dtype="<f4", count=rows * cols).reshape(rows, cols)

        if f.read(8) != b"RBCACHE1" or pod("III") != (4, 4, 8):
            raise ValueError("requires RabbitBin cache v4 / float32 / size_t64")
        n, ns, dimensions, samples = pod("QQQQ")
        seed, total, total_short = pod("QQQ")
        pod("ddQQB")
        names, short_names = strings(), strings()
        lengths, short_lengths = vector("<u8"), vector("<u8")
        ranks, short_depth = matrix(), matrix()
        means, raw, colnorm = vector("<f4"), vector("<f4"), vector("<f8")
        source, target, pmh = vector("<u4"), vector("<u4"), vector("<f4")
        depth = means if means.size else raw
        if len(names) != n or depth.size != n * dimensions:
            raise ValueError("cache feature shape mismatch")
        return dict(names=names, lengths=lengths, depth=depth.reshape(n, dimensions),
                    samples=samples, seed=seed, source=source, target=target)


def load_composition(prefix, names, lengths, gc=False):
    with Path(str(prefix) + ".nodes.tsv").open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    lookup = {r["contig"]: (i, int(r["length_bp"])) for i, r in enumerate(rows)}
    if len(lookup) != len(rows):
        raise ValueError("duplicate FASTA names")
    order = np.empty(len(names), dtype=np.int64)
    for i, name in enumerate(names):
        order[i], length = lookup[name]
        if length != int(lengths[i]):
            raise ValueError("assembly/cache length mismatch: " + name)
    values = np.memmap(str(prefix) + ".f32", dtype="<f4", mode="r",
                       shape=(len(rows), 136))
    counts = np.array(values[order], dtype=np.float64)
    totals = counts.sum(axis=1, keepdims=True)
    counts = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0)
    if gc:
        base_counts = np.memmap(str(prefix) + ".base.f32", dtype="<f4", mode="r",
                                shape=(len(rows), 4))
        probs = np.array(base_counts[order], dtype=np.float64)
        probs /= np.maximum(probs.sum(axis=1, keepdims=True), 1.0)
        canonical = []
        for code in range(256):
            digits = [(code >> (2 * k)) & 3 for k in range(4)]
            rc = sum((3 - b) << (2 * (3-k)) for k, b in enumerate(digits))
            if code <= rc:
                canonical.append((code, rc, digits))
        for col, (code, rc, digits) in enumerate(canonical):
            expected = np.prod(probs[:, digits], axis=1)
            if code != rc:
                expected += np.prod(probs[:, [3-b for b in digits]], axis=1)
            counts[:, col] = np.divide(counts[:, col], expected,
                                       out=np.zeros(len(counts)), where=expected > 1e-14)
        # Same per-base background correction as RabbitBin's PMH input;
        # probability normalization permits a Hellinger feature map.
        counts /= np.maximum(counts.sum(axis=1, keepdims=True), np.finfo(float).tiny)
    return np.sqrt(counts)


def normalize_block(block):
    block -= block.mean(axis=0)
    scale = np.sqrt(np.square(block).sum() / len(block))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("constant or invalid feature block")
    return block / scale, float(scale)


def invsqrt(cov):
    eigenvalues, vectors = eigh(cov, check_finite=False)
    floor = np.finfo(np.float64).eps * len(cov) * max(float(eigenvalues[-1]), 1.0)
    return (vectors * (1.0 / np.sqrt(np.maximum(eigenvalues, floor)))) @ vectors.T


def representation(base, method):
    details = {}
    if method in ("comp", "coverage"):
        features = base[:, :136] if method == "comp" else base[:, 136:]
        return np.ascontiguousarray(features, dtype=np.float32), details
    if method == "concat":
        return np.asarray(base, dtype=np.float32), details
    covariance, shrinkage = oas(base, assume_centered=True)
    details["oas_shrinkage"] = float(shrinkage)
    if method == "whiten":
        embedding = base @ invsqrt(covariance)
    elif method == "cca":
        c = 136
        wc = invsqrt(covariance[:c, :c])
        wd = invsqrt(covariance[c:, c:])
        cross = wc @ covariance[:c, c:] @ wd
        u, correlation, vt = svd(cross, full_matrices=False, check_finite=False)
        correlation = np.clip(correlation, 0, 1)
        keep = correlation > np.finfo(np.float64).eps * max(cross.shape)
        correlation, u, vt = correlation[keep], u[:, keep], vt[keep]
        if not len(correlation):
            raise ValueError("no numerically correlated directions")
        embedding = ((base[:, :c] @ wc @ u) + (base[:, c:] @ wd @ vt.T))
        embedding *= np.sqrt(correlation / 2.0)
        details["canonical_correlations"] = correlation.tolist()
    else:
        raise ValueError(method)
    return np.ascontiguousarray(embedding, dtype=np.float32), details


def neighbours(embedding, k, threads):
    faiss.omp_set_num_threads(threads)
    index = faiss.IndexFlatL2(embedding.shape[1])
    index.add(embedding)
    distance, neighbour = index.search(embedding, k + 1)
    nonself = neighbour != np.arange(len(embedding))[:, None]
    keep = nonself & (np.cumsum(nonself, axis=1) <= k)
    return (np.maximum(distance[keep].reshape(-1, k), 0),
            neighbour[keep].reshape(-1, k).astype(np.uint32))


def mutual_edges(neighbour):
    n, k = neighbour.shape
    row = np.repeat(np.arange(n, dtype=np.uint32), k)
    col = neighbour.ravel()
    directed = np.sort(row.astype(np.uint64) * np.uint64(n) + col)
    keep = row < col
    row, col = row[keep], col[keep]
    reverse = col.astype(np.uint64) * np.uint64(n) + row
    pos = np.searchsorted(directed, reverse)
    mutual = (pos < len(directed))
    mutual[mutual] &= directed[pos[mutual]] == reverse[mutual]
    return row[mutual], col[mutual]


def read_edges(path):
    with path.open("rb") as f:
        magic, n, e = struct.unpack("<8sQQ", f.read(24))
        if magic != b"RBEDGE1\0":
            raise ValueError("bad graph magic")
        edges = np.fromfile(f, dtype=EDGE_DTYPE, count=e)
        if len(edges) != e or f.read(1):
            raise ValueError("graph length mismatch")
    return n, edges["i"].copy(), edges["j"].copy()


def affinity(embedding, scale, source, target):
    weights = np.empty(len(source), dtype=np.float32)
    floor = max(float(np.max(scale)), 1.0) * np.finfo(np.float64).eps
    for start in range(0, len(source), 8192):
        stop = min(start + 8192, len(source))
        i, j = source[start:stop], target[start:stop]
        delta = embedding[i].astype(np.float64) - embedding[j]
        distance = np.einsum("ij,ij->i", delta, delta)
        bandwidth = np.sqrt(np.maximum(scale[i], floor) * np.maximum(scale[j], floor))
        weights[start:stop] = np.exp(-distance / bandwidth)
    return np.maximum(weights, np.finfo(np.float32).tiny)


def write_graph(path, names, lengths, source, target, weights):
    n = len(names)
    order = np.argsort(source.astype(np.uint64) * np.uint64(n) + target)
    values = np.empty(len(source), dtype=EDGE_DTYPE)
    values["i"], values["j"], values["w"] = source[order], target[order], weights[order]
    with path.open("wb") as f:
        f.write(struct.pack("<8sQQ", b"RBEDGE1\0", n, len(values)))
        values.tofile(f)
    with Path(str(path) + ".nodes.tsv").open("w") as f:
        f.write("index\tcontig\tlength_bp\n")
        for i, (name, length) in enumerate(zip(names, lengths)):
            f.write(f"{i}\t{name}\t{length}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--retained", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--methods", nargs="+", choices=("concat", "whiten", "cca", "comp", "coverage",
                                                        "gc_concat", "gc_whiten", "gc_cca", "gc_comp"),
                        default=("concat", "whiten", "cca"))
    parser.add_argument("--threads", type=int, default=64)
    parser.add_argument("--neighbours", type=int, default=200)
    args = parser.parse_args()
    start = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    cache = read_cache(args.cache)
    composition = load_composition(args.features, cache["names"], cache["lengths"])
    coverage = cache["depth"].astype(np.float64)
    column_mean = coverage.mean(axis=0)
    coverage /= np.maximum(column_mean, np.finfo(np.float64).tiny)
    coverage = np.log1p(coverage)
    composition, cs = normalize_block(composition)
    coverage, ds = normalize_block(coverage)
    base = np.concatenate((composition, coverage), axis=1)
    del composition, coverage
    gc_base = None
    if any(m.startswith("gc_") for m in args.methods):
        gc_composition = load_composition(args.features, cache["names"], cache["lengths"], gc=True)
        gc_composition, gc_scale = normalize_block(gc_composition)
        gc_base = np.concatenate((gc_composition, base[:, 136:]), axis=1)
        del gc_composition
    n, original_i, original_j = read_edges(args.retained)
    if n != len(base):
        raise ValueError("retained graph/cache mismatch")
    preprocess_s = time.monotonic() - start
    for method in args.methods:
        begin = time.monotonic()
        with threadpool_limits(limits=args.threads):
            embedding, details = representation(gc_base if method.startswith("gc_") else base,
                                                 method.removeprefix("gc_"))
        projected = time.monotonic()
        print(f"{method}: embedding {embedding.shape}, model_s={projected-begin:.3f}", flush=True)
        # BLAS and OpenMP parallelism are explicitly bounded by the same budget.
        with threadpool_limits(limits=args.threads):
            distance, neighbour = neighbours(embedding, args.neighbours, args.threads)
        searched = time.monotonic()
        print(f"{method}: exact search_s={searched-projected:.3f}", flush=True)
        scale = distance[:, -1].astype(np.float64)
        source, target = mutual_edges(neighbour)
        del distance, neighbour
        weights = affinity(embedding, scale, source, target)
        write_graph(args.output / f"{method}.joint.rbedge", cache["names"],
                    cache["lengths"], source, target, weights)
        old_weights = affinity(embedding, scale, original_i, original_j)
        write_graph(args.output / f"{method}.retained.rbedge", cache["names"],
                    cache["lengths"], original_i, original_j, old_weights)
        details.update(method=method, nodes=n, dimensions=embedding.shape[1],
                       mutual_edges=len(source), original_edges=len(original_i),
                       model_s=projected-begin, search_s=searched-projected,
                       graph_s=time.monotonic()-searched,
                       total_s=time.monotonic()-begin, preprocess_s=preprocess_s,
                       composition_scale=gc_scale if method.startswith("gc_") else cs, coverage_scale=ds,
                       neighbours=args.neighbours, threads=args.threads,
                       cache=str(args.cache), features=str(args.features),
                       faiss=faiss.__version__, numpy=np.__version__)
        with (args.output / f"{method}.json").open("w") as f:
            json.dump(details, f, indent=2)
        print(f"{method}: edges={len(source)}, total_s={details['total_s']:.3f}", flush=True)


if __name__ == "__main__":
    main()
