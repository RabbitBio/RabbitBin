#!/usr/bin/env python3
"""Offline AMBER + gold-labelled graph audit; never a binner input."""
import argparse
import csv
from pathlib import Path
import re
import struct
import numpy as np
from scipy.stats import rankdata

EDGE = np.dtype([("i", "<u4"), ("j", "<u4"), ("w", "<f4")])


def read_pmh(path):
    """Read candidate composition scores; skip unrelated v4 cache arrays."""
    with path.open("rb") as f:
        def pod(fmt):
            return struct.unpack("<" + fmt, f.read(struct.calcsize("<" + fmt)))
        def vec(dtype):
            count, = pod("Q")
            array = np.fromfile(f, dtype=dtype, count=count)
            assert len(array) == count
            return array
        assert f.read(8) == b"RBCACHE1" and pod("III") == (4, 4, 8)
        n, _, _, _ = pod("QQQQ")
        pod("QQQddQQB")
        for _ in range(2):
            count, = pod("Q")
            for _ in range(count):
                length, = pod("I")
                f.seek(length, 1)
        for _ in range(2):
            count, = pod("Q")
            f.seek(8 * count, 1)
        for _ in range(2):
            nr, nc = pod("QQ")
            f.seek(4 * nr * nc, 1)
        for width in (4, 4, 8):
            count, = pod("Q")
            f.seek(width * count, 1)
        i, j, pmh = vec("<u4"), vec("<u4"), vec("<f4")
        assert len(i) == len(j) == len(pmh)
    keys = i.astype(np.uint64) * n + j
    order = np.argsort(keys)
    return keys[order], pmh[order]


def conditional_auc(values, same, known):
    y = same[known]
    pos = int(y.sum())
    neg = len(y) - pos
    if not pos or not neg: return float("nan")
    ranks = rankdata(values[known]).astype(np.float64, copy=False)
    return float((ranks[y].sum() - pos * (pos + 1) / 2) / (pos * neg))


def read_graph(path):
    with path.open("rb") as f:
        magic, n, e = struct.unpack("<8sQQ", f.read(24))
        assert magic == b"RBEDGE1\0"
        edges = np.fromfile(f, dtype=EDGE, count=e)
        assert len(edges) == e and not f.read(1)
    keys = edges["i"].astype(np.uint64) * n + edges["j"]
    order = np.argsort(keys)
    assert np.all(np.diff(keys[order]) > 0)
    return n, edges[order], keys[order]


def amber(path):
    text = path.read_text()
    result = {}
    for key, regex in (("HQ", r"^HQ .*: (\d+)$"), ("MQ", r"^MQ .*: (\d+)$"),
                       ("LQ", r"^LQ .*: (\d+)$"),
                       ("purity", r"^weighted purity .*: ([.\d]+)$"),
                       ("completeness", r"^avg completeness/genome: ([.\d]+)$")):
        result[key] = re.search(regex, text, re.M).group(1)
    timing = path.with_suffix(".time")
    if timing.exists():
        result.update(re.findall(r"(wall_s|peak_kb)=([.\d]+)", timing.read_text()))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    rows, grows = [], []
    for dataset in ("marine", "plant_associated", "strain_madness"):
        folder = args.root / dataset
        for path in sorted(folder.glob("*.amber")):
            row = dict(dataset=dataset, method=path.stem, **amber(path))
            rows.append(row)
            print(row, flush=True)
        basepath = folder / "baseline.r1.rbedge"
        newpath = folder / "coverage_first.r1.rbedge"
        if not all(p.exists() for p in (basepath, newpath,
                   folder / "baseline.r1.amber", folder / "coverage_first.r1.amber")): continue
        assert Path(str(basepath) + ".nodes.tsv").read_bytes() == Path(str(newpath) + ".nodes.tsv").read_bytes()
        n, base, bkeys = read_graph(basepath)
        nnew, new, nkeys = read_graph(newpath)
        assert n == nnew
        slots = np.searchsorted(nkeys, bkeys)
        assert np.all(slots < len(nkeys)) and np.array_equal(nkeys[slots], bkeys)
        assert np.array_equal(new["w"][slots], base["w"])
        labels, genome_ids = {}, {}
        goldpath = Path("/home/bigssd/zt/runs/cami2_benchmark") / dataset / "prep/gold_len.binning"
        with goldpath.open() as f:
            for line in f:
                if line.startswith(("@", "#")): continue
                parts = line.rstrip().split("\t")
                if len(parts) < 4: continue
                genome_ids.setdefault(parts[1], len(genome_ids))
                labels[parts[0]] = genome_ids[parts[1]]
        with Path(str(basepath) + ".nodes.tsv").open() as f:
            nodes = list(csv.DictReader(f, delimiter="\t"))
        lab = np.array([labels.get(r["contig"], -1) for r in nodes])
        length = np.array([int(r["length_bp"]) for r in nodes], dtype=np.int64)
        represented = np.bincount(lab[lab >= 0])
        eligible = lab >= 0
        eligible[eligible] &= represented[lab[eligible]] >= 2
        added = np.ones(len(new), dtype=bool)
        added[slots] = False
        bpmhk, bpmhv = read_pmh(folder / "baseline.r1.cache")
        npmhk, npmhv = read_pmh(folder / "coverage_first.r1.cache")
        bpos, npos = np.searchsorted(bpmhk, bkeys), np.searchsorted(npmhk, nkeys)
        assert np.array_equal(bpmhk[bpos], bkeys) and np.array_equal(npmhk[npos], nkeys)
        bpmh, npmh = bpmhv[bpos], npmhv[npos]
        assert np.array_equal(npmh[slots], bpmh)
        for name, edges, pmh in (("baseline", base, bpmh), ("coverage_first", new, npmh), ("added", new[added], npmh[added])):
            i, j = edges["i"], edges["j"]
            known = (lab[i] >= 0) & (lab[j] >= 0)
            same = known & (lab[i] == lab[j])
            degree = np.bincount(np.concatenate((i[same], j[same])), minlength=n)
            row = dict(dataset=dataset, graph=name, edges=len(edges),
                       true_edges=int(same.sum()), false_edges=int((known & ~same).sum()),
                       same_neighbour_bp_frac=float(length[eligible & (degree > 0)].sum() / length[eligible].sum()))
            for feature, values in (("pmh", pmh), ("coverage", edges["w"])):
                row[feature + "_median_same"] = float(np.median(values[same]))
                row[feature + "_median_cross"] = float(np.median(values[known & ~same]))
                row[feature + "_auc_selected_edges"] = conditional_auc(values, same, known)
            grows.append(row)
            print(row, flush=True)
    for name, data in (("summary.tsv", rows), ("graphs.tsv", grows)):
        if not data: continue
        with (args.root / name).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(data)


if __name__ == "__main__": main()
