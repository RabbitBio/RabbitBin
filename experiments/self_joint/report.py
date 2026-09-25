#!/usr/bin/env python3
"""Post-run, gold-labelled reporting only; not imported by training/binning."""
import argparse
import csv
import os
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/joint_representation"))
from joint_space import read_cache
from summarize import amber, graph_stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method-prefix", nargs="+", help="report only matching method prefixes plus baseline")
    parser.add_argument("--output", type=Path, help="summary path; default is OUTROOT/all_results.tsv")
    args = parser.parse_args()
    out = Path(os.environ.get("OUTROOT", str(ROOT / "experiments/self_joint/results")))
    base_root = Path(os.environ.get("BASE_ROOT", str(ROOT / "experiments/candidate_stage/cluster_results")))
    gold_root = Path(os.environ.get("GOLD_ROOT", "/home/bigssd/zt/runs/cami2_benchmark"))
    records = []
    for dataset in ("marine", "plant_associated", "strain_madness"):
        folder = out / dataset
        base = base_root / dataset
        cache = read_cache(base / "baseline.cache")
        gold, genome_ids = {}, {}
        with (gold_root / dataset / "prep/gold_len.binning").open() as f:
            for line in f:
                if line.startswith(("@", "#")): continue
                cols = line.rstrip().split("\t")
                if len(cols) < 4: continue
                genome_ids.setdefault(cols[1], len(genome_ids))
                gold[cols[0]] = genome_ids[cols[1]]
        labels = np.array([gold.get(name, -1) for name in cache["names"]], dtype=np.int32)
        report = ROOT / "experiments/candidate_stage/results_subgates" / dataset / "baseline.amber"
        entries = [("baseline", report, base / "retained.rbedge")]
        entries += [(p.stem, p, folder / (p.stem + ".scored.rbedge"))
                    for p in sorted(folder.glob("*.amber"))
                    if not args.method_prefix or p.stem.startswith(tuple(args.method_prefix))]
        for name, report, graph in entries:
            row = dict(dataset=dataset, method=name, **amber(report),
                       **graph_stats(graph, labels, cache["lengths"]))
            records.append(row)
            print("\t".join(str(v) for v in row.values()), flush=True)
    with (args.output or out / "all_results.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t")
        w.writeheader(); w.writerows(records)
    by = {(r["dataset"], r["method"]): r for r in records}
    baseline = {ds: by[ds, "baseline"] for ds in ("marine", "plant_associated", "strain_madness")}
    common = set.intersection(*[{r["method"] for r in records if r["dataset"] == ds} for ds in baseline])
    passes = []
    for name in sorted(common - {"baseline", "identity"}):
        if all(int(by[ds, name]["HQ"]) >= int(baseline[ds]["HQ"]) and
               float(by[ds, name]["purity"]) >= float(baseline[ds]["purity"]) for ds in baseline):
            passes.append(name)
    print("Non-identity configurations preserving HQ and weighted purity on all three:", passes, flush=True)


if __name__ == "__main__": main()
