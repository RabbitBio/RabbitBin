"""Synthetic ordering / unchanged-weight / deterministic-graph regression."""
import csv
import os
from pathlib import Path
import random
import struct
import subprocess
import sys


def graph(path):
    data = path.read_bytes()
    magic, nodes, count = struct.unpack_from("<8sQQ", data)
    assert magic == b"RBEDGE1\0" and len(data) == 24 + 12 * count
    return {(i, j): w for i, j, w in struct.iter_unpack("<IIf", data[24:])}


def main():
    binary = str(Path(sys.argv[1]).resolve())
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    seq = "".join(rng.choice("ACGT") for _ in range(2800))
    assembly = work / "assembly.fa"
    depth = work / "depth.tsv"
    with assembly.open("w") as fa, depth.open("w") as dep:
        dep.write("contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\ts2\ts2-var\ts3\ts3-var\n")
        for i in range(25003):
            fa.write(f">c{i}\n{seq}\n")
            # c0's first 200 PMH ties have incompatible magnitude; its later
            # neighbours match both shape and magnitude. Other rows are inert.
            vals = [10, 20, 30] if i == 0 or 203 <= i < 406 else \
                   [100, 200, 300] if i < 203 else [7, 7, 7]
            fields = [f"c{i}", "2800", str(sum(vals) / 3)]
            for v in vals: fields.extend((str(v), "0"))
            dep.write("\t".join(fields) + "\n")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("RABBIT_", "RB_", "OMP_"))}
    env["RABBIT_PMH_BASE"] = "0"
    graphs = []
    for label, threads, early in (("base", 4, False), ("early1", 1, True), ("early4", 4, True)):
        prefix = work / label
        run_env = dict(env)
        if early: run_env["RABBIT_CANDIDATE_COVERAGE"] = "1"
        cmd = [binary, "bin", "--fasta", str(assembly), "--depth", str(depth),
               "--output", str(prefix), "--seed", "42", "--threads", str(threads),
               "--max-edges", "200", "--min-bin-size", "5600", "--no-bin-fasta",
               "--export-retained-graph", str(prefix) + ".rbedge"]
        with Path(str(prefix) + ".log").open("w") as log:
            subprocess.run(cmd, env=run_env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=120)
        graphs.append(graph(Path(str(prefix) + ".rbedge")))
    base, early1, early4 = graphs
    assert (0, 203) not in base and (0, 203) in early1
    assert all(early1.get(pair) == weight for pair, weight in base.items())
    assert early1 == early4
    assert (work / "early1.members.tsv").read_bytes() == (work / "early4.members.tsv").read_bytes()
    print(f"candidate regression passed: {len(base)} -> {len(early1)} edges; "
          "old weighted graph preserved, matching later neighbour recovered, 1/4 threads identical")


if __name__ == "__main__":
    main()
