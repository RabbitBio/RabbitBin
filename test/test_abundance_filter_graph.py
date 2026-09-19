"""Differential graph/output checks against the original per-pair dot path."""

import csv
import os
from pathlib import Path
import random
import subprocess
import sys


def main():
    binary = str(Path(sys.argv[1]).resolve())
    work = Path(sys.argv[2]).resolve()
    baseline_binary = str(Path(sys.argv[3]).resolve()) if len(sys.argv) > 3 else binary
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(741)
    # The production fused graph is selected only above 25000 contigs. Keep
    # 203 varied rows plus constant-rank rows to cross that boundary without
    # making the edge dump large. This covers partial batches and tile edges.
    count = 25003
    names = [f"contig_{i}" for i in range(count)]
    sequence = "".join(rng.choice("ACGT") for _ in range(2800))
    fasta = work / "assembly.fa"
    fasta.write_text("".join(f">{name}\n{sequence}\n" for name in names))
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("RABBIT_", "RB_", "OMP_"))}
    # Uncorrected PMH exposes identical-score neighbour ties. Disable the
    # pre-existing ISA-specific filters so this exercises the portable batch
    # filter on all test hosts, including ones with VNNI.
    env.update(RABBIT_PMH_BASE="0", RABBIT_ABD_RANK8="0", RABBIT_ABD_Q8="0")
    cases = 0
    for dims in (3, 46, 129):
        patterns = [[rng.randrange(1, 40) for _ in range(dims)] for _ in range(5)]
        depth = work / f"depth-{dims}.tsv"
        lines = ["\t".join(["contigName", "contigLen", "totalAvgDepth"] +
                           [v for k in range(dims) for v in (f"s{k}", f"s{k}-var")])]
        for i, name in enumerate(names):
            values = [max(0, x + rng.choice((-1, 0, 0, 1))) for x in patterns[i % 5]]
            if i % 17 == 0: values = [0] * dims
            if i % 19 == 0: values = [7] * dims
            if i % 11 == 0: values = patterns[0][:]  # exact depth/PMH ties
            if i >= 203: values = [7] * dims
            lines.append("\t".join([name, "2800", str(sum(values) / dims)] +
                                    [v for x in values for v in (str(x), "0")]))
        depth.write_text("\n".join(lines) + "\n")
        for cutoff, threads, edges in ((71.53318629591614, 4, 7), (99, 2, 200)):
            outputs = []
            for batched in (False, True):
                label = f"d{dims}-c{cutoff}-t{threads}-batch{int(batched)}"
                prefix = work / label
                dump = work / f"{label}.edges.tsv"
                run_env = dict(env, RB_PAIR_DUMP=str(dump))
                if not batched: run_env["RABBIT_NO_ABD_BATCH"] = "1"
                cmd = [binary if batched else baseline_binary,
                       "bin", "--fasta", str(fasta), "--depth", str(depth),
                       "--output", str(prefix), "--seed", "42", "--threads", str(threads),
                       "--max-edges", str(edges), "--min-edge-score", str(cutoff),
                       "--min-bin-size", "5600", "--no-bin-fasta"]
                result = subprocess.run(cmd, env=run_env, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, text=True, timeout=90)
                (work / f"{label}.log").write_text(result.stdout)
                if result.returncode:
                    raise AssertionError(f"{label}: {result.stdout}")
                if batched and "conservative batched dot" not in result.stdout:
                    raise AssertionError(f"{label}: batch filter was not exercised")
                with dump.open() as stream:
                    graph = sorted(tuple(sorted(row.items())) for row in
                                   csv.DictReader(stream, delimiter="\t"))
                with Path(str(prefix) + ".bins.tsv").open() as stream:
                    bins = [{k: v for k, v in row.items() if k != "FileName"}
                            for row in csv.DictReader(stream, delimiter="\t")]
                outputs.append((graph, bins, Path(str(prefix) + ".members.tsv").read_bytes()))
            if outputs[0] != outputs[1]:
                raise AssertionError(f"graph, bin stats or members differ: {label}")
            cases += 1
    print(f"{cases} graph/bin/member comparisons match the original dot path")


if __name__ == "__main__":
    main()
