"""BAM oracle for split coverage, including CIGARs, MAPQ, shards and fallback."""
from array import array
import csv
import os
from pathlib import Path
import random
import re
import struct
import subprocess
import sys


def main():
    binary = str(Path(sys.argv[1]).resolve())
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    refs = [(f"c{i}", (4999, 5000, 5001, 7000, 10000)[i % 5]) for i in range(30)]
    fa = work / "assembly.fa"
    with fa.open("w") as f:
        for name, length in refs:
            f.write(f">{name}\n" + "".join(rng.choice("ACGT") for _ in range(length)) + "\n")
    reads, bams = [], []
    patterns = (("100M", 0), ("40=2X58=", 2), ("40M2I58M", 2),
                ("50M2D50M", 2), ("50M10N50M", 0), ("100M", 30))
    for sample in range(3):
        sam, bam = work / f"s{sample}.sam", work / f"s{sample}.bam"
        with sam.open("w") as f:
            f.write("@HD\tVN:1.6\tSO:coordinate\n")
            for name, length in refs: f.write(f"@SQ\tSN:{name}\tLN:{length}\n")
            for tid, (name, length) in enumerate(refs):
                positions = sorted([0, 10, length//2-60, length//2, length-115] +
                                   [rng.randrange(length-115) for _ in range(400)])
                for read, pos in enumerate(positions):
                    cigar, nm = patterns[read % len(patterns)]
                    mq = (0, 3, 5, 60)[(read + sample + tid) % 4]
                    f.write(f"r{sample}_{tid}_{read}\t0\t{name}\t{pos+1}\t{mq}\t{cigar}\t*\t0\t0\t" +
                            "A"*100 + "\t" + "I"*100 + f"\tNM:i:{nm}\n")
                    if nm < 30: reads.append((sample, tid, pos, cigar, mq))
            for i in range(6):
                f.write(f"unmapped{i}\t4\t*\t0\t0\t*\t*\t0\t0\t" + "A"*100 + "\t" + "I"*100 + "\n")
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
        bams.append(bam)
    listing = work / "bam.list"
    listing.write_text("".join(str(p) + "\n" for p in bams))
    clean = {k: v for k, v in os.environ.items() if not k.startswith(("RABBIT_", "RB_", "OMP_"))}
    exports = []
    baseline = None
    for label, threads, shards, export in (("base", 4, 17, False), ("one", 1, 1, True),
                                            ("many", 4, 17, True), ("fallback", 4, 4096, True)):
        prefix = work / label
        env = dict(clean, RABBIT_DEPTH_B2_SHARDS=str(shards), RABBIT_FUSE_DUMP_DEPTH=str(prefix)+".whole.tsv")
        cmd = [binary, "bin", "--assembly", str(fa), "--bam-list", str(listing),
               "--output", str(prefix), "--threads", str(threads), "--seed", "42",
               "--dual-depth", "5", "--min-bin-size", "0", "--no-bin-fasta"]
        if export: cmd.extend(("--export-fragment-depth", str(prefix)))
        with Path(str(prefix)+".log").open("w") as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
        whole = (Path(str(prefix)+".whole.tsv").read_bytes(), Path(str(prefix)+".members.tsv").read_bytes())
        if baseline is None: baseline = whole
        assert whole == baseline, f"ordinary depth/members changed: {label}"
        if not export: continue
        with Path(str(prefix)+".depth.nodes.tsv").open() as f: nodes = list(csv.DictReader(f, delimiter="\t"))
        with Path(str(prefix)+".depth.samples.tsv").open() as f: samples = list(csv.DictReader(f, delimiter="\t"))
        indices = {row["contig"]: i for i, row in enumerate(nodes)}
        counts = [0] * (len(nodes) * 2 * 6)
        for sample, tid, pos, cigar, mq in reads:
            name, length = refs[tid]
            if name not in indices: continue
            edge = int(samples[sample]["edge_trim"])
            regions = ((edge, length//2-edge), (length//2+edge, length-edge))
            for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                size = int(size)
                if op in "M=X":
                    for half, (lo, hi) in enumerate(regions):
                        overlap = max(0, min(pos+size, hi)-max(pos, lo))
                        offset = (indices[name]*2+half)*6+sample
                        counts[offset] += overlap
                        if mq >= 5: counts[offset+3] += overlap
                if op in "MDN=X": pos += size
        expected = []
        for r, node in enumerate(nodes):
            for half in range(2):
                length = int(node["right_bp" if half else "left_bp"])
                for col in range(6):
                    edge = int(samples[col % 3]["edge_trim"])
                    expected.append(counts[(r*2+half)*6+col] / (length-2*edge))
        oracle = struct.pack("<"+"f"*len(expected), *expected)
        data = Path(str(prefix)+".depth.f32").read_bytes()
        assert data == oracle, f"region depth oracle mismatch: {label}"
        exports.append(data)
        print(f"PASS {label}: {len(nodes)} split parents, exact BAM region oracle, unchanged whole depth/members", flush=True)
    assert exports[0] == exports[1] == exports[2]
    log = (work / "fallback.log").read_text()
    assert "B2 re-sync failed on" in log, "fallback case was not exercised"
    print("PASS forced fallback clears partial fragment sums; 1/4 threads and shard counts identical")


if __name__ == "__main__": main()
