"""Check multi-sample coverage weights and low-sample mean coverage ratios."""

import csv
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import zlib


DEFAULT_EDGE_CUTOFF = 0.7153318629591614

PROFILES = {
    "anchor": (10, 20, 30, 40),
    "identical": (10, 20, 30, 40),
    "scaled": (20, 40, 60, 80),
    "above_cut": (13, 27, 41, 55),
    "just_below_cut": (14, 28, 42, 56),
    "below_cut": (15, 30, 45, 60),
    "swapped": (10, 20, 40, 30),
    "reverse": (40, 30, 20, 10),
    "tied": (10, 10, 40, 40),
    "constant": (25, 25, 25, 25),
    "zero": (0, 0, 0, 0),
    "near_left": (10, 9, 20, 30),
    "near_right": (9, 10, 20, 30),
    "ratio_left": (100, 10, 20, 30),
    "ratio_right": (90, 5, 20, 30),
    "sparse_left": (10, 0, 20, 30),
    "sparse_right": (9, 0, 20, 30),
}


def ranks(values):
    return [sum(v < x for v in values) + (sum(v == x for v in values) + 1) / 2
            for x in values]


def spearman(a, b):
    a, b = ranks(a), ranks(b)
    a = [x - sum(a) / len(a) for x in a]
    b = [x - sum(b) / len(b) for x in b]
    scale = math.sqrt(sum(x * x for x in a) * sum(x * x for x in b))
    return sum(x * y for x, y in zip(a, b)) / scale if scale else 0.0


def near(actual, expected, label):
    if not math.isfinite(actual) or abs(actual - expected) > 2e-6:
        raise AssertionError(f"{label}: got {actual}, expected {expected}")


def check_edges(rows, samples, power=1.0, cutoff=DEFAULT_EDGE_CUTOFF,
                feature_profiles=None):
    profiles = (feature_profiles if feature_profiles is not None else
                {name: p[:samples] for name, p in PROFILES.items()})
    dimensions = len(next(iter(profiles.values())))
    means = [sum(p[s] for p in profiles.values()) / len(profiles) or 1.0
             for s in range(dimensions)]
    if not rows:
        raise AssertionError("fixture produced no candidate edges")
    for key, row in rows.items():
        a, b = (profiles[row[n]] for n in ("name_i", "name_j"))
        rho = spearman(a, b)
        denominator = sum(max(x, y) / m for x, y, m in zip(a, b, means))
        jcov = (sum(min(x, y) / m for x, y, m in zip(a, b, means)) /
                denominator) if denominator else 0.0
        raw = min(max(rho, 0.0), jcov)
        if samples >= 3:
            weight = raw ** power if raw >= cutoff else 0.0
            near(float(row["rho"]), rho, f"{key} rho")
        else:
            informative = [min(x, y) / max(x, y) for x, y in zip(a, b)
                           if max(x, y) > 0]
            raw = sum(informative) / len(informative) if informative else 0.0
            weight = raw ** power if raw >= cutoff else 0.0
            near(float(row["Amean"]), raw, f"{key} mean coverage ratio")
            if not math.isnan(float(row["rho"])):
                raise AssertionError("low-sample path must not calculate correlation")
        near(float(row["Jcov"]), jcov, f"{key} Jcov")
        near(float(row["dterm"]), raw, f"{key} raw weight")
        near(float(row["weight"]), weight, f"{key} scored weight")

    if feature_profiles is not None:
        return

    if samples <= 2:
        scaled = rows[tuple(sorted(("anchor", "scaled")))]
        near(float(scaled["weight"]), 0.0,
             "PMH similarity must not rescue mismatched coverage magnitudes")
        opposed = rows[tuple(sorted(("near_left", "near_right")))]
        near(float(opposed["weight"]), 0.9 ** power,
             "opposite ranks must not reject a close low-sample match")
        sparse = rows[tuple(sorted(("sparse_left", "sparse_right")))]
        near(float(sparse["weight"]), 0.9 ** power,
             "joint absence must not penalize a matching informative sample")
        if samples == 2:
            example = rows[tuple(sorted(("ratio_left", "ratio_right")))]
            near(float(example["dterm"]), 0.7, "equal sample weighting")
        return

    # Equal rank patterns alone cannot rescue mismatched coverage magnitudes.
    # Conversely, identical coverage must give weight 1 even if PMH differs.
    for other, positive in [("identical", True), ("above_cut", True),
                            ("just_below_cut", 5 / 7 >= cutoff),
                            ("below_cut", False), ("scaled", False),
                            ("reverse", False), ("constant", False),
                            ("zero", False)]:
        row = rows[tuple(sorted(("anchor", other)))]
        if (float(row["weight"]) > 0) != positive:
            raise AssertionError(f"unexpected survival for anchor/{other}: {row}")


def write_bam(path, sequence, sample):
    """Full-length alignments give known constant depth; MAPQ adds a distinct block."""
    names = list(PROFILES)
    length = len(sequence)
    header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name in names)
    data = bytearray(b"BAM\x01" + struct.pack("<i", len(header)) + header.encode())
    data += struct.pack("<i", len(names))
    for name in names:
        encoded = name.encode() + b"\0"
        data += struct.pack("<i", len(encoded)) + encoded + struct.pack("<i", length)
    alphabet = {"A": 1, "C": 2, "G": 4, "T": 8}
    packed = bytes((alphabet[sequence[i]] << 4) | alphabet[sequence[i + 1]]
                   for i in range(0, length, 2))
    for ref, name in enumerate(names):
        for read in range(PROFILES[name][sample]):
            qname = f"r{ref}_{read}".encode() + b"\0"
            mapq = 60 if read % 3 == 0 else 0
            core = struct.pack("<iiIIiiii", ref, 0,
                               (4681 << 16) | (mapq << 8) | len(qname),
                               1, length, -1, -1, 0)
            record = (core + qname + struct.pack("<I", length << 4) + packed +
                      bytes([30]) * length + b"NMi" + struct.pack("<i", 0))
            data += struct.pack("<i", len(record)) + record
    with path.open("wb") as stream:
        # A BGZF stream consists of independently compressed blocks and EOF.
        chunks = [data[i:i + 60000] for i in range(0, len(data), 60000)] + [b""]
        for chunk in chunks:
            encoder = zlib.compressobj(wbits=-15)
            payload = encoder.compress(chunk) + encoder.flush()
            header = bytes.fromhex("1f8b08040000000000ff060042430200")
            stream.write(header + struct.pack("<H", len(payload) + 25) + payload +
                         struct.pack("<II", zlib.crc32(chunk), len(chunk)))


def main():
    binary = Path(sys.argv[1]).resolve()
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("RABBIT_", "RB_", "OMP_"))}
    # Uncorrected PMH keeps the small fixture's candidate pairs visible,
    # including negative and zero coverage support.
    env["RABBIT_PMH_BASE"] = "0"
    rng = random.Random(42)
    sequence = "".join(rng.choice("ACGT") for _ in range(3200))
    fasta = work / "same.fa"
    fasta.write_text("".join(f">{name}\n{sequence}\n" for name in PROFILES))
    varied = work / "varied.fa"
    varied.write_text("".join(
        f">{name}\n" + "".join(rng.choice("ACGT") for _ in range(3200)) + "\n"
        for name in PROFILES))
    for samples in (1, 2, 3, 4):
        header = ["contigName", "contigLen", "totalAvgDepth"]
        for s in range(samples):
            header += [f"s{s}", f"s{s}-var"]
        lines = ["\t".join(header)]
        for name, profile in PROFILES.items():
            values = list(profile[:samples])
            fields = [name, "3200", str(sum(values) / samples)]
            for value in values:
                fields += [str(value), "0"]
            lines.append("\t".join(fields))
        (work / f"depth{samples}.tsv").write_text("\n".join(lines) + "\n")

    def run(label, samples=4, assembly=fasta, extra=(), overrides=None,
            cache=None, expect_ok=True, threads=1, bams=None,
            expected_error="alpha:power mixing was removed"):
        prefix = work / label
        dump = work / f"{label}.pairs.tsv"
        cmd = [str(binary), "bin", "--output", str(prefix), "--seed", "42",
               "--threads", str(threads), "--min-bin-size", "0", "--no-recruit",
               "--no-split", "--no-bin-fasta"]
        if cache:
            cmd += ["--load-cache", str(cache)]
        elif bams:
            cmd += ["--assembly", str(assembly)]
            for bam in bams:
                cmd += ["--bam", str(bam)]
        else:
            cmd += ["--assembly", str(assembly), "--depth",
                    str(work / f"depth{samples}.tsv")]
        cmd += list(extra)
        run_env = dict(env, RB_PAIR_DUMP=str(dump))
        run_env.update(overrides or {})
        completed = subprocess.run(cmd, env=run_env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, timeout=60)
        (work / f"{label}.log").write_text(completed.stdout)
        if expect_ok != (completed.returncode == 0):
            raise AssertionError(f"{label}: exit {completed.returncode}\n{completed.stdout}")
        if not expect_ok:
            if expected_error not in completed.stdout:
                raise AssertionError(f"{label}: unexpected error\n{completed.stdout}")
            return {}, completed.stdout
        with dump.open() as stream:
            rows = {tuple(sorted((r["name_i"], r["name_j"]))): r
                    for r in csv.DictReader(stream, delimiter="\t")}
        return rows, completed.stdout

    cache = work / "coverage.cache"
    baseline, _ = run("four", extra=("--save-cache", str(cache)))
    check_edges(baseline, 4)
    explicit, _ = run("explicit_cutoff", extra=(
        "--min-edge-score", repr(100 * DEFAULT_EDGE_CUTOFF)))
    if explicit != baseline:
        raise AssertionError("default and explicit fractional cutoff differ")
    lower, _ = run("lower_cutoff", extra=("--min-edge-score", "70"))
    check_edges(lower, 4, cutoff=0.7)
    # This magnitude ratio lies between 0.70 and the Fisher-derived default.
    boundary = tuple(sorted(("anchor", "just_below_cut")))
    if float(baseline[boundary]["weight"]) != 0 or float(lower[boundary]["weight"]) <= 0:
        raise AssertionError("fractional default and explicit override must separate the boundary pair")
    varied_rows, _ = run("varied", assembly=varied)
    check_edges(varied_rows, 4)
    if baseline.keys() != varied_rows.keys():
        raise AssertionError("fixture must keep the same candidate pairs")
    if not any(abs(float(baseline[k]["sComp"]) - float(varied_rows[k]["sComp"])) > 0.01
               for k in baseline):
        raise AssertionError("fixture did not change composition similarities")
    for k in baseline:
        near(float(varied_rows[k]["weight"]), float(baseline[k]["weight"]),
             f"{k}: weight must not depend on composition")

    for label, kwargs in [
        ("one", {"samples": 1}),
        ("one_varied", {"samples": 1, "assembly": varied}),
        ("two", {"samples": 2}),
        ("two_varied", {"samples": 2, "assembly": varied}),
        ("three", {"samples": 3}),
        ("threads", {"threads": 4}),
        ("legacy_weight", {"overrides": {"RABBIT_W_COMP": "0.75"}}),
        ("cached", {"cache": cache}),
        ("reuse", {"overrides": {"RABBIT_REUSE_SWEEP": "1.0"}}),
        ("reuse_power", {"overrides": {"RABBIT_REUSE_SWEEP": "2.0"}}),
        ("certify", {"extra": ("--certify",)}),
    ]:
        rows, log = run(label, **kwargs)
        check_edges(rows, kwargs.get("samples", 4), 2.0 if label == "reuse_power" else 1.0)
        if label == "legacy_weight" and "RABBIT_W_COMP is ignored" not in log:
            raise AssertionError("legacy composition override must be reported")
        if label in ("threads", "legacy_weight", "cached", "reuse", "certify"):
            if rows != baseline:
                raise AssertionError(f"{label} changed the scored edge table")

    for samples in (1, 2):
        low_cache = work / f"low{samples}.cache"
        direct, log = run(f"low{samples}_save", samples=samples,
                          extra=("--save-cache", str(low_cache)))
        restored, _ = run(f"low{samples}_restore", cache=low_cache)
        if direct != restored:
            raise AssertionError("low-sample cache changed the scored edges")
        if "Calculated spearman" in log or "Abundance-first prune:" in log:
            raise AssertionError("low-sample build used a correlation path")
        for suffix, options in [
            ("threads", {"threads": 4}),
            ("reuse", {"overrides": {"RABBIT_REUSE_SWEEP": "1.0"}}),
            ("certify", {"extra": ("--certify",)})
        ]:
            checked, _ = run(f"low{samples}_{suffix}", samples=samples, **options)
            if direct != checked:
                raise AssertionError(f"low-sample {suffix} changed the scored edges")

    # Adding a second sample must update the branch without ranking its values.
    # The appended profile duplicates sample 0, so all ratio weights stay equal.
    one, _ = run("incremental_one", samples=1)
    added, log = run("incremental_two", cache=work / "low1.cache",
                     extra=("--add-depth", str(work / "depth1.tsv")))
    if "Edge weighting: S=2; mean coverage-feature ratio" not in log:
        raise AssertionError("incremental depth lost the actual sample count")
    for key in one:
        near(float(added[key]["weight"]), float(one[key]["weight"]),
             "incremental low-sample ratio")

    # Moving from one to three samples must rebuild normalization and ranks.
    added, log = run("incremental_three", cache=work / "low1.cache",
                     extra=("--add-depth", str(work / "depth2.tsv")))
    if "Edge weighting: S=3; coverage only" not in log:
        raise AssertionError("incremental depth did not enter the multi-sample branch")
    check_edges(added, 3, feature_profiles={
        name: (p[0], p[0], p[1]) for name, p in PROFILES.items()})

    old_cache = work / "old.cache"
    contents = bytearray(cache.read_bytes())
    struct.pack_into("<I", contents, 8, 3)
    old_cache.write_bytes(contents)
    run("old_cache", cache=old_cache, expect_ok=False,
        expected_error="Rebuild the cache from the original inputs")

    help_text = subprocess.run([str(binary), "bin", "--help"],
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, timeout=60).stdout
    if "--bam" in help_text:
        bams = [work / f"sample{s}.bam" for s in (0, 1)]
        for s, bam in enumerate(bams):
            write_bam(bam, sequence, s)
        for samples in (1, 2):
            bam_cache = work / f"bam{samples}.cache"
            direct, log = run(f"bam{samples}", bams=bams[:samples],
                              extra=("--save-cache", str(bam_cache)))
            features = {name: list(p[:samples]) +
                        [math.ceil(p[s] / 3) for s in range(samples)]
                        for name, p in PROFILES.items()}
            check_edges(direct, samples, feature_profiles=features)
            opposed = direct[tuple(sorted(("near_left", "near_right")))]
            near(float(opposed["weight"]), 0.825,
                 "ordinary and MAPQ-filtered features must both contribute")
            if f"Edge weighting: S={samples}; mean coverage-feature ratio" not in log:
                raise AssertionError("BAM sample count did not select mean coverage ratios")
            if f"coverage columns={2 * samples};" not in log:
                raise AssertionError("BAM fixture did not exercise dual-depth columns")
            if "Calculated spearman" in log or "Abundance-first prune:" in log:
                raise AssertionError("low-sample BAMs used a correlation path")
            restored, _ = run(f"bam{samples}_restore", cache=bam_cache)
            if restored != direct:
                raise AssertionError("BAM cache lost input sample/coverage column distinction")
            ordinary, ordinary_log = run(f"bam{samples}_ordinary", bams=bams[:samples],
                                         extra=("--dual-depth", "0"))
            check_edges(ordinary, samples)
            if f"coverage columns={samples};" not in ordinary_log:
                raise AssertionError("disabling dual depth must keep one feature per sample")
            if samples == 1:
                added, log = run("bam_incremental_two", cache=bam_cache,
                                 extra=("--add-depth", str(work / "depth1.tsv")))
                if "Edge weighting: S=2; mean coverage-feature ratio" not in log:
                    raise AssertionError("BAM incremental sample count is wrong")
                check_edges(added, 2, feature_profiles={
                    name: (p[0], p[0], math.ceil(p[0] / 3))
                    for name, p in PROFILES.items()})
                added, log = run("bam_incremental_three", cache=bam_cache,
                                 extra=("--add-depth", str(work / "depth2.tsv")))
                if "Edge weighting: S=3; coverage only" not in log:
                    raise AssertionError("BAM incremental depth did not enter the multi-sample branch")
                check_edges(added, 3, feature_profiles={
                    name: (p[0], p[0], p[1], math.ceil(p[0] / 3))
                    for name, p in PROFILES.items()})

    rows, log = run("default_reuse", extra=("--no_gold",))
    if "configs=5" not in log or "a=" in log:
        raise AssertionError("default reuse must sweep five edge powers without alpha")
    # Verify the selected transform against the actual reported configuration.
    selected = next(line for line in log.splitlines() if "[REUSE] SELECTED" in line)
    power = float(selected.split(" p=", 1)[1].split()[0])
    check_edges(rows, 4, power)
    run("legacy_sweep", overrides={"RABBIT_REUSE_SWEEP": "0.5:1.0"}, expect_ok=False)
    print("Edge checks passed: S>=3 coverage formula and composition invariance, "
          "S<=2 mean coverage ratios, 1/2/3/4 samples, threads, cache, graph reuse, "
          "certification and legacy controls")


if __name__ == "__main__":
    main()
