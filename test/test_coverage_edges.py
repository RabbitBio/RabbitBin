"""Check S>=3 coverage edge weights and the original S<=2 fallback."""

import csv
import math
import os
from pathlib import Path
import random
import subprocess
import sys


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


def check_edges(rows, samples, power=1.0, cutoff=DEFAULT_EDGE_CUTOFF):
    means = [sum(p[s] for p in PROFILES.values()) / len(PROFILES)
             for s in range(samples)]
    if not rows:
        raise AssertionError("fixture produced no candidate edges")
    for key, row in rows.items():
        if samples == 1:
            near(float(row["weight"]), float(row["sComp"]) ** power,
                 f"{key} single-sample composition fallback")
            continue
        a, b = (PROFILES[row[n]][:samples] for n in ("name_i", "name_j"))
        rho = spearman(a, b)
        denominator = sum(max(x, y) / m for x, y, m in zip(a, b, means))
        jcov = (sum(min(x, y) / m for x, y, m in zip(a, b, means)) /
                denominator) if denominator else 0.0
        raw = min(max(rho, 0.0), jcov)
        if samples >= 3:
            weight = raw ** power if raw >= cutoff else 0.0
        else:
            # The original two-sample path uses composition after its signed
            # depth gate, then applies the same edge threshold and power.
            composition = float(row["sComp"]) if rho >= -0.3 else 0.0
            weight = composition ** power if composition >= cutoff else 0.0
        near(float(row["rho"]), rho, f"{key} rho")
        near(float(row["Jcov"]), jcov, f"{key} Jcov")
        near(float(row["dterm"]), raw, f"{key} raw weight")
        near(float(row["weight"]), weight, f"{key} scored weight")

    if samples <= 2:
        scaled = rows[tuple(sorted(("anchor", "scaled")))]
        near(float(scaled["weight"]), float(scaled["sComp"]) ** power,
             "low-sample composition fallback must retain magnitude-mismatched pair")
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
            cache=None, expect_ok=True, threads=1):
        prefix = work / label
        dump = work / f"{label}.pairs.tsv"
        cmd = [str(binary), "bin", "--output", str(prefix), "--seed", "42",
               "--threads", str(threads), "--min-bin-size", "0", "--no-recruit",
               "--no-split", "--no-bin-fasta"]
        if cache:
            cmd += ["--load-cache", str(cache)]
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
            if "alpha:power mixing was removed" not in completed.stdout:
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

    rows, log = run("default_reuse", extra=("--no_gold",))
    if "configs=5" not in log or "a=" in log:
        raise AssertionError("default reuse must sweep five edge powers without alpha")
    # Verify the selected transform against the actual reported configuration.
    selected = next(line for line in log.splitlines() if "[REUSE] SELECTED" in line)
    power = float(selected.split(" p=", 1)[1].split()[0])
    check_edges(rows, 4, power)
    run("legacy_sweep", overrides={"RABBIT_REUSE_SWEEP": "0.5:1.0"}, expect_ok=False)
    print("Edge checks passed: S>=3 coverage formula and composition invariance, "
          "S<=2 fallback, 1/2/3/4 samples, threads, cache, graph reuse, "
          "certification and legacy controls")


if __name__ == "__main__":
    main()
