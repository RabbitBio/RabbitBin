"""Smoke-test an installed prefix without depending on the checkout cwd."""
import argparse
import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--dependency-prefix", type=Path)
    args = parser.parse_args()
    if bool(args.build_dir) != bool(args.dependency_prefix):
        parser.error("--build-dir and --dependency-prefix must be used together")
    if args.build_dir:
        cache = {}
        for line in (args.build_dir / "CMakeCache.txt").read_text().splitlines():
            if line.startswith(("#", "//")) or "=" not in line or ":" not in line:
                continue
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
        keys = ["Boost_INCLUDE_DIR", "HTSlib_INCLUDE_DIR", "HTSlib_LIBRARY",
                "LIBDEFLATE_INCLUDE_DIR", "LIBDEFLATE_LIBRARY", "ZLIB_INCLUDE_DIR",
                "ZLIB_LIBRARY_RELEASE"]
        keys += [k for k in cache if k.startswith("Boost_") and k.endswith("_LIBRARY_RELEASE")]
        for key in keys:
            if not cache.get(key) or cache[key].endswith("-NOTFOUND"):
                raise AssertionError("missing selected dependency: " + key)
            path = Path(cache[key]).resolve()
            try:
                path.relative_to(args.dependency_prefix.resolve())
            except ValueError:
                raise AssertionError("dependency outside requested prefix: " + str(path))
    bindir, data = args.prefix.resolve() / "bin", args.source.resolve() / "test"
    binary = bindir / "rabbitbin"
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("RABBIT_", "RB_", "OMP_"))
           and k not in ("LD_LIBRARY_PATH", "LD_PRELOAD", "DYLD_LIBRARY_PATH",
                         "BADMAP", "PCTID", "MINDEPTH")}

    with tempfile.TemporaryDirectory(prefix="rabbitbin installed ") as tmp:
        work = Path(tmp)

        def run(command):
            result = subprocess.run([str(x) for x in command], cwd=str(work), env=env,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    universal_newlines=True, timeout=60)
            if result.returncode:
                raise AssertionError(result.stdout)
            return result.stdout

        version = run([binary, "--version"]).strip()
        if not version.startswith("RabbitBin ") or "(commit " not in version:
            raise AssertionError("missing version metadata: " + version)
        run([binary, "--help"])
        # Separate executables are part of the public installation too.
        run([bindir / "rabbit_depth"])
        for filename in ("rabbit_overlap", "run_rabbitbin.sh"):
            if not os.access(str(bindir / filename), os.X_OK):
                raise AssertionError("missing installed executable: " + filename)

        assembly, bam = work / "assembly space.fa", work / "sample one.bam"
        shutil.copyfile(str(data / "contigs.fa"), str(assembly))
        shutil.copyfile(str(data / "contigs-1000.fastq.bam"), str(bam))
        compressed = work / "assembly space.fa.gz"
        with assembly.open("rb") as source, gzip.open(str(compressed), "wb") as dest:
            shutil.copyfileobj(source, dest)
        common = ["--threads", "2", "--min-bin-size", "0", "--seed", "42"]
        run([binary, "bin", "--assembly", assembly, "--bam", bam,
             "--output", "bam-out"] + common)
        run([binary, "bin", "--assembly", compressed,
             "--depth", data / "contigs_depth.txt", "--output", "gz-out"] + common)
        env["RB_LABEL"] = "installed test"
        run([bindir / "run_rabbitbin.sh"] + common + [assembly, bam])
        outputs = [work / "bam-out.members.tsv", work / "gz-out.members.tsv"]
        outputs += list(work.glob("*.rabbitbin-*/out.members.tsv"))
        if len(outputs) != 3 or any(len(p.read_text().splitlines()) != 2 for p in outputs):
            raise AssertionError("missing or empty smoke-test binning output")
        # TotalDepth is not comparable across these input routes: fused BAM
        # adds coverage features, while the stored TSV has one coverage column.
        # Check the actual membership and nonzero coverage, not byte identity.
        sequence = assembly.read_text().splitlines()[0][1:].split()[0]
        for path in outputs:
            with path.open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            if {(r["BinNum"], r["SequenceName"]) for r in rows} != {("1", sequence)}:
                raise AssertionError("incorrect smoke-test members: " + str(path))
            if any(not float(r["TotalDepth"]) > 0 for r in rows):
                raise AssertionError("missing coverage: " + str(path))
        print(version)
        print("installed BAM, gzip and wrapper smoke tests passed")


if __name__ == "__main__":
    main()
