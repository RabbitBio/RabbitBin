# RabbitBin precompiled Linux package

This package is for **Linux x86-64 with glibc 2.28 or newer** (not Alpine/musl,
macOS, Windows, or ARM). It uses the portable CPU build, not `-march=native`.
The OS supplies glibc and the dynamic loader; the other required shared
libraries are included in `lib/`. Keep `bin/` and `lib/` together when moving
the package. No root access, Conda, compiler, external samtools or Python is
needed for the main BAM-binning executable.

## Quick start

Download the `.tar.gz` and matching `.sha256` files from the same release:

```bash
sha256sum -c rabbitbin-1.0.0-linux-x86_64-glibc2.28.tar.gz.sha256
tar -xzf rabbitbin-1.0.0-linux-x86_64-glibc2.28.tar.gz
cd rabbitbin-1.0.0-linux-x86_64-glibc2.28
./bin/rabbitbin --version
./bin/rabbitbin bin --assembly contigs.fa --bam sample.bam --output result --threads 16
```

Inputs remain the user's assembly and corresponding sorted BAM(s). No datasets
or pre-trained models are downloaded. Use `./bin/rabbitbin bin --help` for all
options. The release does not change binning defaults or scientific thresholds.

## Self-test

The included tiny FASTA/BAM fixture lets a reviewer test the extracted package:

```bash
bash test/binary_smoke.sh
# More checks, optionally, when Python 3.6+ is installed:
python3 test/test_installation.py --prefix . --source .
```

These tests use `--min-bin-size 0` only to accommodate the 5,404 bp fixture;
they are installation tests, not CAMI2 benchmarks. Shell self-test output is
left in the temporary directory printed at completion.

## Provenance, optional features and licenses

`BUILDINFO.txt` records the full source commit, dependency revisions, compiler,
build time and ABI requirement. The matching source archive and build recipe
are provided with the release. See `LICENSE`, `license.txt` and `licenses/` for
project and dependency notices. GCC runtime libraries are distributed under
their accompanying licenses and GCC Runtime Library Exception; source is
available from https://gcc.gnu.org/releases.html and the matching Rocky Linux
source RPMs at https://dl.rockylinux.org/vault/rocky/ . Boost source is available
from https://www.boost.org/releases/ . Pinned zlib, HTSlib and libdeflate sources
and revisions are recorded in the repository's `cmake/` build recipes.

Optional marker workflows still need Prodigal, HMMER and user-supplied HMMs;
auxiliary Python/Perl scripts need their respective interpreters. This archive
does not bundle those optional tools. The inherited non-default
`rabbit_depth --unmappedFastq` / wrapper `BADMAP=1` path has a known crash and
is not supported by this release; leave `BADMAP=0` (the default). See
`RELEASE_NOTES.md` for validation scope and other known limitations.
