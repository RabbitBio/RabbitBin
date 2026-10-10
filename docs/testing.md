# Optional installation checks and regression tests

Tests are not needed to run RabbitBin. The normal installation commands use
`-DNO_TESTING=ON` and do not request Python or external `samtools`.

## Small installation check (no Python or samtools)

For an extracted binary release, from its top-level directory:

```bash
bash test/binary_smoke.sh
```

For a source build, from the checkout, pass the installed executable's path:

```bash
# Conda install:
bash test/core_smoke.sh "$CONDA_PREFIX/bin/rabbitbin"
# Or a per-user system install:
bash test/core_smoke.sh "$HOME/.local/bin/rabbitbin"
```

The source-build check uses Bash, ordinary Unix utilities, awk, and gzip. It
tests a single BAM, a BAM list, gzip FASTA with stored depth, and a newly
generated depth table, including paths containing spaces. It checks that the
expected contig is binned, coverage is positive, and bin FASTA files exist.
The script clears RabbitBin tuning/loader overrides, runs in a temporary
directory outside the source tree, and prints the retained output directory.
The `--min-bin-size 0` setting is **only for the tiny supplied fixture**; it is
not a recommended setting for real datasets. No CAMI2 experiment is launched.

## Full regression suite (for development)

Only use this section if you want to test or modify the implementation.
Python 3.6+ drives several regression checks; `samtools` prepares reference
alignments for BAM/depth equivalence tests. Neither is called by normal
RabbitBin binning. In the existing Conda source-build environment:

```bash
conda activate rabbitbin
conda install --yes --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda 'python>=3.6' 'samtools>=1.13'
cmake -S . -B build-tests -DCMAKE_BUILD_TYPE=Release -DNO_TESTING=OFF \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
cmake --build build-tests --parallel 4
(cd build-tests && ctest --output-on-failure)
cmake --install build-tests
python3 test/test_installation.py --prefix "$CONDA_PREFIX" --source . \
  --build-dir build-tests --dependency-prefix "$CONDA_PREFIX"
```

Use a fresh test-build directory after changing dependencies. For a system
build, install optional `python3` and `samtools` through your package manager,
omit `CMAKE_PREFIX_PATH`, choose your installation prefix, and omit the last
command's `--build-dir`/`--dependency-prefix` pair (it specifically checks a
single dependency prefix).

If prerequisites are missing, CMake reports which conditional tests it skips.
That reduced suite is not the full regression suite. CI first checks a minimal
Conda build, then installs the extra tools explicitly for the full tests.
Release-packaging jobs likewise use test tools during validation; those tools
are not dependencies of the distributed binary.
