# Building RabbitBin

Conda and Docker are optional. A normal system compiler build is the primary
installation route.

## Precompiled Linux release

Reviewers can use the Linux x86-64 / glibc 2.28+ archive from
[GitHub Releases](https://github.com/RabbitBio/RabbitBin/releases). It bundles
non-glibc runtime libraries and does not require compilation or Conda. Keep
`bin/` and `lib/` together, verify the accompanying SHA-256 file, then run
`bin/rabbitbin --version`. The package includes tiny test inputs and
`bash test/binary_smoke.sh` for an offline installation check.

See [binary instructions](packaging/README.binary.md) and
[known limitations](packaging/RELEASE_NOTES.md). Source builds remain available
for other environments. Maintainers can reproduce packaging with
`bash scripts/package_linux.sh BUILD_DIR OUTPUT_DIR` after a clean, tested
Rocky Linux 8 build with all three bundled dependency options enabled; the
release workflow records the complete recipe and tests fresh runtime images.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler with OpenMP (GCC 7+ or recent Clang)
- Boost 1.66 or newer (`program_options`, `filesystem`, `system`, `graph`,
  `serialization`, `iostreams`, and `regex`)
- Git, Make, Autoconf, Automake, Libtool, and pkg-config
- zlib 1.2.11+, HTSlib 1.13+, and libdeflate
- Python 3.6+ and samtools for all regression tests; the core executable does
  not require either interpreter or external samtools at runtime

RabbitBin first uses installed zlib, HTSlib and libdeflate development
packages. Missing copies are downloaded at pinned revisions and built locally;
they are not installed system-wide. An offline build must provide all three
development packages in advance.

## Source build

```bash
git clone https://github.com/RabbitBio/RabbitBin.git
cd RabbitBin
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build --parallel
(cd build && ctest --output-on-failure)
cmake --install build
"$HOME/.local/bin/rabbitbin" --version
```

The uninstalled main binary is `build/src/rabbitbin`.

Installed binaries retain runtime search paths for dependencies supplied from
a non-system prefix. If dependencies are later moved or removed, rebuild or
make their `lib` directory available to the dynamic loader.

## Build options

- `RABBITBIN_NATIVE_ARCH=OFF` (default): portable CPU baseline, suitable for
  releases, containers, clusters, and reviewer machines.
- `RABBITBIN_NATIVE_ARCH=ON`: tune for the build host with `-march=native`;
  use only when the binary stays on compatible hardware.
- `RABBITBIN_USE_SYSTEM_ZLIB=OFF`, `RABBITBIN_USE_SYSTEM_HTSLIB=OFF`, or
  `RABBITBIN_USE_SYSTEM_LIBDEFLATE=OFF`: force the pinned dependency build.
- `RABBITBIN_STATIC_BOOST=ON`: link the Boost components statically when static
  Boost libraries are available.

## Optional Conda environment

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate rabbitbin
cmake -S . -B build-conda -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
cmake --build build-conda --parallel
(cd build-conda && ctest --output-on-failure)
cmake --install build-conda
rabbitbin --version
```

`environment.yml` includes `conda-forge`, `bioconda` and `nodefaults`; HTSlib
and samtools come from bioconda. It installs C/C++ compilers and all test
dependencies. Always configure a fresh build directory after changing the
dependency environment: CMake caches previously selected library locations.
`CMAKE_PREFIX_PATH` is also supported for non-Conda dependency prefixes;
`HTSLIB_ROOT` can locate an HTSlib-only installation.

## Source archives and version provenance

Git source archives have no `.git` directory. `git archive` substitutes the
commit into `SOURCE_REVISION` via `.gitattributes`, and CMake uses it when Git
checkout metadata is absent. Build and test the extracted directory with the
same commands as above; cloning is not required for this route.

Plain source copies with no revision metadata report `commit unknown` and
still build and test. A release packager can explicitly supply a known source
revision with `-DRABBITBIN_SOURCE_REVISION=<commit>`. Do not label modified
sources with the revision of an unmodified release.

## Installation smoke test

After installing, check the installed executables from outside the source
directory. From the checkout, for a per-user install:

```bash
python3 test/test_installation.py --prefix "$HOME/.local" --source .
```

For the Conda install, use `--prefix "$CONDA_PREFIX"`. The test creates a
temporary working directory, clears RabbitBin tuning and dynamic-loader
overrides, and checks version/help, BAM binning, compressed FASTA + depth,
and the convenience wrapper with space-containing paths. It uses only the
repository's small fixtures and does not run CAMI2 benchmarks. Python tests
are skipped with a configure message when their prerequisites are absent;
install Python and samtools to run the complete suite.

## Optional Docker image

```bash
docker build \
  --build-arg RABBITBIN_SOURCE_REVISION="$(git rev-parse --short=12 HEAD)" \
  -t rabbitbin:local .
docker run --rm rabbitbin:local rabbitbin --version
```

## Version and license

The release version comes from `VERSION`; Git checkout builds additionally
embed the source commit shown by `rabbitbin --version`. RabbitBin uses the BSD
3-Clause license in `LICENSE`; inherited and third-party notices are in
`license.txt`.
