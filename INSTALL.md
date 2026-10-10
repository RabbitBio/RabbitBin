# Installing RabbitBin

Linux is the tested platform. Choose a precompiled binary, a Conda-assisted
source build, a system-compiler build, or Docker. You do not need all four.

## What is actually required?

| Purpose | Dependencies |
|---------|--------------|
| Run the precompiled release | Linux x86-64, glibc 2.28+; non-glibc runtime libraries are bundled |
| Compile from source | CMake 3.16+, C++17 compiler with OpenMP (GCC 7+ or recent Clang), Make, Git, pkg-config; Boost 1.66+, zlib 1.2.11+, HTSlib 1.13+, libdeflate development libraries |
| Build missing libraries automatically | Also Autoconf, Automake, Libtool, and network access |
| Run the optional full regression suite | Also Python 3.6+ and external `samtools`; see [testing](docs/testing.md) |
| Generate a marker map for optional QC/purification | Prodigal, HMMER (`hmmsearch`), and a single-copy-marker HMM collection |

Normal BAM/depth binning requires **neither Python nor external `samtools`**.
HTSlib is a C library linked into RabbitBin, not the `samtools` executable.
Optional Python/Perl helper scripts need their respective interpreters only
when invoked. Marker tools and databases are not downloaded or used by the
normal binning pipeline.

## Precompiled Linux release

Follow the [README download commands](README.md#option-a-precompiled-binary-no-compilation).
Keep `bin/` and `lib/` together, verify the accompanying SHA-256 file, then run
`bin/rabbitbin --version`. The package includes tiny test inputs and
`bash test/binary_smoke.sh` for an optional offline check requiring only
Bash, ordinary Unix utilities, and gzip.

See [binary instructions](packaging/README.binary.md) and
[known limitations](packaging/RELEASE_NOTES.md). Use a source build for a
different platform or for changes newer than the release tag.

## Conda-assisted source build

Run these commands with Git and Conda/Miniforge available:

```bash
git clone https://github.com/RabbitBio/RabbitBin.git
cd RabbitBin
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate rabbitbin
cmake -S . -B build-conda -DCMAKE_BUILD_TYPE=Release -DNO_TESTING=ON \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
cmake --build build-conda --parallel 4
cmake --install build-conda
rabbitbin --version
```

`environment.yml` uses `conda-forge`, `bioconda`, and `nodefaults`. It requests
source-build dependencies, not optional regression tools. Package managers
may install their own transitive dependencies; these are distinct from
RabbitBin's runtime requirements. Activate this environment in each new
shell. `NO_TESTING=ON` skips test targets only; it does not change the algorithm
or disable BAM input.

An optional small check, from the source directory:

```bash
bash test/core_smoke.sh "$CONDA_PREFIX/bin/rabbitbin"
```

## System-compiler source build (no Conda)

On Ubuntu/Debian, install the needed development packages:

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends build-essential cmake git pkg-config \
  libboost-program-options-dev libboost-filesystem-dev libboost-system-dev \
  libboost-graph-dev libboost-serialization-dev libboost-iostreams-dev \
  libboost-regex-dev zlib1g-dev libhts-dev libdeflate-dev
```

This avoids the all-components Boost package and does not request Python or
`samtools`. On Rocky/RHEL, one route is to install the compiler, Boost and
Autotools, and allow CMake to build pinned zlib/HTSlib/libdeflate copies:

```bash
sudo dnf install -y gcc gcc-c++ cmake git make autoconf automake libtool \
  pkgconf-pkg-config boost-devel
```

Repository availability depends on the distribution. If packages are
unavailable or you lack administrator access, use the binary or Conda route.
Then build and install under your own account:

```bash
git clone https://github.com/RabbitBio/RabbitBin.git
cd RabbitBin
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNO_TESTING=ON \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build --parallel 4
cmake --install build
export PATH="$HOME/.local/bin:$PATH"
rabbitbin --version
# Optional; only Bash, Unix utilities, gzip, and the supplied fixtures:
bash test/core_smoke.sh "$HOME/.local/bin/rabbitbin"
```

The uninstalled executable is `build/src/rabbitbin`. Make the PATH addition
persistent in your shell configuration, or call the full installed path.

RabbitBin first uses installed zlib, HTSlib and libdeflate development
packages. Missing copies are downloaded at pinned revisions and built
locally, not installed system-wide. An offline source build must provide all
three development packages in advance. The fallback needs Autotools; on
Debian/Ubuntu these are `autoconf automake libtool`.

## Optional Docker image

From the checkout:

```bash
docker build \
  --build-arg RABBITBIN_SOURCE_REVISION="$(git rev-parse --short=12 HEAD)" \
  -t rabbitbin:local .
docker run --rm rabbitbin:local rabbitbin --version
mkdir -p results
docker run --rm -v "$PWD:/work" -w /work rabbitbin:local \
  rabbitbin bin --assembly contigs.fa --bam sample1.sorted.bam \
  --output results/out --threads 8 --seed 1 --bin-fasta
```

Input paths must exist **inside the container**. With this mount, use paths
under `/work` (including paths inside a BAM list). The runtime image does not
explicitly install Python, Perl, or external `samtools`; optional helper
scripts needing those tools are outside the default image workflow.

## Troubleshooting

- **`rabbitbin: command not found`:** activate the Conda environment, add the
  installation's `bin/` directory to PATH, or use the full executable path.
- **Missing HTSlib/Boost or a loader error:** use a fresh build directory and
  the correct `CMAKE_PREFIX_PATH`. Installed binaries preserve non-system
  library locations; do not move/remove those dependencies after installing.
  A precompiled release must keep its sibling `lib/` directory.
- **No bin FASTA files:** add `--bin-fasta`; the default output is TSV tables.
- **No bins for a small assembly:** the default minimum output bin size is
  200,000 bp. A tiny installation fixture needs a smaller test-only limit.
- **BAM/depth contigs do not match:** check that the exact same assembly was
  used for mapping and that IDs match its first whitespace-delimited headers.
- **Full tests are skipped:** the documented normal build uses `NO_TESTING=ON`;
  follow [docs/testing.md](docs/testing.md) if you want the regression suite.

Use a **fresh build directory** after changing dependency environments or
pulling code that changes the build configuration. CMake caches library
locations. `CMAKE_PREFIX_PATH` also supports non-Conda prefixes;
`HTSLIB_ROOT` can locate an HTSlib-only installation.

## Advanced build options and provenance

- `RABBITBIN_NATIVE_ARCH=OFF` (default): portable CPU baseline.
- `RABBITBIN_NATIVE_ARCH=ON`: tune for the build host; use only on compatible CPUs.
- `RABBITBIN_USE_SYSTEM_ZLIB=OFF`, `RABBITBIN_USE_SYSTEM_HTSLIB=OFF`, or
  `RABBITBIN_USE_SYSTEM_LIBDEFLATE=OFF`: force a pinned dependency build.
- `RABBITBIN_STATIC_BOOST=ON`: use static Boost libraries, if available.
- `NO_TESTING=ON`: omit test targets; `OFF` enables test discovery.

Git source archives have no `.git` directory. `git archive` substitutes the
commit into `SOURCE_REVISION` via `.gitattributes`; CMake uses that value when
Git checkout metadata is absent. Build an extracted archive with the same
commands as a checkout. Plain source copies with no revision metadata report
`commit unknown` and still build. A packager can explicitly supply a known
revision with `-DRABBITBIN_SOURCE_REVISION=<commit>`; do not label modified
sources as an unmodified release.

The version comes from `VERSION`; `rabbitbin --version` also reports the
source commit. For release packaging, see [packaging](packaging/README.binary.md)
and `.github/workflows/release.yml`. RabbitBin uses the [BSD 3-Clause license](LICENSE);
inherited and third-party notices are in [license.txt](license.txt).
