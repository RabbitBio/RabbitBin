# RabbitBin 1.0.0 — portable Linux reviewer build

The binary archive targets Linux x86-64 / glibc 2.28+, includes non-glibc runtime
libraries, and requires no compilation or Conda environment. Download both the
archive and its SHA-256 file; extraction and self-test instructions are inside.
The exact source commit is embedded in `rabbitbin --version` and `BUILDINFO.txt`.

Installation fixes include explicit Conda channels and an environment file,
GCC 15 header compatibility, dependency-prefix discovery, Git source-archive
revision metadata, CTest fixture dependencies, and shell quoting/error handling.
No binning algorithms, default thresholds or dataset-specific tuning changed.

Release CI builds on Rocky Linux 8 with portable CPU settings, runs the full
CTest suite, and gates publication on unpacked-binary tests in fresh Rocky
Linux 8, Ubuntu 22.04 and Ubuntu 24.04 containers. The binary is relocated before
testing; BAM input, gzip FASTA input and space-containing paths are exercised.
These small installation tests do not replace CAMI2 accuracy/performance tests.

Known limitations:

- The inherited, non-default `rabbit_depth --unmappedFastq` / `BADMAP=1`
  export path can crash. Leave `BADMAP=0`, the default. Standard BAM binning
  and the default wrapper are tested.
- An existing RabbitBAM atomic-memory-order compiler warning remains under
  separate review; no concurrency or BAM-processing algorithms were changed.
- Optional marker workflows require their documented external tools/models.
- No native macOS, Windows, ARM or musl binary is provided; use the source build
  on other supported environments. Full map/refine/marker workflows and
  cross-CPU performance are not comprehensively validated by these smoke tests.

The source archive, binary archive and checksum are versioned together.
