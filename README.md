# RabbitBin

RabbitBin groups assembled metagenomic contigs into genome bins using sequence
composition and coverage across samples. Start with a contig FASTA and BAM
files mapped to that assembly, or with a precomputed coverage table. RabbitBin
computes coverage internally when given BAM files; external `samtools` and
Python are **not required for normal binning**.

[Install](#installation) · [Prepare inputs](#input-requirements) ·
[Run](#quick-start) · [Outputs](#outputs) · [Troubleshooting](INSTALL.md#troubleshooting)

## Installation

Choose **one** of the following routes. Linux is the tested platform.

### Option A: precompiled binary (no compilation)

For Linux x86-64 with glibc 2.28 or newer, download the binary archive and its
matching `.sha256` file from
[GitHub Releases](https://github.com/RabbitBio/RabbitBin/releases).
For v1.0.0:

```bash
curl -fLO https://github.com/RabbitBio/RabbitBin/releases/download/v1.0.0/rabbitbin-1.0.0-linux-x86_64-glibc2.28.tar.gz
curl -fLO https://github.com/RabbitBio/RabbitBin/releases/download/v1.0.0/rabbitbin-1.0.0-linux-x86_64-glibc2.28.tar.gz.sha256
sha256sum -c rabbitbin-1.0.0-linux-x86_64-glibc2.28.tar.gz.sha256
tar -xzf rabbitbin-1.0.0-linux-x86_64-glibc2.28.tar.gz
cd rabbitbin-1.0.0-linux-x86_64-glibc2.28
./bin/rabbitbin --version
export PATH="$PWD/bin:$PATH"
```

Keep the extracted `bin/` and `lib/` directories together. This route needs no
compiler, Conda, Python, or external `samtools`. An optional offline check is:

```bash
bash test/binary_smoke.sh
```

The `export` applies to the current shell; in later sessions, use the full
path to `bin/rabbitbin` or add that directory to your shell's `PATH`.
See [binary details and limitations](packaging/README.binary.md).

### Option B: build the latest source using Conda

Use this route for the current GitHub version. You need Git and an existing
Conda/Miniforge installation; no administrator access is needed.

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

The environment contains the compiler and libraries needed to build RabbitBin;
it does not request Python or `samtools` for testing. HTSlib is a required
library for reading alignments, and is **not** the external `samtools` command.
Activate `rabbitbin` again in a new shell before running the program.
An optional check using only Bash and the small supplied inputs is:

```bash
bash test/core_smoke.sh "$CONDA_PREFIX/bin/rabbitbin"
```

### Option C: build without Conda, or use Docker

See [INSTALL.md](INSTALL.md) for system packages, Docker, source archives,
and troubleshooting. Full regression tests and their optional dependencies
are documented in [docs/testing.md](docs/testing.md), not required for installation.

## Input requirements

### Contig assembly

Supply assembled contigs in FASTA format (`.fa`, `.fna`, or gzip-compressed
FASTA), not raw FASTQ reads. Contig IDs must be unique. By default, RabbitBin
uses the first whitespace-delimited word of each FASTA header as the ID:
`>contig_1 description` is identified as `contig_1`.

### Coverage: choose BAM files or a depth table

**BAM files (recommended when already available):**

- Map each sample's reads to the **same assembly** supplied to RabbitBin.
  BAM reference names and lengths must match that assembly.
- Supply one **coordinate-sorted** BAM per sample. Paired-end mates belong
  in the same sample BAM, not in two separate sample files.
- A `.bai` index is not required. RabbitBin reads BAMs directly; no external
  `samtools` is needed when the sorted BAMs already exist.
- Mapping and sorting are upstream preparation steps; `rabbitbin bin` does
  not perform them. CRAM input additionally requires `--reference contigs.fa`
  with the matching reference.

For several samples, `bams.txt` is a plain-text file with one path per line:

```text
/data/project/sample1.sorted.bam
/data/project/sample2.sorted.bam
```

Do not add a header, comments, or shell quotation marks to the list. Absolute
paths are easiest to reuse; relative paths are resolved from the directory
where you run RabbitBin, **not** the directory containing `bams.txt`.
Alternatively, pass paths directly with `--bam sample1.bam sample2.bam`.

**Precomputed depth table:**

Use a tab-separated MetaBAT/JGI-style table with a header: three leading
columns (`contigName`, `contigLen`, `totalAvgDepth`), followed by a mean-depth
column and a variance column for each sample. Contig IDs must match the FASTA.
For example, a two-sample table is:

```text
contigName	contigLen	totalAvgDepth	sample1.bam	sample1.bam-var	sample2.bam	sample2.bam-var
contig_1	5000	18	10	2	8	1.5
contig_2	3200	11	6	1	5	0.8
```

Columns are separated by **tabs**, not spaces or commas. RabbitBin can generate
a compatible table without external `samtools`:

```bash
rabbitbin depth --bam-list bams.txt --out depth.tsv --threads 8
```

Use either `--bam`/`--bam-list` **or** `--depth` in one binning command, not both.
Direct BAM input also uses MAPQ-filtered coverage by default, so it is not
identical to supplying an ordinary mean-depth table; see
[coverage modes](docs/usage.md#coverage-modes).

## Quick start

The examples assume `rabbitbin --version` works. Replace the input paths with
your own. `--threads 8` is a resource choice; `--seed 1` fixes the random seed
for reproducibility. Neither changes the default thresholds.

### From one BAM

```bash
mkdir -p results
rabbitbin bin --assembly contigs.fa --bam sample1.sorted.bam \
  --output results/out --threads 8 --seed 1 --bin-fasta
```

### From several BAMs

```bash
mkdir -p results
rabbitbin bin --assembly contigs.fa --bam-list bams.txt \
  --output results/out --threads 8 --seed 1 --bin-fasta
```

### From a precomputed depth table

```bash
mkdir -p results
rabbitbin bin --assembly contigs.fa --depth depth.tsv \
  --output results/out --threads 8 --seed 1 --bin-fasta
```

Default length rules are the same for all three commands:

- Contigs **at least 2,500 bp** form the initial bins.
- Contigs **1,000–2,499 bp** can be recruited into those bins using coverage;
  contigs shorter than 1,000 bp are discarded.
- Only bins with a total length of **at least 200,000 bp** are emitted.

These limits are controlled by `--min-contig`, `--min-small-contig`, and
`--min-bin-size`, respectively. `--bin-fasta` only requests extra output files;
it does not change the binning.

Assembly-only binning is possible by omitting coverage flags, but disables
coverage-based splitting and recruitment. It is **not** the mode used for the
reported coverage-based benchmarks.

## Outputs

`--output results/out` is a **filename prefix**, not a directory of bins.

| Output | Contents |
|--------|----------|
| `results/out.members.tsv` | Contig-to-bin assignments |
| `results/out.bins.tsv` | Per-bin statistics |
| `results/out_bin_001.fa`, etc. | Bin FASTA files, only with `--bin-fasta` |
| `results/out.unbinned.fa` | Unbinned contigs, only with `--unbinned` |

Without `--bin-fasta`, the membership and statistics tables are still written.
If no bins are emitted for a tiny input, check the 200 kb minimum bin size;
the smoke tests lower this limit **only for their tiny test fixture**.

For a reproducible record, save the full command, `rabbitbin --version`, input
files, and console output (for example, append `> run.log 2>&1`). The default
seed is time-based (`0`), so specify a nonzero `--seed` for reported runs.

## Further documentation

- [Installation, dependencies, and troubleshooting](INSTALL.md)
- [Optional tests](docs/testing.md)
- [Additional commands and options](docs/usage.md): depth, AMBER evaluation,
  markers, and the convenience wrapper
- [Default pipeline and PMH validation](docs/pipeline.md)
- [Performance diagnostics and optimization checks](docs/performance.md)

Run `rabbitbin bin --help` or `rabbitbin <command> --help` for all options.
Marker-based QC/purification needs extra tools and marker data only when
explicitly requested; default binning does not need them.
RabbitBin uses [RabbitBAM](https://github.com/RabbitBio/RabbitBAM/tree/sortedbam)
for parallel BAM I/O.

## License

RabbitBin is distributed under the [BSD 3-Clause license](LICENSE). Attribution
and third-party notices are recorded in [license.txt](license.txt).
