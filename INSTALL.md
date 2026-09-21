# Building RabbitBin

Conda and Docker are optional. A normal system compiler build is the primary
installation route.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler with OpenMP (GCC 7+ or recent Clang)
- Boost 1.66 or newer (`program_options`, `filesystem`, `system`, `graph`,
  `serialization`, `iostreams`, and `regex`)
- Git, Make, Autoconf, Automake, Libtool, and pkg-config
- zlib 1.2.11+, HTSlib 1.13+, and libdeflate

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
ctest --test-dir build --output-on-failure
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
conda create -n rabbitbin -c conda-forge \
  cmake make compilers boost-cpp zlib htslib libdeflate
conda activate rabbitbin
cmake -S . -B build-conda -DCMAKE_BUILD_TYPE=Release
cmake --build build-conda --parallel
ctest --test-dir build-conda --output-on-failure
```

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
