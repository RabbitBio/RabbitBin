# Building RabbitBin

## Dependencies

- CMake ≥ 3.5
- C++17 compiler with OpenMP
- Boost (program_options, filesystem, system, graph, serialization, iostreams)
- zlib
- htslib (only for `rabbit_depth` and `rabbit_overlap`)

## Build

```bash
git clone <your-rabbitbin-repo>
cd RabbitBin
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
make rabbitbin -j$(nproc)
make install   # optional
```

The main binary is `build/src/rabbitbin`.

## Version

Release version is read from the top-level `VERSION` file (currently `1.0.0`). No git metadata is embedded in builds.

## License

See `license.txt` for licensing terms.
