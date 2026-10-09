#!/usr/bin/env bash
# Package a tested, bundled-dependency build. No root access is needed here.
set -euo pipefail
if [ "$#" -ne 2 ]; then
  echo "Usage: bash scripts/package_linux.sh BUILD_DIR OUTPUT_DIR" >&2
  exit 2
fi
source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=$(cd -- "$1" && pwd)
mkdir -p -- "$2"
output_dir=$(cd -- "$2" && pwd)
patchelf=${PATCHELF:-patchelf}
version=$(tr -d '\r\n' < "$source_dir/VERSION")
revision=$(git -C "$source_dir" rev-parse HEAD)
if [ -n "$(git -C "$source_dir" status --porcelain)" ]; then
  echo "Refusing to package a dirty checkout" >&2
  exit 1
fi
test "$(uname -m)" = x86_64
grep -q '^RABBITBIN_NATIVE_ARCH:BOOL=OFF$' "$build_dir/CMakeCache.txt"
for dep in ZLIB HTSLIB LIBDEFLATE; do
  grep -q "^RABBITBIN_USE_SYSTEM_${dep}:BOOL=OFF$" "$build_dir/CMakeCache.txt"
done
name="rabbitbin-${version}-linux-x86_64-glibc2.28"
stage_root=$(mktemp -d "$output_dir/.package.XXXXXXXX")
package="$stage_root/$name"
cmake --install "$build_dir" --prefix "$package"
test "$("$package/bin/rabbitbin" --version)" = "RabbitBin $version (commit ${revision:0:12})"
mkdir -p "$package/lib" "$package/licenses" "$package/test"

# Keep glibc and the ELF loader supplied by the operating system. Bundle all
# other dependencies, including transitive Boost/C++/OpenMP runtime libraries.
is_system_lib() {
  case "$1" in
    libc.so.6|libm.so.6|libpthread.so.0|libdl.so.2|librt.so.1|libresolv.so.2|libutil.so.1|ld-linux*) return 0 ;;
    *) return 1 ;;
  esac
}
declare -A dependency_paths=()
for exe in rabbitbin rabbit_depth rabbit_overlap; do
  dependencies=$(ldd "$package/bin/$exe")
  if [[ "$dependencies" == *'not found'* ]]; then
    echo "$dependencies" >&2
    exit 1
  fi
  while read -r soname path; do
    if ! is_system_lib "$soname"; then
      cp -L -- "$path" "$package/lib/$soname"
      dependency_paths["$soname"]=$path
    fi
  done < <(printf '%s\n' "$dependencies" | awk '$2 == "=>" && $3 ~ /^\// {print $1, $3}')
done
for library in "$package"/lib/*; do
  "$patchelf" --set-rpath '$ORIGIN' "$library"
done
for exe in rabbitbin rabbit_depth rabbit_overlap; do
  # Older binutils can corrupt relocated .dynstr sections if strip is run
  # after patchelf. Strip the original ELF first; patching must be the last edit.
  strip --strip-unneeded "$package/bin/$exe"
  "$patchelf" --set-rpath '$ORIGIN/../lib' "$package/bin/$exe"
done

# A newer build host or accidentally selected runtime must not silently raise
# the advertised ABI requirement.
max_glibc=$(readelf --version-info "$package"/bin/rabbitbin "$package"/bin/rabbit_depth \
  "$package"/bin/rabbit_overlap "$package"/lib/* | \
  grep -oE 'GLIBC_[0-9]+\.[0-9]+(\.[0-9]+)?' | sed 's/GLIBC_//' | sort -Vu | tail -n 1)
test "$(printf '%s\n' "$max_glibc" 2.28 | sort -V | tail -n 1)" = 2.28

cp "$source_dir/LICENSE" "$source_dir/license.txt" "$package/"
cp "$source_dir/packaging/README.binary.md" "$package/README.md"
cp "$source_dir/packaging/RELEASE_NOTES.md" "$package/RELEASE_NOTES.md"
while IFS= read -r notice; do
  mkdir -p -- "$package/licenses/project/$(dirname -- "$notice")"
  cp -- "$source_dir/$notice" "$package/licenses/project/$notice"
done < <(git -C "$source_dir" ls-files '*LICENSE*' '*COPYING*' '*license*')
# Several vendored extensions carry their notices inside source headers, not
# separate LICENSE files. Preserve those originals as well.
cp -R "$source_dir/src/align/strobe/ext" "$package/licenses/vendored-strobe-extensions"
# boost-devel may install individual component RPMs without the boost meta-RPM.
boost_notice=
for candidate in /usr/share/licenses/boost*/LICENSE_1_0.txt; do
  if [ -f "$candidate" ]; then
    boost_notice=$candidate
    break
  fi
done
test -n "$boost_notice"
cp "$boost_notice" "$package/licenses/Boost.txt"
mkdir -p "$package/licenses/GCC"
cp /usr/share/licenses/gcc/COPYING* "$package/licenses/GCC/"
cp "$build_dir/contrib/zlib-prefix/src/zlib/LICENSE" "$package/licenses/zlib.txt"
cp "$build_dir/contrib/libdeflate-prefix/src/libdeflate_external/COPYING" "$package/licenses/libdeflate.txt"
cp "$build_dir/contrib/htslib-prefix/src/htslib/LICENSE" "$package/licenses/HTSlib.txt"
cp "$build_dir/contrib/htslib-prefix/src/htslib/htscodecs/LICENSE.md" "$package/licenses/htscodecs.md"
# Some platform packages (e.g. optional NUMA support) have additional notices.
if command -v rpm >/dev/null; then
  for library in "$package"/lib/*; do
    soname=${library##*/}
    system_path=${dependency_paths[$soname]}
    if [ -n "$system_path" ]; then
      rpm_name=$(rpm -qf "$system_path" 2>/dev/null || true)
      if rpm -q "$rpm_name" >/dev/null 2>&1; then
        while IFS= read -r notice; do
          if [ -f "$notice" ]; then
            mkdir -p "$package/licenses/rpm/$rpm_name"
            cp -- "$notice" "$package/licenses/rpm/$rpm_name/"
          fi
        done < <(rpm -ql "$rpm_name" | grep -E '/licenses/|/(COPYING|LICENSE)' || true)
      fi
    fi
  done
fi
cp "$source_dir/test/contigs.fa" "$source_dir/test/contigs-1000.fastq.bam" \
   "$source_dir/test/contigs_depth.txt" "$source_dir/test/test_installation.py" \
   "$source_dir/test/binary_smoke.sh" "$package/test/"
{
  printf 'version=%s\ncommit=%s\nbuild_utc=%s\n' "$version" "$revision" "$(date -u +%FT%TZ)"
  printf 'platform=linux-x86_64\nminimum_glibc=2.28\nobserved_glibc_symbols=%s\n' "$max_glibc"
  printf 'native_arch=OFF\n'
  gcc --version | head -n 1
  cmake --version | head -n 1
  "$patchelf" --version
  strip --version | head -n 1
  for dep in zlib-prefix/src/zlib htslib-prefix/src/htslib libdeflate-prefix/src/libdeflate_external; do
    printf '%s=' "$dep"
    git -C "$build_dir/contrib/$dep" rev-parse HEAD
  done
  if command -v rpm >/dev/null; then
    rpm -q gcc libgcc libstdc++ libgomp boost boost-devel || true
  fi
} > "$package/BUILDINFO.txt"

# Verify relocation, not just execution from the original installation prefix.
relocated="$stage_root/relocated package with spaces"
mv -- "$package" "$relocated"
env -u LD_LIBRARY_PATH -u LD_PRELOAD bash "$relocated/test/binary_smoke.sh"
env -u LD_LIBRARY_PATH -u LD_PRELOAD python3 "$relocated/test/test_installation.py" \
  --prefix "$relocated" --source "$relocated"
for exe in rabbitbin rabbit_depth rabbit_overlap; do
  linkage=$(env -u LD_LIBRARY_PATH -u LD_PRELOAD ldd "$relocated/bin/$exe")
  test "${linkage#*not found}" = "$linkage"
  while read -r soname path; do
    if ! is_system_lib "$soname" && [[ "$path" != "$relocated/lib/"* && "$path" != "$relocated/bin/../lib/"* ]]; then
      echo "Dependency escaped package: $soname => $path" >&2
      exit 1
    fi
  done < <(printf '%s\n' "$linkage" | sed -nE 's/^[[:space:]]*([^ ]+) => (.*) \(0x[0-9a-f]+\)$/\1 \2/p')
done
mv -- "$relocated" "$package"
archive="$output_dir/$name.tar.gz"
test ! -e "$archive"
tar -C "$stage_root" -czf "$archive" "$name"
(cd "$output_dir" && sha256sum "$name.tar.gz" > "$name.tar.gz.sha256")
printf 'Binary archive: %s\n' "$archive"
