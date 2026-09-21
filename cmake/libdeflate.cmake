include(ExternalProject)

set(libdeflate_PREFIX "${CMAKE_BINARY_DIR}/contrib/libdeflate-prefix")
set(libdeflate_INSTALL "${CMAKE_BINARY_DIR}/contrib/libdeflate-install")

# Pin the resolved v1.26 commit rather than a movable branch or tag so a fresh
# reviewer build uses exactly the dependency revision tested by CI.
ExternalProject_Add(libdeflate_external
  PREFIX "${libdeflate_PREFIX}"
  GIT_REPOSITORY "https://github.com/ebiggers/libdeflate.git"
  GIT_TAG "92e6a0db9fa848d742f9eb286c92afc60f2c3dda"
  UPDATE_COMMAND ""
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=${libdeflate_INSTALL}
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DLIBDEFLATE_BUILD_SHARED_LIB=OFF
    -DLIBDEFLATE_BUILD_GZIP=OFF
    -DLIBDEFLATE_BUILD_TESTS=OFF
  BUILD_BYPRODUCTS "${libdeflate_INSTALL}/lib/libdeflate.a"
)

set(LIBDEFLATE_INCLUDE_DIR "${libdeflate_INSTALL}/include")
set(LIBDEFLATE_LIBRARY "${libdeflate_INSTALL}/lib/libdeflate.a")
set(LIBDEFLATE_FOUND TRUE)
