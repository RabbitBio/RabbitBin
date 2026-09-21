cmake_minimum_required(VERSION 3.16)
project (RabbitBin)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Installing ${CMAKE_BUILD_TYPE} RabbitBin into ${CMAKE_INSTALL_PREFIX}")

# Preserve non-system dependency locations (for example a user-managed
# HTSlib/libdeflate prefix) in installed executables.  System library
# directories are omitted automatically by CMake.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules )
include(GetVersion)

set(EXTRA_DEPENDENCIES "")

option(RABBITBIN_USE_SYSTEM_ZLIB "Try to find and use system ZLIB" ON)
if(RABBITBIN_USE_SYSTEM_ZLIB)
  find_package(ZLIB 1.2.11)
  if(NOT ZLIB_FOUND)
    message(STATUS "Could not find a system version of zlib")
  endif()
else()
  message(STATUS "Skipping system detection of system install of zlib")
endif()
if (NOT ZLIB_FOUND)
  message(STATUS "Including zlib via add_external_project")
  include(${PROJECT_SOURCE_DIR}/cmake/zlib.cmake)
  set(EXTRA_DEPENDENCIES ${EXTRA_DEPENDENCIES} zlib)
else()
  message(STATUS "Using system zlib")
endif()

# libdeflate is used by the main executable for fast .gz FASTA input and by the
# vendored RabbitBAM reader.  It is therefore a real build dependency, not only
# an optional BAM acceleration.  Download a pinned copy when no development
# installation is available.
option(RABBITBIN_USE_SYSTEM_LIBDEFLATE "Try to find and use system libdeflate" ON)
if(RABBITBIN_USE_SYSTEM_LIBDEFLATE)
  find_path(LIBDEFLATE_INCLUDE_DIR
    NAMES libdeflate.h
    HINTS ${LIBDEFLATE_ROOT} $ENV{LIBDEFLATE_ROOT}
    PATH_SUFFIXES include)
  find_library(LIBDEFLATE_LIBRARY
    NAMES deflate
    HINTS ${LIBDEFLATE_ROOT} $ENV{LIBDEFLATE_ROOT}
    PATH_SUFFIXES lib lib64)
endif()
if(LIBDEFLATE_INCLUDE_DIR AND LIBDEFLATE_LIBRARY)
  set(LIBDEFLATE_FOUND TRUE)
  message(STATUS "Using system libdeflate: ${LIBDEFLATE_LIBRARY}")
else()
  message(STATUS "System libdeflate development files not found; downloading a pinned copy")
  include(${PROJECT_SOURCE_DIR}/cmake/libdeflate.cmake)
  list(APPEND EXTRA_DEPENDENCIES libdeflate_external)
endif()

option(RABBITBIN_USE_SYSTEM_HTSLIB "Try to find and use system HTSLIB" ON)
if(RABBITBIN_USE_SYSTEM_HTSLIB)
  find_package(HTSlib 1.13)
  if (NOT HTSlib_FOUND)
    message(STATUS "Could not find a system installation of htslib")
  endif()
else()
  message(STATUS "Skipping detection of system install of htslib")
endif()
if (NOT HTSlib_FOUND)
  message(STATUS "Including htslib via add_external_project")
  include(${PROJECT_SOURCE_DIR}/cmake/htslib.cmake)
  set(EXTRA_DEPENDENCIES ${EXTRA_DEPENDENCIES} htslib)
else()
  message(STATUS "using system htslib")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# C11 (not C17): the only C source here is the vendored ssw.c, which needs no
# C17 feature, and CMake 3.20 lacks the GNU C17 dialect-flag mapping (added in
# 3.21) — forcing C17 breaks configure on the system toolchain.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_definitions(-D_XOPEN_SOURCE=700)

option(RABBITBIN_NATIVE_ARCH
  "Tune for the build host CPU (faster locally, non-portable binaries)" OFF)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using AppleClang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (RABBITBIN_NATIVE_ARCH AND "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mtune=native")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=native")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
endif()

include_directories(${CMAKE_SOURCE_DIR})

get_version(RabbitBin)

add_subdirectory(src)

if (NOT NO_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

install(PROGRAMS run_rabbitbin.sh filter_sam_by_min_len.py aggregateBinDepths.pl aggregateContigOverlapsByBin.pl
                 merge_abundances.py scripts/rabbitbin_markers.sh
        DESTINATION bin/
       )
