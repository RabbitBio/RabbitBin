# * Try to find htslib Once done, this will define
#
# htslib_FOUND - system has htslib htslib_INCLUDE_DIRS - the htslib include directories
# htslib_LIBRARIES - link these to use htslib
#
# This code was modified from
# https://github.com/genome/build-common/blob/master/cmake/FindHTSlib.cmake Source:
# https://github.com/luntergroup/octopus/blob/develop/build/cmake/modules/FindHTSlib.cmake

# A simple wrapper to make pkg-config searches a bit easier. Works the same as CMake's internal
# pkg_check_modules but is always quiet.
macro(libfind_pkg_check_modules)
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(${ARGN} QUIET)
  endif()
endmacro()

macro(libfind_package PREFIX)
  set(LIBFIND_PACKAGE_ARGS ${ARGN})
  if(${PREFIX}_FIND_QUIETLY)
    set(LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} QUIET)
  endif(${PREFIX}_FIND_QUIETLY)
  if(${PREFIX}_FIND_REQUIRED)
    set(LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} REQUIRED)
  endif(${PREFIX}_FIND_REQUIRED)
  find_package(${LIBFIND_PACKAGE_ARGS})
endmacro(libfind_package)

macro(libfind_process PREFIX)
  # Skip processing if already processed during this run
  if(NOT ${PREFIX}_FOUND)
    # Start with the assumption that the library was found
    set(${PREFIX}_FOUND TRUE)

    # Process all includes and set _FOUND to false if any are missing
    foreach(i ${${PREFIX}_PROCESS_INCLUDES})
      if(${i})
        set(${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIRS} ${${i}})
        mark_as_advanced(${i})
      else(${i})
        set(${PREFIX}_FOUND FALSE)
      endif(${i})
    endforeach(i)

    # Process all libraries and set _FOUND to false if any are missing
    foreach(i ${${PREFIX}_PROCESS_LIBS})
      if(${i})
        set(${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARIES} ${${i}})
        mark_as_advanced(${i})
      else(${i})
        set(${PREFIX}_FOUND FALSE)
      endif(${i})
    endforeach(i)

    # Print message and/or exit on fatal error
    if(${PREFIX}_FOUND)
      if(NOT ${PREFIX}_FIND_QUIETLY)
        message(STATUS "Found ${PREFIX} ${${PREFIX}_VERSION}")
      endif(NOT ${PREFIX}_FIND_QUIETLY)
    else(${PREFIX}_FOUND)
      if(${PREFIX}_FIND_REQUIRED)
        foreach(i ${${PREFIX}_PROCESS_INCLUDES} ${${PREFIX}_PROCESS_LIBS})
          message("${i}=${${i}}")
        endforeach(i)
        message(
          FATAL_ERROR
            "Required library ${PREFIX} NOT FOUND.\nInstall the library (dev version) and try again. If the library is already installed, use ccmake to set the missing variables manually."
        )
      endif(${PREFIX}_FIND_REQUIRED)
    endif(${PREFIX}_FOUND)
  endif(NOT ${PREFIX}_FOUND)
endmacro(libfind_process)

set(_htslib_roots ${HTSLIB_ROOT} $ENV{HTSLIB_ROOT} ${HTSLIB_SEARCH_DIRS})
set(_htslib_find_options "")
if(HTSlib_NO_SYSTEM_PATHS)
  # Keep explicitly supplied prefixes usable even when system lookup is off.
  file(TO_CMAKE_PATH "$ENV{CMAKE_PREFIX_PATH}" _htslib_env_prefixes)
  list(APPEND _htslib_roots ${CMAKE_PREFIX_PATH} ${_htslib_env_prefixes})
  set(_htslib_find_options NO_DEFAULT_PATH)
endif()

set(_htslib_ver_path "htslib-${HTSlib_FIND_VERSION}")

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(HTSLIB_PKGCONF htslib)

# Include dir
find_path(
  HTSlib_INCLUDE_DIR
  NAMES ${HTSLIB_ADDITIONAL_HEADERS} htslib/sam.h
  HINTS ${_htslib_roots} ${HTSLIB_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES include htslib/${_htslib_ver_path}
  ${_htslib_find_options}
)

if(HTSlib_USE_STATIC_LIBS)
  # Dependencies
  set(ZLIB_ROOT ${HTSLIB_ROOT})
  libfind_package(HTSlib ZLIB)
  set(BZip2_ROOT ${HTSLIB_ROOT})
  libfind_package(HTSlib BZip2)
  set(LibLZMA_ROOT ${HTSLIB_ROOT})
  libfind_package(HTSlib LibLZMA)
  set(CURL_ROOT ${HTSLIB_ROOT})
  libfind_package(HTSlib CURL)
  if(NOT APPLE)
    set(OpenSSL_ROOT ${HTSLIB_ROOT})
    libfind_package(HTSlib OpenSSL)
  endif()
  set(HTSlib_LIBRARY_names libhts.a)
else()
  set(HTSlib_LIBRARY_names libhts.so libhts.so.2 libhts.dylib libhts.2.dylib)
endif()

# Finally the library itself
find_library(
  HTSlib_LIBRARY
  NAMES ${HTSlib_LIBRARY_names}
  HINTS ${_htslib_roots} ${HTSLIB_PKGCONF_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu ${_htslib_ver_path}
  ${_htslib_find_options}
)

# Set the include dir variables and the libraries and let libfind_process do the rest. NOTE:
# Singular variables for this library, plural for libraries this lib depends on.
set(HTSlib_PROCESS_INCLUDES HTSlib_INCLUDE_DIR)
set(HTSlib_PROCESS_LIBS HTSlib_LIBRARY)

if(HTSlib_USE_STATIC_LIBS)
  set(HTSlib_PROCESS_INCLUDES ${HTSlib_PROCESS_INCLUDES} ZLIB_INCLUDE_DIR BZIP2_INCLUDE_DIR
                              LIBLZMA_INCLUDE_DIRS CURL_INCLUDE_DIRS
  )
  set(HTSlib_PROCESS_LIBS ${HTSlib_PROCESS_LIBS} ZLIB_LIBRARIES BZIP2_LIBRARIES LIBLZMA_LIBRARIES
                          CURL_LIBRARIES
  )
  if(NOT APPLE)
    set(HTSlib_PROCESS_INCLUDES ${HTSlib_PROCESS_INCLUDES} OPENSSL_INCLUDE_DIR)
    set(HTSlib_PROCESS_LIBS ${HTSlib_PROCESS_LIBS} OPENSSL_LIBRARIES)
  endif()
endif()

# Read the selected headers, not a possibly unrelated system pkg-config file.
# HTS_VERSION encodes X.Y.Z as XYYYZZ (including development snapshots).
set(HTSlib_VERSION "")
if(EXISTS "${HTSlib_INCLUDE_DIR}/htslib/hts.h")
  file(STRINGS "${HTSlib_INCLUDE_DIR}/htslib/hts.h" _htslib_version_line
       REGEX "^[ \t]*#[ \t]*define[ \t]+HTS_VERSION[ \t]+[0-9]+")
  if(_htslib_version_line MATCHES "HTS_VERSION[ \t]+([0-9]+)")
    set(_htslib_version_number "${CMAKE_MATCH_1}")
    math(EXPR _htslib_major "${_htslib_version_number} / 100000")
    math(EXPR _htslib_minor "(${_htslib_version_number} / 100) % 1000")
    math(EXPR _htslib_patch "${_htslib_version_number} % 100")
    set(HTSlib_VERSION "${_htslib_major}.${_htslib_minor}.${_htslib_patch}")
  endif()
endif()
libfind_process(HTSlib)

if(HTSlib_FOUND AND HTSlib_FIND_VERSION AND HTSlib_VERSION AND
   HTSlib_VERSION VERSION_LESS HTSlib_FIND_VERSION)
  message(STATUS "Found HTSlib ${HTSlib_VERSION}, but ${HTSlib_FIND_VERSION} or newer is required")
  set(HTSlib_FOUND FALSE)
endif()

message(STATUS "   HTSlib include dirs: ${HTSlib_INCLUDE_DIRS}")
message(STATUS "   HTSlib libraries: ${HTSlib_LIBRARIES}")
message(STATUS " HTSLIB FOUND: ${HTSlib_FOUND} ")
