set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib-prefix)
set(htslib_INSTALL ${CMAKE_BINARY_DIR}/contrib/htslib-install)

if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
    # invoke make explicitly
    # in this case, we assume the parent build system is running in parallel already so no -j flag is added
    find_program(MAKE_COMMAND NAMES make gmake)
endif()
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
set(HTSLIB_CPPFLAGS "-I${LIBDEFLATE_INCLUDE_DIR}")
get_filename_component(LIBDEFLATE_LIBRARY_DIR "${LIBDEFLATE_LIBRARY}" DIRECTORY)
set(HTSLIB_LDFLAGS "-L${LIBDEFLATE_LIBRARY_DIR}")

if ((NOT ZLIB_FOUND))
  if(NOT TARGET zlib)
    include(${CMAKE_SOURCE_DIR}/cmake/zlib.cmake)
  endif()
  message(STATUS "Using local build of zlib at ${zlib_INSTALL} with htslib")
  string(APPEND HTSLIB_CPPFLAGS " -I${zlib_INSTALL}/include")
  string(APPEND HTSLIB_LDFLAGS " -L${zlib_INSTALL}/lib")
endif()

ExternalProject_Add(htslib
    PREFIX ${htslib_PREFIX}
    GIT_REPOSITORY "https://github.com/samtools/htslib.git"
    # Resolved commit for HTSlib 1.20; pinning makes reviewer builds stable.
    GIT_TAG "0cadce238af0c6398751999bad703d4b19615860"
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    #CONFIGURE_COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/contrib/htslib-prefix/src/htslib/configure"
    #CONFIGURE_COMMAND "autoheader"
    #CONFIGURE_COMMAND "autoconf"
    # No --without-libdeflate: htslib's configure probes for libdeflate on its own
    # and silently falls back to zlib when it is absent. BGZF inflate is the single
    # hottest thing RabbitBin does -- the BAM->depth stage runs at the decompression
    # roofline -- so forcing the zlib path here roughly halves whole-pipeline
    # throughput on any machine that has to use this vendored build.
    CONFIGURE_COMMAND autoheader && autoconf && autoreconf --install &&
      ${CMAKE_COMMAND} -E env
        "CPPFLAGS=${HTSLIB_CPPFLAGS}"
        "LDFLAGS=${HTSLIB_LDFLAGS}"
        ./configure --disable-bz2 --disable-lzma --disable-libcurl
    BUILD_COMMAND ${MAKE_COMMAND} lib-static
    INSTALL_COMMAND ${MAKE_COMMAND} install prefix=${htslib_INSTALL}
    LOG_DOWNLOAD 1
    BUILD_BYPRODUCTS ${htslib_INSTALL}/lib/libhts.a
    )

if (NOT ZLIB_FOUND)
  message(STATUS "Adding zlib built via external project as dependency for htslib")
  add_dependencies(htslib zlib)
endif()

if(TARGET libdeflate_external)
  add_dependencies(htslib libdeflate_external)
endif()

include_directories(${htslib_INSTALL}/include)
set(HTSlib_INCLUDE_DIR ${htslib_INSTALL}/include)
set(HTSlib_INCLUDE_DIRS ${htslib_INSTALL}/include)
set(HTSlib_LIBRARIES ${htslib_INSTALL}/lib/libhts.a)
set(HTSlib_FOUND TRUE)
