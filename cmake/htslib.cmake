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
set(HTSLIB_ZLIB_INC)
set(HTSLIB_ZLIB_LIB)

if ((NOT ZLIB_FOUND))
  if(NOT EXISTS "${zlib_INSTALL}")
    include(${CMAKE_SOURCE_DIR}/cmake/zlib.cmake)
  endif()
  message(STATUS "Using local build of zlib at ${zlib_INSTALL} with htslib")
  set(HTSLIB_ZLIB_INC "CPPFLAGS=-I${zlib_INSTALL}/include")
  set(HTSLIB_ZLIB_LIB "LDFLAGS=-L${zlib_INSTALL}/lib/")
endif()

ExternalProject_Add(htslib
    PREFIX ${htslib_PREFIX}
    GIT_REPOSITORY "https://github.com/samtools/htslib.git"
    GIT_TAG "1.20"
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    #CONFIGURE_COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/contrib/htslib-prefix/src/htslib/configure"
    #CONFIGURE_COMMAND "autoheader"
    #CONFIGURE_COMMAND "autoconf"
    CONFIGURE_COMMAND autoheader && autoconf && autoreconf --install && ./configure --disable-bz2 --disable-lzma --disable-libcurl --without-libdeflate ${HTSLIB_ZLIB_INC} ${HTSLIB_ZLIB_LIB}
    BUILD_COMMAND ${MAKE_COMMAND} lib-static
    INSTALL_COMMAND ${MAKE_COMMAND} install prefix=${htslib_INSTALL}
    LOG_DOWNLOAD 1
    )

if (NOT ZLIB_FOUND)
  message(STATUS "Adding zlib built via external project as dependency for htslib")
  add_dependencies(htslib zlib)
endif()

include_directories(${htslib_INSTALL}/include)
set(HTSlib_LIBRARIES ${htslib_INSTALL}/lib/libhts.a)
