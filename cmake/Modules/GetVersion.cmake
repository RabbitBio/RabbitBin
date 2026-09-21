# Release version comes from VERSION.  When configuring from a Git checkout,
# also record the exact source commit so benchmark and reviewer builds can be
# identified unambiguously.  Release tarballs remain buildable without Git.
function(GET_VERSION PREFIX)
  file(READ "${CMAKE_SOURCE_DIR}/VERSION" _VER_RAW LIMIT_COUNT 1)
  string(STRIP "${_VER_RAW}" _VER_STRIPPED)
  set(RB_VERSION_STRING "${_VER_STRIPPED}")
  string(TIMESTAMP RB_BUILD_DATE "%Y%m%d_%H%M%S")
  set(RB_GIT_COMMIT "unknown")
  if(DEFINED RABBITBIN_SOURCE_REVISION AND
     NOT "${RABBITBIN_SOURCE_REVISION}" STREQUAL "")
    set(RB_GIT_COMMIT "${RABBITBIN_SOURCE_REVISION}")
  else()
    find_package(Git QUIET)
  endif()
  if(RB_GIT_COMMIT STREQUAL "unknown" AND GIT_FOUND AND
     EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE _GIT_RESULT
      OUTPUT_VARIABLE _GIT_COMMIT
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(_GIT_RESULT EQUAL 0 AND NOT "${_GIT_COMMIT}" STREQUAL "")
      set(RB_GIT_COMMIT "${_GIT_COMMIT}")
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" status --porcelain
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE _GIT_DIRTY
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      if(NOT "${_GIT_DIRTY}" STREQUAL "")
        string(APPEND RB_GIT_COMMIT "-dirty")
      endif()
    endif()
  endif()

  set(VERSION_DIR "${PROJECT_BINARY_DIR}/makeVersionFile")
  file(MAKE_DIRECTORY "${VERSION_DIR}")

  message(STATUS "Building ${PREFIX} version ${RB_VERSION_STRING} (${RB_GIT_COMMIT})")

  configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/version.cpp.in"
    "${VERSION_DIR}/version.cpp"
    @ONLY
  )

  include_directories(${VERSION_DIR})
  add_library(${PREFIX}_VERSION OBJECT "${VERSION_DIR}/version.cpp")
  add_library(${PREFIX}_VERSION_LIB $<TARGET_OBJECTS:${PREFIX}_VERSION>)
  set_target_properties(${PREFIX}_VERSION_LIB PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
endfunction()
