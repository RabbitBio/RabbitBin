# Static version from VERSION file (no git dependency).
function(GET_VERSION PREFIX)
  file(READ "${CMAKE_SOURCE_DIR}/VERSION" _VER_RAW LIMIT_COUNT 1)
  string(STRIP "${_VER_RAW}" _VER_STRIPPED)
  set(RB_VERSION_STRING "${_VER_STRIPPED}")
  string(TIMESTAMP RB_BUILD_DATE "%Y%m%d_%H%M%S")

  set(VERSION_DIR "${PROJECT_BINARY_DIR}/makeVersionFile")
  file(MAKE_DIRECTORY "${VERSION_DIR}")

  message(STATUS "Building ${PREFIX} version ${RB_VERSION_STRING}")

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
