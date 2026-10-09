file(MAKE_DIRECTORY "${WORK_DIR}/source/cmake")
file(WRITE "${WORK_DIR}/source/VERSION" "1.0.0\n")
configure_file("${MODULE_DIR}/../version.cpp.in"
               "${WORK_DIR}/source/cmake/version.cpp.in" COPYONLY)
file(WRITE "${WORK_DIR}/source/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.16)
project(VersionMetadataTest CXX)
include("${MODULE_DIR}/GetVersion.cmake")
get_version(RabbitBin)
if(NOT RabbitBin_GIT_COMMIT STREQUAL EXPECTED)
  message(FATAL_ERROR "Expected '${EXPECTED}', got '${RabbitBin_GIT_COMMIT}'")
endif()
]=])

foreach(CASE IN ITEMS placeholder invalid archive override)
  set(EXPECTED unknown)
  set(OPTIONS "")
  if(CASE STREQUAL "placeholder")
    file(WRITE "${WORK_DIR}/source/SOURCE_REVISION" "$Format:%H$\n")
  elseif(CASE STREQUAL "invalid")
    file(WRITE "${WORK_DIR}/source/SOURCE_REVISION" "bad\n")
  else()
    file(WRITE "${WORK_DIR}/source/SOURCE_REVISION"
         "0123456789abcdef0123456789abcdef01234567\n")
    set(EXPECTED 0123456789ab)
    if(CASE STREQUAL "override")
      set(EXPECTED fedcba987654-dirty)
      set(OPTIONS "-DRABBITBIN_SOURCE_REVISION=${EXPECTED}")
    endif()
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${WORK_DIR}/source" -B "${WORK_DIR}/${CASE}"
            "-DMODULE_DIR=${MODULE_DIR}" "-DEXPECTED=${EXPECTED}" ${OPTIONS}
    RESULT_VARIABLE RC OUTPUT_VARIABLE OUT ERROR_VARIABLE ERR)
  if(NOT RC EQUAL 0)
    message(FATAL_ERROR "Version metadata ${CASE} failed\n${OUT}\n${ERR}")
  endif()
endforeach()
