# No real HTSlib, pkg-config, compiler or network is needed for these lookups.
file(MAKE_DIRECTORY "${WORK_DIR}/prefix with spaces/include/htslib"
                    "${WORK_DIR}/prefix with spaces/lib" "${WORK_DIR}/project")
set(PREFIX "${WORK_DIR}/prefix with spaces")
file(WRITE "${PREFIX}/include/htslib/sam.h" "/* finder fixture */\n")
file(WRITE "${PREFIX}/include/htslib/hts.h" "#define HTS_VERSION 102400\n")
file(WRITE "${PREFIX}/lib/libhts.so" "")
file(WRITE "${WORK_DIR}/project/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.16)
project(HTSlibPrefixTest NONE)
list(PREPEND CMAKE_MODULE_PATH "${MODULE_DIR}")
set(CMAKE_DISABLE_FIND_PACKAGE_PkgConfig TRUE)
find_package(HTSlib 1.13 QUIET)
if(EXPECT_OLD)
  if(HTSlib_FOUND)
    message(FATAL_ERROR "Accepted HTSlib below the minimum version")
  endif()
else()
  if(NOT HTSlib_FOUND OR NOT HTSlib_INCLUDE_DIR STREQUAL "${EXPECTED_PREFIX}/include"
     OR NOT HTSlib_LIBRARY STREQUAL "${EXPECTED_PREFIX}/lib/libhts.so"
     OR NOT HTSlib_VERSION STREQUAL "1.24.0")
    message(FATAL_ERROR "Did not select the requested HTSlib prefix/version")
  endif()
endif()
]=])

foreach(CASE IN ITEMS prefix root restricted)
  set(OPTIONS "-DCMAKE_PREFIX_PATH=${PREFIX}")
  if(CASE STREQUAL "root")
    set(OPTIONS "-DHTSLIB_ROOT=${PREFIX}")
  elseif(CASE STREQUAL "restricted")
    list(APPEND OPTIONS -DHTSlib_NO_SYSTEM_PATHS=ON)
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${WORK_DIR}/project" -B "${WORK_DIR}/${CASE}"
            "-DMODULE_DIR=${MODULE_DIR}" "-DEXPECTED_PREFIX=${PREFIX}" ${OPTIONS}
    RESULT_VARIABLE RC OUTPUT_VARIABLE OUT ERROR_VARIABLE ERR)
  if(NOT RC EQUAL 0)
    message(FATAL_ERROR "HTSlib ${CASE} lookup failed\n${OUT}\n${ERR}")
  endif()
endforeach()

file(WRITE "${PREFIX}/include/htslib/hts.h" "#define HTS_VERSION 101200\n")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -S "${WORK_DIR}/project" -B "${WORK_DIR}/old"
          "-DMODULE_DIR=${MODULE_DIR}" "-DEXPECTED_PREFIX=${PREFIX}"
          "-DCMAKE_PREFIX_PATH=${PREFIX}" -DEXPECT_OLD=ON
  RESULT_VARIABLE RC OUTPUT_VARIABLE OUT ERROR_VARIABLE ERR)
if(NOT RC EQUAL 0)
  message(FATAL_ERROR "HTSlib minimum-version check failed\n${OUT}\n${ERR}")
endif()
