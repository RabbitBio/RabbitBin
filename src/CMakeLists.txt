# find and use Boost
set(Boost_NO_BOOST_CMAKE ON) # for backwards comptibility of cmake and boost >= 1.70
set(Boost_USE_STATIC_LIBS   OFF)
find_package(Boost 1.66.0 COMPONENTS program_options filesystem system graph serialization iostreams REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )

# find and use OpenMP
find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    message("Using OpenMP")
else()
    message(WARNING "OpenMP must be available to your compiler or it will be REALLLLLLLLLLY slow: ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER}")
endif()

if(APPLE)
	#full Single Unix Standard v3 (SUSv3) conformance (the Unix API)
	add_definitions(-D_DARWIN_C_SOURCE)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -qreport -qhot -qmaxmem=819200 -qoffload -qtgtarch=sm_70 -g")
set(EXTRALIBS ${EXTRALIBS} )
endif()

message("Boost cmake dir: ${Boost_DIR}")
message("Boost include: ${Boost_INCLUDE_DIRS}")
message("Boost libraries:" ${Boost_LIBRARIES})
message("Boost lib ver:" ${Boost_LIB_VERSION})
message("Boost lib dir:" ${Boost_LIBRARY_DIRS})

# RabbitBin: sketch-based metagenome binning
add_executable(rabbitbin rabbitbin.cpp probmh.cpp)
if(NOT "${EXTRA_DEPENDENCIES}" STREQUAL "")
  add_dependencies(rabbitbin ${EXTRA_DEPENDENCIES})
endif()
target_link_libraries(rabbitbin ${ZLIB_LIBRARIES} ${Boost_LIBRARIES} ${EXTRALIBS} RabbitBin_VERSION_LIB deflate)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  target_compile_options(rabbitbin PRIVATE -march=native)
  message("rabbitbin: -march=native enabled (AVX-512 VPOPCNTDQ if available)")
else()
  check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
  if(COMPILER_SUPPORTS_AVX2)
    target_compile_options(rabbitbin PRIVATE -mavx2)
    message("rabbitbin: AVX2 SIMD enabled")
  else()
    message("rabbitbin: AVX2 not supported, using scalar fallback")
  endif()
endif()
# Link-time optimization: rabbitbin is split across rabbitbin.cpp + probmh.cpp,
# so the hot ProbMinHash methods (addHashFromRng, tracker update, jaccard) live
# in a different translation unit from their per-contig sketch caller and cannot
# be inlined without LTO. Enabling IPO lets the compiler inline + specialise the
# sketch hot path across the TU boundary. No -ffast-math is implied, so floating
# point results are unchanged. Only applied if the toolchain supports it.
check_cxx_compiler_flag("-flto" COMPILER_SUPPORTS_LTO)
if(COMPILER_SUPPORTS_LTO)
  # Apply LTO explicitly to both compile and link so the objects are emitted as
  # LTO bytecode and the final codegen sees both TUs. -march=native is re-asserted
  # on the link line because LTO performs code generation at link time.
  target_compile_options(rabbitbin PRIVATE -flto)
  target_link_options(rabbitbin PRIVATE -flto)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    target_link_options(rabbitbin PRIVATE -march=native)
  endif()
  message("rabbitbin: LTO enabled")
else()
  message("rabbitbin: LTO not supported by compiler")
endif()

install(TARGETS rabbitbin DESTINATION bin/)

set(targets rabbit_depth rabbit_overlap)
foreach(target ${targets})
   add_executable(${target} ${target}.cpp)
   if (NOT "${EXTRA_DEPENDENCIES}" STREQUAL "")
     add_dependencies(${target} ${EXTRA_DEPENDENCIES})
   endif()
   target_link_libraries(${target} ${HTSlib_LIBRARIES} ${ZLIB_LIBRARIES} ${Boost_LIBRARIES} RabbitBin_VERSION_LIB) 
endforeach()

install(TARGETS ${targets}
        DESTINATION bin/
       )
