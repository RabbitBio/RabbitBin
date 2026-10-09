
set(CTEST_OUTPUT_ON_FAILURE 1)

add_executable(test_silhouette test_silhouette.cpp)
add_test(NAME test_silhouette COMMAND test_silhouette)

add_executable(test_recruitment test_recruitment.cpp)
add_test(NAME test_recruitment COMMAND test_recruitment)

add_executable(test_abundance_filter test_abundance_filter.cpp)
add_test(NAME test_abundance_filter COMMAND test_abundance_filter)

add_executable(test_bam_aux test_bam_aux.cpp)
find_package(Threads REQUIRED)
target_include_directories(test_bam_aux PRIVATE ${HTSlib_INCLUDE_DIRS})
target_link_libraries(test_bam_aux ${HTSlib_LIBRARIES} ${ZLIB_LIBRARIES}
                     Threads::Threads ${CMAKE_DL_LIBS})
if(DEFLATE_LIB)
  target_link_libraries(test_bam_aux ${DEFLATE_LIB})
endif()
if(NOT "${EXTRA_DEPENDENCIES}" STREQUAL "")
  add_dependencies(test_bam_aux ${EXTRA_DEPENDENCIES})
endif()
add_test(NAME test_bam_aux COMMAND test_bam_aux)

find_package(Python3 3.6 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME test_coverage_edges
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_coverage_edges.py
                   ${CMAKE_BINARY_DIR}/src/rabbitbin
                   ${CMAKE_BINARY_DIR}/test-coverage-edges)
  add_test(NAME test_abundance_filter_graph
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_abundance_filter_graph.py
                   ${CMAKE_BINARY_DIR}/src/rabbitbin
                   ${CMAKE_BINARY_DIR}/test-abundance-filter-graph)
  find_program(SAMTOOLS_EXECUTABLE samtools)
  if(SAMTOOLS_EXECUTABLE AND DEFLATE_LIB)
    add_test(NAME test_bam_depth_equivalence
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/test_bam_depth_equivalence.py
                     ${CMAKE_BINARY_DIR}/src/rabbitbin
                     ${CMAKE_BINARY_DIR}/test-bam-depth-equivalence
                     ${SAMTOOLS_EXECUTABLE})
  else()
    message(STATUS "samtools/libdeflate unavailable; BAM depth equivalence test not registered")
  endif()
  add_test(NAME test_wrapper_paths
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_wrapper_paths.py
                   ${CMAKE_BINARY_DIR}/src/run_rabbitbin.sh
                   ${CMAKE_CURRENT_SOURCE_DIR}
                   ${CMAKE_BINARY_DIR}/test-wrapper-paths)
else()
  message(STATUS "Python 3.6+ unavailable; Python integration tests not registered")
endif()

file(MAKE_DIRECTORY  ${CMAKE_BINARY_DIR}/test-each)
add_test(NAME test_depths
         COMMAND ${CMAKE_BINARY_DIR}/src/rabbit_depth --referenceFasta ${CMAKE_CURRENT_SOURCE_DIR}/contigs.fa --outputDepth contigs_depth.txt ${CMAKE_CURRENT_SOURCE_DIR}/contigs-1000.fastq.bam
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-each
        )

file(MAKE_DIRECTORY  ${CMAKE_BINARY_DIR}/test-each)
add_test(NAME test_depths_gc_window_longer_than_contig
         COMMAND ${CMAKE_BINARY_DIR}/src/rabbit_depth --referenceFasta ${CMAKE_CURRENT_SOURCE_DIR}/contigs.fa --outputDepth contigs_depth_gc_window.txt --outputGC contigs_gc_window.txt --outputReadStats contigs_read_stats.txt --outputKmers contigs_kmers.txt --gcWindow 6000 ${CMAKE_CURRENT_SOURCE_DIR}/contigs-1000.fastq.bam
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-each
        )

file(MAKE_DIRECTORY  ${CMAKE_BINARY_DIR}/test-each)
add_test(NAME test_depths_malformed_softclip_only
         COMMAND ${CMAKE_BINARY_DIR}/src/rabbit_depth --referenceFasta ${CMAKE_CURRENT_SOURCE_DIR}/contigs.fa --outputDepth malformed_depth.txt --outputGC malformed_gc.txt --outputReadStats malformed_readstats.txt --outputKmers malformed_kmers.txt ${CMAKE_CURRENT_SOURCE_DIR}/malformed_softclip_only.bam
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-each
        )

file(MAKE_DIRECTORY  ${CMAKE_BINARY_DIR}/test-each)
add_test(NAME test_rabbitbin
         COMMAND ${CMAKE_BINARY_DIR}/src/rabbitbin --min-bin-size 0 --assembly ${CMAKE_CURRENT_SOURCE_DIR}/contigs.fa --depth contigs_depth.txt --output rabbitbin/contigs_bins -v
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-each
        )
set_tests_properties(test_depths PROPERTIES FIXTURES_SETUP contig_depth)
set_tests_properties(test_rabbitbin PROPERTIES FIXTURES_REQUIRED contig_depth)

add_test(NAME test_min_bin_size_strict
         COMMAND ${CMAKE_COMMAND}
                 -DRABBITBIN=${CMAKE_BINARY_DIR}/src/rabbitbin
                 -DWORK_DIR=${CMAKE_BINARY_DIR}/test-min-bin-size
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/test_min_bin_size.cmake
        )

file(COPY ${CMAKE_SOURCE_DIR}/run_rabbitbin.sh DESTINATION ${CMAKE_BINARY_DIR}/src)
file(MAKE_DIRECTORY  ${CMAKE_BINARY_DIR}/test-runRabbitBin)
add_test(NAME run_rabbitbin.sh
         COMMAND  ${CMAKE_BINARY_DIR}/src/run_rabbitbin.sh --min-bin-size 0 ${CMAKE_CURRENT_SOURCE_DIR}/contigs.fa ${CMAKE_CURRENT_SOURCE_DIR}/contigs-1000.fastq.bam
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-runRabbitBin
        )

add_test(NAME rabbitbin_version
         COMMAND ${CMAKE_COMMAND}
                 "-DRABBITBIN=${CMAKE_BINARY_DIR}/src/rabbitbin"
                 "-DEXPECTED_VERSION=${RabbitBin_VERSION_STRING}"
                 "-DEXPECTED_REVISION=${RabbitBin_GIT_COMMIT}"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/test_version.cmake)

add_test(NAME test_htslib_prefix
         COMMAND ${CMAKE_COMMAND}
                 "-DMODULE_DIR=${CMAKE_SOURCE_DIR}/cmake/Modules"
                 "-DWORK_DIR=${CMAKE_BINARY_DIR}/test-htslib-prefix"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/test_htslib_prefix.cmake)

add_test(NAME test_version_metadata
         COMMAND ${CMAKE_COMMAND}
                 "-DMODULE_DIR=${CMAKE_SOURCE_DIR}/cmake/Modules"
                 "-DWORK_DIR=${CMAKE_BINARY_DIR}/test-version-metadata"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/test_version_metadata.cmake)
