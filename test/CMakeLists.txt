
set(CTEST_OUTPUT_ON_FAILURE 1)

add_executable(test_silhouette test_silhouette.cpp)
add_test(NAME test_silhouette COMMAND test_silhouette)

add_executable(test_recruitment test_recruitment.cpp)
add_test(NAME test_recruitment COMMAND test_recruitment)

find_package(PythonInterp 3 QUIET)
if(PYTHONINTERP_FOUND)
  add_test(NAME test_coverage_edges
           COMMAND ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_coverage_edges.py
                   ${CMAKE_BINARY_DIR}/src/rabbitbin
                   ${CMAKE_BINARY_DIR}/test-coverage-edges)
else()
  message(STATUS "Python 3 unavailable; coverage edge integration test not registered")
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
