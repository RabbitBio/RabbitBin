if(NOT DEFINED RABBITBIN OR NOT DEFINED WORK_DIR)
  message(FATAL_ERROR "RABBITBIN and WORK_DIR are required")
endif()

file(MAKE_DIRECTORY "${WORK_DIR}")
set(FASTA "${WORK_DIR}/min_bin_size.fa")
set(DEPTH "${WORK_DIR}/min_bin_size.depth.tsv")
set(PREFIX "${WORK_DIR}/min_bin_size")

function(assert_strict_floor PREFIX_VALUE FLOOR_VALUE)
  set(BINS_FILE "${PREFIX_VALUE}.bins.tsv")
  set(MEMBERS_FILE "${PREFIX_VALUE}.members.tsv")
  foreach(OUTPUT_FILE IN ITEMS "${BINS_FILE}" "${MEMBERS_FILE}")
    if(NOT EXISTS "${OUTPUT_FILE}")
      message(FATAL_ERROR "Missing output: ${OUTPUT_FILE}")
    endif()
  endforeach()

  file(STRINGS "${BINS_FILE}" BIN_LINES)
  list(LENGTH BIN_LINES BIN_LINE_COUNT)
  if(BIN_LINE_COUNT LESS 1)
    message(FATAL_ERROR "Missing header in ${BINS_FILE}")
  endif()
  list(REMOVE_AT BIN_LINES 0)
  foreach(LINE IN LISTS BIN_LINES)
    string(REPLACE "\t" ";" FIELDS "${LINE}")
    list(GET FIELDS 2 BIN_LENGTH)
    if(BIN_LENGTH LESS FLOOR_VALUE)
      message(FATAL_ERROR
        "${BINS_FILE} emitted a ${BIN_LENGTH}-bp bin below --min-bin-size ${FLOOR_VALUE}")
    endif()
  endforeach()
endfunction()

# The three coverage patterns form internal communities below 20 kb. The output
# layer must suppress every bin below the requested floor.
set(SEQUENCE "")
foreach(REPEAT_INDEX RANGE 1 100)
  string(APPEND SEQUENCE "ACGTTGCATGTCAGTACCGATGCACTGATCGA")
endforeach()
set(FASTA_TEXT "")
set(DEPTH_TEXT
    "contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\ts2\ts2-var\ts3\ts3-var\ts4\ts4-var\n")
foreach(I RANGE 1 13)
  string(APPEND FASTA_TEXT ">contig_${I}\n${SEQUENCE}\n")
  if(I LESS_EQUAL 4)
    string(APPEND DEPTH_TEXT "contig_${I}\t3200\t2.5\t1\t0\t2\t0\t3\t0\t4\t0\n")
  elseif(I LESS_EQUAL 7)
    string(APPEND DEPTH_TEXT "contig_${I}\t3200\t2.5\t1\t0\t2\t0\t4\t0\t3\t0\n")
  else()
    string(APPEND DEPTH_TEXT "contig_${I}\t3200\t2.5\t4\t0\t3\t0\t2\t0\t1\t0\n")
  endif()
endforeach()
file(WRITE "${FASTA}" "${FASTA_TEXT}")
file(WRITE "${DEPTH}" "${DEPTH_TEXT}")

execute_process(
  COMMAND "${RABBITBIN}" bin
          --assembly "${FASTA}"
          --depth "${DEPTH}"
          --output "${PREFIX}"
          --min-bin-size 20000
          --min-contig 2500
          --threads 1
          --seed 42
          --no-bin-fasta
          --quiet
  RESULT_VARIABLE RC
  OUTPUT_VARIABLE STDOUT
  ERROR_VARIABLE STDERR
)
if(NOT RC EQUAL 0)
  message(FATAL_ERROR "RabbitBin failed (${RC})\nstdout:\n${STDOUT}\nstderr:\n${STDERR}")
endif()
assert_strict_floor("${PREFIX}" 20000)

# With one sample the graph is composition-only, so the identical contigs form
# one 41.6-kb parent. The marker multiplicity then forces a two-way abundance
# split into 19.2-kb and 22.4-kb children, both below the requested 30-kb floor.
set(MARKER_DEPTH "${WORK_DIR}/marker_min_bin_size.depth.tsv")
set(MARKER_SEED "${WORK_DIR}/marker_min_bin_size.seed.tsv")
set(MARKER_PREFIX "${WORK_DIR}/marker_min_bin_size")
set(MARKER_DEPTH_TEXT
    "contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\n")
foreach(I RANGE 1 13)
  if(I LESS_EQUAL 6)
    string(APPEND MARKER_DEPTH_TEXT "contig_${I}\t3200\t1\t1\t0\n")
  else()
    string(APPEND MARKER_DEPTH_TEXT "contig_${I}\t3200\t10\t10\t0\n")
  endif()
endforeach()
file(WRITE "${MARKER_DEPTH}" "${MARKER_DEPTH_TEXT}")
file(WRITE "${MARKER_SEED}" "marker_A\tcontig_1\tcontig_7\n")

execute_process(
  COMMAND "${RABBITBIN}" bin
          --assembly "${FASTA}"
          --depth "${MARKER_DEPTH}"
          --output "${MARKER_PREFIX}"
          --min-bin-size 30000
          --min-contig 2500
          --threads 1
          --seed 42
          --marker-seed "${MARKER_SEED}"
          --no-bin-fasta
          --quiet
  RESULT_VARIABLE MARKER_RC
  OUTPUT_VARIABLE MARKER_STDOUT
  ERROR_VARIABLE MARKER_STDERR
)
if(NOT MARKER_RC EQUAL 0)
  message(FATAL_ERROR
    "Marker-guided RabbitBin failed (${MARKER_RC})\nstdout:\n${MARKER_STDOUT}\nstderr:\n${MARKER_STDERR}")
endif()
if(NOT "${MARKER_STDOUT}${MARKER_STDERR}" MATCHES "Marker-guided split:[^\n]*1 split")
  message(FATAL_ERROR
    "Fixture did not exercise marker-guided splitting\nstdout:\n${MARKER_STDOUT}\nstderr:\n${MARKER_STDERR}")
endif()
assert_strict_floor("${MARKER_PREFIX}" 30000)
