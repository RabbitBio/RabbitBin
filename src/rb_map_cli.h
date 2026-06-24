#ifndef RB_MAP_CLI_H_
#define RB_MAP_CLI_H_

// Subcommand entry points for the align -> sort -> index pipeline.
// Each is defined in its own translation unit (src/align, src/bamsort,
// src/bamindex, src/rb_map.cpp) and registered in rabbitbin.cpp main().
// Only compiled/linked when RABBITBIN_ENABLE_MAP is defined (needs htslib).

int rb_cmd_bwa(int ac, char *av[]);      // src/align/rb_align.cpp
int rb_cmd_sortbam(int ac, char *av[]);  // src/bamsort/rb_sort.cpp
int rb_cmd_bai(int ac, char *av[]);      // src/bamindex/rb_index.cpp
int rb_cmd_map(int ac, char *av[]);      // src/rb_map.cpp

#endif  // RB_MAP_CLI_H_
