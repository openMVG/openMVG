##
# SCOTCH library for hypergraph partitioning
##
#URL http://gforge.inria.fr/frs/download.php/zip/9109/scotch-6.0.4.zip

project(scotch CXX)

set(scotch_out_dir ${PROJECT_BINARY_DIR})

file(
  GLOB
  myScotch_headers
  src/libscotch/*.h
  )

file(
  GLOB
  myScotch_sources
  src/libscotch/arch*.c
  src/libscotch/bgraph*.c
  src/libscotch/common*.c
  src/libscotch/fibo.c
  src/libscotch/gain.c
  src/libscotch/geom.c
  src/libscotch/graph*.c
  src/libscotch/hall*.c
  src/libscotch/hgraph*.c
  src/libscotch/hmesh*.c
  src/libscotch/kgraph*.c
  src/libscotch/library_arch*.c
  src/libscotch/library_common_f.c
  src/libscotch/library_geom*.c
  src/libscotch/library_graph*.c
  src/libscotch/library_mapping*.c
  src/libscotch/library_memory*.c
  src/libscotch/library_mesh*.c
  src/libscotch/library_order.c
  src/libscotch/library_parser*.c
  src/libscotch/library_random*.c
  src/libscotch/library_strat*.c
  src/libscotch/library_version*.c
  src/libscotch/mapping*.c
  src/libscotch/mesh*.c
  src/libscotch/order*.c
  src/libscotch/parser*.c
  src/libscotch/vgraph*.c
  src/libscotch/vmesh*.c
  src/libscotch/wgraph*.c
  src/libscotch/library_error.c
  )

file(GLOB_RECURSE REMOVENONSOURCEDOTCFILES
  bgraph_bipart_df_loop.c
  common_sort.c
  graph_coarsen_edge.c
  graph_match_scan.c
  hgraph_induce_edge.c
  kgraph_map_df_loop.c
  )
list(REMOVE_ITEM myScotch_sources ${REMOVENONSOURCEDOTCFILES})

add_library(myScotch STATIC GLOBAL)

# generate scotch.h stuff
add_executable(dummysize src/libscotch/dummysizes.c)
target_compile_options(
  dummysize PRIVATE
  "-DSCOTCH_VERSION=6"
  "-DSCOTCH_RELEASE=0"
  "-DSCOTCH_PATCHLEVEL=4"
  )

set(outParserYYFile ${scotch_out_dir}/src/libscotch/parser_yy.c)
set(outParserLLFile ${scotch_out_dir}/src/libscotch/parser_ll.c)
set(outParserLYFile ${scotch_out_dir}/src/libscotch/parser_ly.h)
set(outScotchHFile ${scotch_out_dir}/include/scotch.h)

#FIND_PACKAGE(BISON)
#FIND_PACKAGE(FLEX)

# do without bison/flex for now (does not work yet)
set(BISON_FOUND FALSE)
set(FLEX_FOUND FALSE)

file(MAKE_DIRECTORY ${scotch_out_dir}/include)
file(MAKE_DIRECTORY ${scotch_out_dir}/src)
file(MAKE_DIRECTORY ${scotch_out_dir}/src/libscotch)

if(BISON_FOUND AND FLEX_FOUND)
  SET(linux_yacc_flags "--defines=${outParserLYFile} -pscotchyy -y -b y -v")

  message("Bison/yacc and Flex/lex found ! use them for file generation")
  BISON_TARGET(YaccParser ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_yy.y ${outParserYYFile} COMPILE_FLAGS ${linux_yacc_flags})
  #ADD_CUSTOM_COMMAND(
  #  OUTPUT ${outParserYYFile}
  #  COMMAND ${BISON_EXECUTABLE}
  #    ${linux_yacc_flags} -d -v
  #    ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_yy.y
  #  COMMENT "YACC creating parser_yy.c and parser_ly.h"
  #  )
  #ADD_CUSTOM_COMMAND(
  #  OUTPUT ${outParserYYFile}
  #  COMMAND ${BISON_EXECUTABLE}
  #    --output=${outParserYYFile}
  #    --defines=${outParserLYFile}
  #    ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_yy.y
  #  COMMENT "YACC creating parser_yy.c and parser_ly.h"
  #  )

  SET(linux_lex_flags "-Pscotchyy -olex.yy.c")
  FLEX_TARGET(LexParser ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_ll.l ${outParserLLFile} COMPILE_FLAGS ${linux_lex_flags})
  #ADD_CUSTOM_COMMAND(
  #  OUTPUT ${outParserLLFile}
  #  COMMAND ${FLEX_EXECUTABLE}
  #    ${linux_lex_flags}
  #    #--outfile=${outParserLLFile}
  #    ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_ll.l
  #  COMMENT "LEX creating parser_ll.c"
  #  )

else(BISON_FOUND AND FLEX_FOUND)
  # if lex and yacc are not found we use the files provided by scotch in the last_resort dir
  SET(parserFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/last_resort/parser_yy.c)
  SET(parserLLFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/last_resort/parser_ll.c)
  SET(parserHeaderFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/last_resort/parser_ly.h)

  # first time cmake is run we need to copy the files so that it
  # finds them when adding them to the myScotch target...
  file(COPY ${parserFile} DESTINATION ${scotch_out_dir}/src/libscotch)
  file(COPY ${parserHeaderFile} DESTINATION ${scotch_out_dir}/src/libscotch)
  file(COPY ${parserLLFile} DESTINATION ${scotch_out_dir}/src/libscotch)
  add_custom_command(
    OUTPUT ${outParserYYFile}
    COMMAND ${CMAKE_COMMAND} -E copy ${parserFile} ${outParserYYFile}
    DEPENDS ${parserFile}
    COMMENT "copy parser_yy.c from libscotch/last_resort/"
    )
  add_custom_target(copyParserYYDone
    DEPENDS ${outParserYYFile}
    COMMENT "check if re-copying parser_yy.c is required")
  add_custom_command(
    OUTPUT ${outParserLYFile}
    COMMAND ${CMAKE_COMMAND} -E copy ${parserHeaderFile} ${outParserLYFile}
    DEPENDS ${parserHeaderFile}
    COMMENT "copy parser_ly.h from libscotch/last_resort/"
    )
  add_custom_target(copyParserLYDone
    DEPENDS ${outParserLYFile}
    COMMENT "check if re-copying parser_ly.h is required")
  add_custom_command(
    OUTPUT ${outParserLLFile}
    COMMAND ${CMAKE_COMMAND} -E copy ${parserLLFile} ${outParserLLFile}
    DEPENDS ${parserLLFile}
    COMMENT "copy parser_ll.c from libscotch/last_resort/"
    )
  add_custom_target(copyParserLLDone
    DEPENDS ${outParserLLFile}
    COMMENT "check if re-copying parser_ll.c is required")
  add_dependencies(myScotch copyParserYYDone copyParserLYDone copyParserLLDone)
endif(BISON_FOUND AND FLEX_FOUND)

list(APPEND myScotch_sources ${outParserYYFile} ${outParserLLFile})

include_directories(${scotch_out_dir}/include)
include_directories(src/libscotch)

set_property(TARGET myScotch PROPERTY SOURCES ${myScotch_headers} ${myScotch_sources})

# these compile flags come from the scotch make.inc files...
set(scotch_compile_flags "-DHAVE_STDINT_H;-DHAVE_UINT_T;-DCOMMON_RANDOM_FIXED_SEED;-DSCOTCH_RENAME;-Drestrict=__restrict")
if(UNIX OR APPLE)
list(APPEND scotch_compile_flags "-DCOMMON_FILE_COMPRESS_GZ;-DCOMMON_PTHREAD")
if(APPLE)
  list(APPEND scotch_compile_flags "-DCOMMON_TIMING_OLD;-DCOMMON_PTHREAD_BARRIER")
endif(APPLE)
else(UNIX OR APPLE)
list(APPEND scotch_compile_flags "-DCOMMON_WINDOWS;-DYY_NO_UNISTD_H")
endif(UNIX OR APPLE)

 target_compile_options(
   myScotch PRIVATE
   ${scotch_compile_flags})

target_compile_options(
  myScotch PRIVATE
  ${scotch_compile_flags})
target_compile_options(
  dummysize PRIVATE
  ${scotch_compile_flags})

set(inScotchHFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/library.h)

add_custom_command(
  OUTPUT ${outScotchHFile}
  COMMAND dummysize ${inScotchHFile} ${outScotchHFile}
  DEPENDS dummysize ${inScotchHFile}
  COMMENT "run dummysize to generate scotch.h"
  )

add_custom_target(dummysize_run
  DEPENDS ${outScotchHFile}
  COMMENT "check if re-generation of scotch.h is required")

add_dependencies(myScotch dummysize_run)

install(TARGETS myScotch DESTINATION lib EXPORT openMVG-targets)
target_include_directories(myScotch INTERFACE
    $<BUILD_INTERFACE:${scotch_out_dir}/include>
)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${outParserYYFile};${outParserLLFile};${outParserLYFile};${outScotchHFile}")
