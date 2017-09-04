##
# SCOTCH library for hypergraph partitioning
##
#URL http://gforge.inria.fr/frs/download.php/zip/9109/scotch-6.0.4.zip

SET(scotch_out_dir ${PROJECT_BINARY_DIR}/third_party/scotch_6.0.4)

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

# generate scotch.h stuff
add_executable(dummysize src/libscotch/dummysizes.c)
target_compile_options(
  dummysize PRIVATE
  "-DSCOTCH_VERSION=6"
  "-DSCOTCH_RELEASE=0"
  "-DSCOTCH_PATCHLEVEL=4"
  )

file(MAKE_DIRECTORY ${scotch_out_dir}/include)
file(MAKE_DIRECTORY ${scotch_out_dir}/src)
file(MAKE_DIRECTORY ${scotch_out_dir}/src/libscotch)

# parser stuff (generated with yacc and lex)
SET(USE_YACC FALSE)
if(USE_YACC)
  SET(YaccOutputFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_yy.c)
  SET(YaccOutputHeaderFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_ly.h)

  FIND_PACKAGE(BISON REQUIRED)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${YaccOutputFile}
    COMMAND ${BISON_EXECUTABLE}
      --output=${YaccOutputFile}
      --defines=${YaccOutputHeaderFile}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_yy.y
    COMMENT "YACC creating parser.c"
    )

  SET(LexOutputFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_ll.c)
  FIND_PACKAGE(FLEX REQUIRED)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${LexOutputFile}
    COMMAND ${FLEX_EXECUTABLE}
      --outfile=${LexOutputFile}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/parser_ll.l
    COMMENT "LEX creating parser.c"
    )

  list(APPEND myScotch_sources ${YaccOutputFile})
  list(APPEND myScotch_sources ${LexOutputFile})
else(USE_YACC)
  SET(parserFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/last_resort/parser_yy.c)
  SET(parserFileCopy ${scotch_out_dir}/src/libscotch/parser_yy.c)
  SET(parserLLFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/last_resort/parser_ll.c)
  SET(parserLLFileCopy ${scotch_out_dir}/src/libscotch/parser_ll.c)
  SET(parserHeaderFile ${CMAKE_CURRENT_SOURCE_DIR}/src/libscotch/last_resort/parser_ly.h)
  SET(parserHeaderFileCopy ${scotch_out_dir}/src/libscotch/parser_ly.h)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${parserFileCopy}
    COMMAND cp
      ${parserFile}
      ${parserFileCopy}
    COMMENT "copying parser_yy.c"
    )
  ADD_CUSTOM_COMMAND(
    OUTPUT ${parserHeaderFileCopy}
    COMMAND cp
      ${parserHeaderFile}
      ${parserHeaderFileCopy}
    COMMENT "copying parser_ly.h"
    )
  ADD_CUSTOM_COMMAND(
    OUTPUT ${parserLLFileCopy}
    COMMAND cp
      ${parserLLFile}
      ${parserLLFileCopy}
    COMMENT "copying parser_ly.h"
    )
  list(APPEND myScotch_sources ${parserFileCopy} ${parserLLFileCopy})
  list(APPEND myScotch_headers ${parserHeaderFileCopy})
endif(USE_YACC)

include_directories(${scotch_out_dir}/include)
include_directories(src/libscotch)

file(WRITE ${scotch_out_dir}/include/scotch.h "#writing empty scotch.h for cmake")
list(APPEND myScotch_headers ${scotch_out_dir}/include/scotch.h)

add_library(myScotch STATIC GLOBAL)
set_property(TARGET myScotch PROPERTY SOURCES ${myScotch_headers} ${myScotch_sources})
target_compile_options(
  myScotch PRIVATE
  "-DSCOTCH_PTHREAD")

add_custom_target(dummysize_run
  COMMAND dummysize src/libscotch/library.h ${scotch_out_dir}/include/scotch.h
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "run dummysize to generate scotch.h"
  SOURCES src/libscotch/dummysizes.c)

add_dependencies(myScotch dummysize_run)

install(TARGETS myScotch DESTINATION lib EXPORT openMVG-targets)
