SET(PACKAGE_NAME ${PROJECT_NAME})
SET(PACKAGE_VERSION ${PROJECT_VERSION})
SET(abs_top_srcdir ${PROJECT_SOURCE_DIR})
SET(abs_top_builddir ${PROJECT_BINARY_DIR})

SET(LEMON_DOC_SOURCE_BROWSER "NO" CACHE STRING "Include source into the doc (YES/NO).")
SET(LEMON_DOC_USE_MATHJAX "NO" CACHE STRING "Use MathJax to display math formulae (YES/NO).")
SET(LEMON_DOC_MATHJAX_RELPATH "http://www.mathjax.org/mathjax" CACHE STRING "MathJax library location.")

SET(LEMON_DOC_LIBSTDC++_URL
  "http://gcc.gnu.org/onlinedocs/gcc-4.7.3/libstdc++/api"
  CACHE STRING "GCC libstdc++ doxygen doc url.")


CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
  ${PROJECT_BINARY_DIR}/doc/Doxyfile
  @ONLY
)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/doc/mainpage.dox.in
  ${PROJECT_BINARY_DIR}/doc/mainpage.dox
  @ONLY
)

# Copy doc from source (if exists)
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/html AND 
    NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  MESSAGE(STATUS "Copy doc from source tree")
  EXECUTE_PROCESS(
    COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/html ${CMAKE_CURRENT_BINARY_DIR}/html
    )
ENDIF()

IF(DOXYGEN_EXECUTABLE AND PYTHONINTERP_FOUND AND GHOSTSCRIPT_EXECUTABLE)
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/)
  SET(GHOSTSCRIPT_OPTIONS -dNOPAUSE -dBATCH -q -dEPSCrop -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sDEVICE=pngalpha)
  ADD_CUSTOM_TARGET(html
    COMMAND ${CMAKE_COMMAND} -E remove_directory gen-images
    COMMAND ${CMAKE_COMMAND} -E make_directory gen-images
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r20 -sOutputFile=gen-images/grid_graph.png ${CMAKE_CURRENT_SOURCE_DIR}/images/grid_graph.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32 -sOutputFile=gen-images/adaptors2.png ${CMAKE_CURRENT_SOURCE_DIR}/images/adaptors2.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32 -sOutputFile=gen-images/connected_components.png ${CMAKE_CURRENT_SOURCE_DIR}/images/connected_components.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32 -sOutputFile=gen-images/strongly_connected_components.png ${CMAKE_CURRENT_SOURCE_DIR}/images/strongly_connected_components.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32 -sOutputFile=gen-images/node_biconnected_components.png ${CMAKE_CURRENT_SOURCE_DIR}/images/node_biconnected_components.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32 -sOutputFile=gen-images/edge_biconnected_components.png ${CMAKE_CURRENT_SOURCE_DIR}/images/edge_biconnected_components.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32 -sOutputFile=gen-images/bipartite_partitions.png ${CMAKE_CURRENT_SOURCE_DIR}/images/bipartite_partitions.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r24 -sOutputFile=gen-images/matching.png ${CMAKE_CURRENT_SOURCE_DIR}/images/matching.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r24 -sOutputFile=gen-images/bipartite_matching.png ${CMAKE_CURRENT_SOURCE_DIR}/images/bipartite_matching.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r40 -sOutputFile=gen-images/planar.png ${CMAKE_CURRENT_SOURCE_DIR}/images/planar.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r24 -sOutputFile=gen-images/tsp.png ${CMAKE_CURRENT_SOURCE_DIR}/images/tsp.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8 -sOutputFile=gen-images/nodeshape_0.png ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_0.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8 -sOutputFile=gen-images/nodeshape_1.png ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_1.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8 -sOutputFile=gen-images/nodeshape_2.png ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_2.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8 -sOutputFile=gen-images/nodeshape_3.png ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_3.eps
    COMMAND ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8 -sOutputFile=gen-images/nodeshape_4.png ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_4.eps
    COMMAND ${CMAKE_COMMAND} -E remove_directory html
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  SET_TARGET_PROPERTIES(html PROPERTIES PROJECT_LABEL BUILD_DOC)

  IF(UNIX)
    INSTALL(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
      DESTINATION share/doc/lemon/html
      COMPONENT html_documentation
    )
  ELSEIF(WIN32)
    INSTALL(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
      DESTINATION doc
      COMPONENT html_documentation
    )
  ENDIF()

ENDIF()

IF(WGET_FOUND)
ADD_CUSTOM_TARGET(update-external-tags
  COMMAND ${WGET_EXECUTABLE} -N ${LEMON_DOC_LIBSTDC++_URL}/libstdc++.tag
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()
