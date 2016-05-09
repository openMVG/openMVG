find_package(Sphinx REQUIRED)

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_BINARY_DIR}/docs/html")

# Install documentation
install(DIRECTORY ${SPHINX_HTML_DIR}
        DESTINATION share/doc/ceres
        COMPONENT Doc
        PATTERN "${SPHINX_HTML_DIR}/*")

# Building using 'make_docs.py' python script
add_custom_target(ceres_docs ALL
                  python
                  "${CMAKE_SOURCE_DIR}/scripts/make_docs.py"
                  "${CMAKE_SOURCE_DIR}"
                  "${CMAKE_BINARY_DIR}/docs"
                  "${SPHINX_EXECUTABLE}"
                  COMMENT "Building HTML documentation with Sphinx")
