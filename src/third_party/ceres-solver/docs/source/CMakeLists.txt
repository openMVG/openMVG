# HTML output directory
set(SPHINX_HTML_DIR "${Ceres_BINARY_DIR}/docs/html")

# Install documentation
install(DIRECTORY ${SPHINX_HTML_DIR}
        DESTINATION ${CMAKE_INSTALL_DOCDIR}
        COMPONENT Doc
        PATTERN "${SPHINX_HTML_DIR}/*")

# Find python
find_package(Python REQUIRED COMPONENTS Interpreter)

# Building using 'make_docs.py' python script
add_custom_target(ceres_docs ALL
                  $<TARGET_FILE:Python::Interpreter>
                  "${Ceres_SOURCE_DIR}/scripts/make_docs.py"
                  "${Ceres_SOURCE_DIR}"
                  "${Ceres_BINARY_DIR}/docs"
                  "${Sphinx_BUILD_EXECUTABLE}"
                  USES_TERMINAL
                  COMMENT "Building HTML documentation with Sphinx")
