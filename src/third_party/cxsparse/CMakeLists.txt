
# Author: TheFrenchLeaf https://github.com/TheFrenchLeaf

cmake_minimum_required(VERSION 2.6)

project(CXSparse C)

file(GLOB CXSparse_SRCS "Source/*.c")
set_source_files_properties(${CXSparse_SRCS} PROPERTIES LANGUAGE C)

add_library(cxsparse STATIC ${CXSparse_SRCS})
target_compile_definitions(cxsparse PUBLIC NCOMPLEX)
target_include_directories(cxsparse PRIVATE ./ ./Include)
if(UNIX)
  target_link_libraries(cxsparse m)
  set_target_properties(cxsparse PROPERTIES COMPILE_FLAGS -fPIC)
endif(UNIX)

if(CXSPARSE_INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR ${CXSPARSE_INCLUDE_INSTALL_DIR})
else()
  set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/cxsparse")
endif()

install(
  DIRECTORY ./Include
  DESTINATION ${INCLUDE_INSTALL_DIR}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h"
)
install(TARGETS cxsparse DESTINATION lib EXPORT openMVG-targets)

