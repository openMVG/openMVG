
# Author: TheFrenchLeaf https://github.com/TheFrenchLeaf

project(CXSparse C)

file(GLOB CXSparse_SRCS "Source/*.c")
set_source_files_properties(${CXSparse_SRCS} PROPERTIES LANGUAGE C)

add_library(cxsparse STATIC ${CXSparse_SRCS})
target_compile_definitions(cxsparse PUBLIC NCOMPLEX)
target_include_directories(cxsparse PRIVATE ./ ./Include)
if(UNIX)
  target_link_libraries(cxsparse m)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif(UNIX)
