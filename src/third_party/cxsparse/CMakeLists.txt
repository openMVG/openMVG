
# Author: TheFrenchLeaf https://github.com/TheFrenchLeaf

project(CXSparse C)

file(GLOB CXSparse_SRCS "Source/*.c")
set_source_files_properties(${CXSparse_SRCS} PROPERTIES LANGUAGE C)

add_library(openMVG_cxsparse STATIC ${CXSparse_SRCS})
target_compile_definitions(openMVG_cxsparse PUBLIC NCOMPLEX)
target_include_directories(openMVG_cxsparse PRIVATE ./ ./Include)
if(UNIX)
  target_link_libraries(openMVG_cxsparse m)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif(UNIX)

install(TARGETS openMVG_cxsparse DESTINATION lib EXPORT openMVG-targets)
