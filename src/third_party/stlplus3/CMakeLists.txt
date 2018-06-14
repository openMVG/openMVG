project(stlplus C CXX)

file(GLOB LIBSLTPLUS_HPP_FILESYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/filesystemSimplified/*.hpp" )
file(GLOB LIBSLTPLUS_CPP_FILESYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/filesystemSimplified/*.cpp" )

list(APPEND LIBSLTPLUS_SRCS
  ${LIBSLTPLUS_HPP_FILESYSTEM} ${LIBSLTPLUS_CPP_FILESYSTEM})

add_library(openMVG_stlplus ${LIBSLTPLUS_SRCS})
if(STLPLUS_INCLUDE_INSTALL_DIR)
  target_include_directories(openMVG_stlplus PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/openMVG>
    $<INSTALL_INTERFACE:${STLPLUS_INCLUDE_INSTALL_DIR}>/filesystemSimplified)
else()
  target_include_directories(openMVG_stlplus PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/stlplus>)
endif()
set_property(TARGET openMVG_stlplus PROPERTY FOLDER OpenMVG/3rdParty)

install(
  TARGETS openMVG_stlplus
  DESTINATION lib
  EXPORT openMVG-targets
)

if(STLPLUS_INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR ${STLPLUS_INCLUDE_INSTALL_DIR})
else()
  set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/stlplus")
endif()

install(
  DIRECTORY ./filesystemSimplified
  DESTINATION ${INCLUDE_INSTALL_DIR}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp"
)
