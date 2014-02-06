cmake_minimum_required(VERSION 2.6)

project(JsonBox)

# glob all sources
file(GLOB_RECURSE JSONBOX_SOURCES "src" "src/*.cpp")
file(GLOB_RECURSE JSONBOX_HEADERS "include" "include/*.h")

include_directories("include")

# setup static archive
set(CMAKE_DEBUG_POSTFIX "_d")
add_library(JsonBox STATIC ${JSONBOX_SOURCES} ${JSONBOX_HEADERS})

# # install
# install(TARGETS JsonBox ARCHIVE DESTINATION "lib" COMPONENT lib)
# install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/" DESTINATION "include" COMPONENT dev)

# # examples
# add_executable(example1 "${CMAKE_SOURCE_DIR}/examples/main.cpp")
# add_dependencies(example1 JsonBox)
# target_link_libraries(example1 JsonBox)
