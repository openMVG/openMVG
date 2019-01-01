project(openMVG_fast)
file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
add_library(openMVG_fast STATIC ${C_FILES})
set_property(TARGET openMVG_fast PROPERTY FOLDER OpenMVG/3rdParty)

install(TARGETS openMVG_fast DESTINATION lib EXPORT openMVG-targets)
