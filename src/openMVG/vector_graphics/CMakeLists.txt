add_executable(main_svgSample main.cpp)
set_property(TARGET main_svgSample PROPERTY FOLDER OpenMVG/Samples)

add_library(openMVG_svg INTERFACE)
set_property(TARGET openMVG_svg APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES "$<INSTALL_INTERFACE:include>")

install(TARGETS openMVG_svg DESTINATION lib EXPORT openMVG-targets)
