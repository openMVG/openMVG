set(PKG_DESC "Fast Library for Approximate Nearest Neighbors")
set(pkg_conf_file ${CMAKE_CURRENT_BINARY_DIR}/flann.pc)
configure_file(flann.pc.in ${pkg_conf_file} @ONLY)
install(FILES ${pkg_conf_file}
    DESTINATION ${FLANN_LIB_INSTALL_DIR}/pkgconfig/ COMPONENT pkgconfig)

