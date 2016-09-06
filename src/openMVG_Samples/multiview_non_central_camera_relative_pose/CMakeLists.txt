
find_package(opengv QUIET)
if (opengv_FOUND)

  message("OPENMVG-OpenGV sample will be compiled")
  # Unit tests
  UNIT_TEST(openMVG solver_non_central_kernel "openMVG_multiview;opengv")

endif (opengv_FOUND)
