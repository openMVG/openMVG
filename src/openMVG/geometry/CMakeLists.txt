
UNIT_TEST(openMVG rigid_transformation3D_srt "openMVG_numeric")

UNIT_TEST(openMVG half_space_intersection "openMVG_numeric;${openMVG_linearProgramming}")

UNIT_TEST(openMVG frustum_intersection "openMVG_numeric;${openMVG_linearProgramming};openMVG_multiview_test_data;openMVG_multiview")
