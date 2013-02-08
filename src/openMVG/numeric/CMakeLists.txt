ADD_LIBRARY(openMVG_numeric
			numeric.cpp
			numeric.h)
			
UNIT_TEST(openMVG numeric openMVG_numeric)
UNIT_TEST(openMVG poly openMVG_numeric)
UNIT_TEST(openMVG lm openMVG_numeric)
