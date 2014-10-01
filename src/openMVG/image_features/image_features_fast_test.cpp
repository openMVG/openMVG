// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/image_features/detector_fast.hpp"

#include "testing/testing.h"

#include <iostream>
using namespace std;

using namespace openMVG;
using namespace openMVG::image_features;

TEST(FAST, Test)
{
  //-- Gray(unsigned char) Image creation
  Image<unsigned char> imaGray(20, 20, true, 127);
  imaGray(12,10) = 255;
  
  // Try to detect corner
  std::vector<SIOPointFeature> feats;
  FastCornerDetector fastCornerDetector;
  fastCornerDetector.detect(imaGray, feats);
    
  // 1 detection must be found
  EXPECT_EQ(1, feats.size());
  EXPECT_EQ(10, feats[0].x());
  EXPECT_EQ(12, feats[0].y());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
