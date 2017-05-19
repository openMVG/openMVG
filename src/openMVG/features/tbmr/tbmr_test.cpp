// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/tbmr/tbmr.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/features/feature.hpp"

#include "testing/testing.h"


using namespace openMVG;
using namespace image;
using namespace features;
using namespace openMVG::features::tbmr;

TEST( TBMR , Extraction )
{
  Image<unsigned char> image(600, 400, true, 0);

  // A TBMR will be detected with 3 monotonic color change.
  image.block<300,400>(50,100).setConstant(50);
  image.block<200,300>(100,150).setConstant(127);
  image.block<100,200>(400/2-50, 600/2-100).setConstant(255);

  WriteImage( "in.png" , image );

  std::vector<features::AffinePointFeature> feats_tbmr;
  // Bright TBMR
  tbmr::Extract_tbmr (image, feats_tbmr, std::less<uint8_t> (), 1, 1);
  EXPECT_EQ(1, feats_tbmr.size());
  EXPECT_EQ(300, std::ceil(feats_tbmr[0].x()));
  EXPECT_EQ(200, std::ceil(feats_tbmr[0].y()));

  feats_tbmr.clear();
  // DARK TBMR
  tbmr::Extract_tbmr (image, feats_tbmr, std::greater<uint8_t> (), 1, 1);
  EXPECT_EQ(1, feats_tbmr.size());
  EXPECT_EQ(300, std::ceil(feats_tbmr[0].x()));
  EXPECT_EQ(200, std::ceil(feats_tbmr[0].y()));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
