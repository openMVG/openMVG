// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/image/image_io.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;

static const std::string png_filename = std::string( THIS_SOURCE_DIR )
  + "/../../../openMVG_Samples/imageData/StanfordMobileVisualSearch/Ace_0.png";

TEST( AKAZE , EmptyImage )
{
  Image<unsigned char> image_in;
  AKAZE akaze_extractor(image_in, AKAZE::Params());
  akaze_extractor.Compute_AKAZEScaleSpace();
  std::vector<AKAZEKeypoint> keypoints;
  akaze_extractor.Feature_Detection(keypoints);
  EXPECT_TRUE(keypoints.empty());
}

TEST( AKAZE , AkazeImageDescriberSurf )
{
  Image<unsigned char> image_in;
  EXPECT_TRUE( ReadImage( png_filename.c_str(), &image_in ) );

  AKAZE_Image_describer_SURF extractor;
  EXPECT_TRUE(extractor.Describe(Image<unsigned char>{})->RegionCount() == 0);
  EXPECT_TRUE(extractor.Describe(image_in)->RegionCount() > 0);
}

TEST( AKAZE , AkazeImageDescriberLiop )
{
  Image<unsigned char> image_in;
  EXPECT_TRUE( ReadImage( png_filename.c_str(), &image_in ) );

  AKAZE_Image_describer_LIOP extractor;
  EXPECT_TRUE(extractor.Describe(Image<unsigned char>{})->RegionCount() == 0);
  EXPECT_TRUE(extractor.Describe(image_in)->RegionCount() > 0);
}

TEST( AKAZE , AKAZE_Image_describer_MLDB )
{
  Image<unsigned char> image_in;
  EXPECT_TRUE( ReadImage( png_filename.c_str(), &image_in ) );

  AKAZE_Image_describer_MLDB extractor;
  EXPECT_TRUE(extractor.Describe(Image<unsigned char>{})->RegionCount() == 0);
  EXPECT_TRUE(extractor.Describe(image_in)->RegionCount() > 0);
}

/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */
