// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/hierarchical_gaussian_scale_space.hpp"
#include "openMVG/features/sift/sift_DescriptorExtractor.hpp"
#include "openMVG/features/sift/sift_keypoint.hpp"
#include "openMVG/features/sift/sift_KeypointExtractor.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/system/timer.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include "testing/testing.h"

#include <sstream>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::features::sift;
using namespace openMVG::system;

TEST( GaussianScaleSpace , OctaveGeneration )
{
  Image<unsigned char> in;

  const std::string png_filename = std::string( THIS_SOURCE_DIR )
    + "/../../../openMVG_Samples/imageData/StanfordMobileVisualSearch/Ace_0.png";
  EXPECT_TRUE( ReadImage( png_filename.c_str(), &in ) );

  Image<float> image;
  image = in.GetMat().cast<float>();

  HierarchicalGaussianScaleSpace octave_gen(6, 3, GaussianScaleSpaceParams(1.6f, 1.0f, 0.5f));
  Octave octave;

  octave_gen.SetImage( image );

  std::cerr << "Octave computation started" << std::endl;
  uint8_t octave_id = 0;
  while (octave_gen.NextOctave( octave ))
  {
    std::cerr << "Computed octave : " << std::to_string(octave_id) << std::endl;
    for (int i = 0; i < octave.slices.size(); ++i )
    {
      std::stringstream str;
      str << "gaussian_octave_" << std::to_string(octave_id) << "_" << i << ".png";
      EXPECT_TRUE( WriteImage( str.str().c_str() , Image<unsigned char>( octave.slices[i].cast<unsigned char>()  ) ) );
    }
    ++octave_id;
  }
}

TEST( Sift_Keypoint , DetectionAndDescription )
{
  Image<unsigned char> in;

  const std::string png_filename = std::string( THIS_SOURCE_DIR )
    + "/../../../openMVG_Samples/imageData/StanfordMobileVisualSearch/Ace_0.png";
  EXPECT_TRUE( ReadImage( png_filename.c_str(), &in ) );

  openMVG::system::Timer timer;

  const int supplementary_images = 3;
  // => in order to ensure each gaussian slice is used in the process 3 extra images are required:
  // +1 for dog computation
  // +2 for 3d discrete extrema definition
  HierarchicalGaussianScaleSpace octave_gen(6, 3, GaussianScaleSpaceParams(1.6f, 1.0f, 0.5f, supplementary_images));
  Octave octave;

  // Convert to float in range [0;1]
  const image::Image<float> image(in.GetMat().cast<float>()/255.0f);
  octave_gen.SetImage( image );

  std::vector<Keypoint> keypoints;
  keypoints.reserve(5000);
  std::cerr << "Octave computation started" << std::endl;
  uint8_t octave_id = 0;
  while (octave_gen.NextOctave( octave ))
  {
    std::cerr << "Computed octave : " << std::to_string(octave_id) << std::endl;
    std::vector< Keypoint > keys;
    SIFT_KeypointExtractor keypointDetector(0.04f / octave_gen.NbSlice(), 10.f, 5);
    keypointDetector(octave, keys);
    Sift_DescriptorExtractor descriptorExtractor;
    descriptorExtractor(octave, keys);

    // Concatenante the keypoints
    std::move(keys.begin(), keys.end(), std::back_inserter(keypoints));
    ++octave_id;
  }
  keypoints.shrink_to_fit();
  EXPECT_TRUE(keypoints.size() > 0);
  std::cout << "\n#Keypoints: " << keypoints.size() << std::endl;
  std::cout << "OpenMVG sift: " << timer.elapsedMs() << "ms" << std::endl;

  //--
  // Export found keypoints as a SVG surface
  // SIFT features drawn as circle and one line for the orientation
  //--
  using namespace svg;
  svgDrawer svgStream( image.Width(), image.Height());
  svgStream.drawImage(png_filename, image.Width(), image.Height());
  for (size_t i = 0; i < keypoints.size(); ++i) {
    const Keypoint & key = keypoints[i];
    svgStream.drawCircle(key.x, key.y, key.sigma, svgStyle().stroke("yellow", 2.0));
    // Orientation
    svgStream.drawLine(key.x, key.y,
      key.x + cos(key.sigma) * key.sigma,
      key.y + sin(key.sigma) * key.sigma, svgStyle().stroke("green", 2.0));
  }
  std::string out_filename = "Sift_Features.svg";
  std::ofstream svgFile( out_filename.c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();
}

/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */
