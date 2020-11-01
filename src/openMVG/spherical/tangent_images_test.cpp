// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Marc Eder.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/spherical/tangent_images.hpp"
#include "openMVG/image/image_io.hpp"
#include "testing/testing.h"

#include "nonFree/sift/SIFT_describer_io.hpp"
#include "openMVG/features/image_describer.hpp"

#include <sstream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::spherical;
using namespace openMVG::features;

/* Test and demo to create tangent images */
TEST(Spherical, EquirectToTangent) {
  // Load the test equirectangular image
  Image<RGBColor> image;
  const std::string png_filename =
      std::string(THIS_SOURCE_DIR) + "/earthmap4k.jpg";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  // Instantiate a TangentImage object that defines the relationship between the
  // dimension of the equirectangular image and those of the tangent images we
  // want to create
  TangentImages tangent_images(0, 9, image.Height(), image.Width());

  // For fun, this is the FOV of each tangent image we're going to create
  std::cout << "FOV: " << tangent_images.FOV() << " degrees" << std::endl;

  // Create the tangent images
  std::vector<Image<RGBColor>> t_images;
  std::vector<Image<unsigned char>> t_mask;
  tangent_images.CreateTangentImages(image, t_images, &t_mask);

  // Now write the create tangent images to file
  for (size_t i = 0; i < t_images.size(); i++) {
    const std::string out_filename =
        ("test_tangent images_" + std::to_string(i) + ".png");
    EXPECT_TRUE(WriteImage(out_filename.c_str(), t_images[i]));
    const std::string mask_filename =
        ("test_tangent images_" + std::to_string(i) + "_mask.png");
    EXPECT_TRUE(WriteImage(mask_filename.c_str(), t_mask[i]));
  }
}

/* Test and demo to detect features on tangent images */
TEST(Spherical, DescribeFeatures) {
  // Load the test equirectangular image
  Image<unsigned char> imageGray;
  const std::string png_filename =
      std::string(THIS_SOURCE_DIR) + "/earthmap4k.jpg";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &imageGray));

  // Instantiate a TangentImage object that defines the relationship between the
  // dimension of the equirectangular image and those of the tangent images we
  // want to create
  TangentImages tangent_images(0, 9, imageGray.Height(), imageGray.Width());

  // Create SIFT image describer
  std::unique_ptr<Image_describer> image_describer;
  image_describer.reset(
      new SIFT_Image_describer(SIFT_Image_describer::Params(), true));

  // Compute features on tangent images
  auto regions = ComputeFeaturesOnTangentImages(tangent_images, image_describer,
                                                imageGray);

  // Write the tangent images to file
  if (regions &&
      !image_describer->Save(regions.get(), "feat.txt", "desc.txt")) {
    std::cerr << "Cannot save regions" << std::endl;
  }
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */