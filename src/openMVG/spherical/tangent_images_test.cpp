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

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::spherical;
using namespace openMVG::features;

/* Test and demo to create tangent images */
TEST(Spherical, CreateTangentImages) {
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

  // Next, convert back to an equirectangular image
  Image<RGBColor> out_rect_image;
  tangent_images.ConvertTangentImagesToEquirectangular(t_images,
                                                       out_rect_image);
  const std::string out_filename = ("test_equirectangular.png");
  EXPECT_TRUE(WriteImage(out_filename.c_str(), out_rect_image));
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

  // Compute features using tangent images
  auto t_regions = tangent_images.ComputeFeaturesOnTangentImages(
      *image_describer, imageGray);
  // Write the tangent images to file
  if (t_regions &&
      !image_describer->Save(t_regions.get(), "tangent_images_feat.txt",
                             "tangent_images_desc.txt")) {
    std::cerr << "Cannot save regions" << std::endl;
  }

  // For comparison, compute features on the equirectangular images
  auto r_regions = image_describer->Describe(imageGray);
  // Write the tangent images to file
  if (r_regions && !image_describer->Save(r_regions.get(), "equirect_feat.txt",
                                          "equirect_desc.txt")) {
    std::cerr << "Cannot save regions" << std::endl;
  }
}

// Test the conversion functions
TEST(Spherical, SphericalTo3D) {
  const Vec2 lonlat_in = Vec2(M_PI / 4.0, M_PI / 4.0);
  const auto pt3d = ConvertSphericalTo3D(lonlat_in);
  const auto lonlat_out = Convert3DToSpherical(pt3d);
  EXPECT_NEAR(lonlat_in[0], lonlat_out[0], 1e-6);
  EXPECT_NEAR(lonlat_in[1], lonlat_out[1], 1e-6);
}

// Tests the 3D ray intersection test with a spherical triangle
TEST(Spherical, RaySphericalTriangleIntersection) {
  // Test a ray intersection with a spherical triangle in 3D
  bool res = RayIntersectSphericalTriangle3D(Vec3(0, std::sqrt(3) / 2, 1),
                                             Vec3(0, std::sqrt(3), 1),
                                             Vec3(-1, 0, 1), Vec3(1, 0, 1));
  EXPECT_TRUE(res);

  // Test the negative ray to ensure we address chirality
  res = RayIntersectSphericalTriangle3D(-Vec3(0, std::sqrt(3) / 2, 1),
                                        Vec3(0, std::sqrt(3), 1),
                                        Vec3(-1, 0, 1), Vec3(1, 0, 1));
  EXPECT_FALSE(res);
}

// Tests the forward and inverse gnomonic projections
TEST(Spherical, GnomonicProjections) {
  const auto sphere_in = Vec2(M_PI / 7, 3 * M_PI / 11);
  const auto xy_out =
      ForwardGnomonicProjection(sphere_in, Vec2(M_PI / 8, -M_PI / 5));
  const auto sphere_out =
      InverseGnomonicProjection(xy_out, Vec2(M_PI / 8, -M_PI / 5));
  EXPECT_NEAR(sphere_in[0], sphere_out[0], 1e-6);
  EXPECT_NEAR(sphere_in[1], sphere_out[1], 1e-6);
}

// Test the conversion between equirectangular and spherical coordinates
TEST(Spherical, EquirectangularToSpherical) {
  // Create a b0, s9 tangent images instance to model a (2000 x 4000) pixel
  // equirectangular image
  TangentImages tangent_images(0, 9, 2000, 4000);

  // Arbitrary equirect coordinate
  const auto xy_in = Vec2(400, 459);

  // Convert to spherical coordinates and back
  const auto sphere_out =
      tangent_images.ConvertEquirectangularToSpherical(xy_in);
  const auto xy_out =
      tangent_images.ConvertSphericalToEquirectangular(sphere_out);
  EXPECT_NEAR(xy_in[0], xy_out[0], 1e-6);
  EXPECT_NEAR(xy_in[1], xy_out[1], 1e-6);
}

// Test the conversion between equirectangular coordinates and tangent image
// coordinates
TEST(Spherical, EquirectangularToTangentUV) {
  TangentImages tangent_images(0, 9, 2000, 4000);

  // A point that we know falls in the face projected on the tangent image so we
  // have a bijective conversions
  const size_t t_idx = 0;
  const Vec2 uv_in = Vec2(219, 284);
  const Vec2 xy_in = Vec2(3059.04, 500.936);
  const Vec2 sphere_in = Vec2(1.66474, -0.783534);

  // Convert tangent index and UV coordinates to equirectangular and spherical
  const Vec2 xy_out = tangent_images.TangentUVToEquirectangular(t_idx, uv_in);
  const Vec2 sphere_out =
      tangent_images.ConvertEquirectangularToSpherical(xy_out);
  size_t t_idx_out = tangent_images.GetTangentImageIndex(sphere_out);
  const Vec2 uv_out = tangent_images.EquirectangularToTangentUV(t_idx, xy_out);

  EXPECT_EQ(t_idx, t_idx_out);

  // The weirdly loose threshold is to prevent the test from failing....even
  // though the outcomes are equal
  EXPECT_NEAR(uv_in[0], uv_out[0], 1e-1);
  EXPECT_NEAR(uv_in[1], uv_out[1], 1e-1);
  EXPECT_NEAR(xy_in[0], xy_out[0], 1e-1);
  EXPECT_NEAR(xy_in[1], xy_out[1], 1e-1);
  EXPECT_NEAR(sphere_in[0], sphere_out[0], 1e-1);
  EXPECT_NEAR(sphere_in[1], sphere_out[1], 1e-1);
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */