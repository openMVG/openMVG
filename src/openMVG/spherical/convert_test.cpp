// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/spherical/cubic_image_sampler.hpp"
#include "openMVG/image/image_io.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::cameras;

// Define some global variable that gonna be shared by the various demo sample
static const int cubic_size = 256;
static const std::vector<Image<RGBColor>> images =
{
  Image<RGBColor>(cubic_size, cubic_size, true, BLUE),
  Image<RGBColor>(cubic_size, cubic_size, true, RED),
  Image<RGBColor>(cubic_size, cubic_size, true, GREEN),
  Image<RGBColor>(cubic_size, cubic_size, true, YELLOW),
  Image<RGBColor>(cubic_size, cubic_size, true, CYAN),
  Image<RGBColor>(cubic_size, cubic_size, true, MAGENTA)
};
static const openMVG::cameras::Pinhole_Intrinsic pinhole_camera =
  spherical::ComputeCubicCameraIntrinsics(cubic_size);

static const std::array<openMVG::Mat3,6> cubic_image_rotations = spherical::GetCubicRotations();
// Project each pinhole on a sphere
static const Intrinsic_Spherical spherical_camera(512, 256);

void SphericalToPinholeIndirect
(
  Image<RGBColor> & spherical_image
)
{
  // Compute the projection of the tiles on the pano (indirect)
  // For each pixel of the spherical image, find back which tile pixel project on it
  {
    for (int x = 0; x < spherical_image.Width(); ++x)
    for (int y = 0; y < spherical_image.Height(); ++y)
    {
      // Build a bearing vector
      const Vec3 bearing_vector = spherical_camera(Vec2(x + .5,y + .5));
      // See if the bearing_vector reproject on a pinhole images
      for (const int i_rot : {0, 1, 2, 3, 4, 5}) // For every rotation
      {
        const Mat3 rotation_matrix = cubic_image_rotations[i_rot];
        const Vec2 proj = pinhole_camera.project(rotation_matrix * bearing_vector);
        const Vec3 pinhole_cam_look_dir = rotation_matrix.transpose() * pinhole_camera(Vec2(spherical_image.Width()/2,spherical_image.Height()/2));
        const auto & pinhole_image = images[i_rot];
        if (// Look if the projection is projecting into the image plane
           proj.x() > 0 && proj.x() < pinhole_image.Width() &&
           proj.y() > 0 && proj.y() < pinhole_image.Height() &&
           // Look if the camera and the sampled sphere vector are looking in the same direction
           pinhole_cam_look_dir.dot(bearing_vector) > 0.0
          )
        {
          spherical_image(y, x) = pinhole_image((int)proj.y(), (int)proj.x());
          break;
        }
      }
    }
  }
}

TEST(Ressampling, PinholeToSpherical_indirect_mapping)
{
  Image<RGBColor> spherical_image(512, 256);
  SphericalToPinholeIndirect(spherical_image);
  WriteImage("test_pinhole_to_spherical_indirect.png", spherical_image);
}

TEST(Ressampling, PinholeToSpherical_direct_mapping)
{
  Image<RGBColor> spherical_image(512, 256);

  for (const int i_rot : {0, 1, 2, 3, 4, 5}) // For every rotation
  {
    const auto & pinhole_image = images[i_rot];
    const Mat3 rotation_matrix = cubic_image_rotations[i_rot];
    // Project each pixel
    Mat2X xy_coords(2, static_cast<int>(cubic_size * cubic_size));
    for (int x = 0; x < cubic_size; ++x)
    for (int y = 0; y < cubic_size; ++y)
      xy_coords.col(x + cubic_size * y ) << x + .5, y +.5;
    // Compute bearing vectors
    const Mat3X bearing_vectors = pinhole_camera(xy_coords);
    // Compute rotated bearings
    const Mat3X rotated_bearings = rotation_matrix.transpose() * bearing_vectors;
    for (int it = 0; it < rotated_bearings.cols(); ++it)
    {
      // Project the bearing vector to the sphere
      const Vec2 sphere_proj = spherical_camera.project(rotated_bearings.col(it));
      // and use the corresponding pixel location in the panorama
      const Vec2 xy = xy_coords.col(it);
      if (spherical_image.Contains(sphere_proj.y(), sphere_proj.x()))
      {
        spherical_image(sphere_proj.y(), sphere_proj.x()) =
          pinhole_image((int)xy.y(),(int)xy.x());
      }
    }
  }
  WriteImage("test_pinhole_to_spherical_direct.png", spherical_image);
  // Notes:
  // Since the projection is direct, we project only the pixel from the pinhole to the spherical domain
  // -> It's gonna create hole in the spherical image for top and bottom image
  // The indirect mapping is sampling from destimation to source
  // Direct mapping is mapping from source to destination and so cannot fill all the pixel if the cubic image is small
}

TEST(Ressampling, PinholeToSpherical_SphericalToPinhole)
{
  Image<RGBColor> spherical_image(512, 256);
  SphericalToPinholeIndirect(spherical_image);

  openMVG::cameras::Pinhole_Intrinsic pinhole_camera = spherical::ComputeCubicCameraIntrinsics(64);
  std::vector<Image<RGBColor>> cube_images;
  spherical::SphericalToCubic(spherical_image, pinhole_camera, cube_images, image::Sampler2d<image::SamplerNearest>());

  WriteImage("cube_0.png", cube_images[0]);
  WriteImage("cube_1.png", cube_images[1]);
  WriteImage("cube_2.png", cube_images[2]);
  WriteImage("cube_3.png", cube_images[3]);
  WriteImage("cube_4.png", cube_images[4]);
  WriteImage("cube_5.png", cube_images[5]);

  //WriteImage("test_pinhole_to_spherical_indirect.png", spherical_image);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
