// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"
#include "openMVG/image/image_io.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <array>
#include <string>

/// Compute a rectilinear camera focal for a given angular desired FoV and image size
double FocalFromPinholeHeight
(
  int h,
  double thetaMax = openMVG::D2R(60) // Camera FoV
)
{
  float f = 1.f;
  while ( thetaMax < atan2( h / (2 * f) , 1))
  {
    ++f;
  }
  return f;
}

const std::array<openMVG::Mat3, 6> GetCubicRotations()
{
  using namespace openMVG;
  return {
    RotationAroundY(D2R(0)),    // front
    RotationAroundY(D2R(90)),   // right
    RotationAroundY(D2R(180)),  // behind
    RotationAroundY(D2R(270)),  // left
    RotationAroundX(D2R(90)),   // up
    RotationAroundX(D2R(-90))   // down
  };
}

template <typename ImageT>
void SphericalToCubic
(
  const ImageT & equirectangular_image,
  std::vector<ImageT> & cube_images,
  openMVG::cameras::Pinhole_Intrinsic & pinhole_camera
)
{
  using namespace openMVG;
  using namespace openMVG::cameras;
  //
  // Initialize a camera model for each image domain
  // - the equirectangular panorama
  const Intrinsic_Spherical sphere_camera(
    equirectangular_image.Width() - 1, equirectangular_image.Height() - 1);
  // - the cube faces
  const double cubic_image_size = equirectangular_image.Width() / 4;
  const double focal = FocalFromPinholeHeight(cubic_image_size, D2R(45));
  const double principal_point_xy = cubic_image_size / 2;
  pinhole_camera = Pinhole_Intrinsic(cubic_image_size,
                                     cubic_image_size,
                                     focal,
                                     principal_point_xy,
                                     principal_point_xy);

  //
  // Perform backward/inverse rendering:
  // - For each cube face (rotation)
  // - Sample the panorama pixel by camera to camera bearing vector projection

  cube_images.resize(6, ImageT(cubic_image_size, cubic_image_size));

  // Get the rotation matrices corresponding to each cube face
  const std::array<Mat3, 6> rot_matrix = GetCubicRotations();

  #pragma omp parallel for
  for (int i_rot = 0; i_rot < rot_matrix.size(); ++i_rot)
  {
    auto & pinhole_image = cube_images[i_rot];
    // For every pinhole image pixels
    for (int x = 0; x < pinhole_image.Width(); ++x)
      for (int y = 0; y < pinhole_image.Height(); ++y)
      { // Project the pinhole bearing vector to the spherical camera
        const Vec3 pinhole_bearing = rot_matrix[i_rot] * pinhole_camera(Vec2(x, y));
        const Vec2 sphere_proj = sphere_camera.project(pinhole_bearing);
        // and use the corresponding pixel location in the panorama
        pinhole_image(y, x) = equirectangular_image(sphere_proj.y(), sphere_proj.x());
      }
  }
}

// Convert a spherical panorama to 6 perspective view (cubic panorama)
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string
    s_input_image,
    s_output_image;

  // required
  cmd.add( make_option('i', s_input_image, "input_image") );
  cmd.add( make_option('o', s_output_image, "output_image") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_image] the path to the spherical panorama\n"
    << "[-o|--output_image] the export directory path \n"
    << std::endl;
  }

  if (s_input_image.empty() || s_output_image.empty())
  {
    std::cerr << "input_image and output_image options must not be empty" << std::endl;
    return EXIT_FAILURE;
  }

  using namespace openMVG;

  image::Image<image::RGBColor> spherical_image;
  if (!ReadImage(s_input_image.c_str(), &spherical_image))
  {
    std::cerr << "Cannot read the input panoramic image: " << s_input_image << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<image::Image<image::RGBColor>> cube_images;
  openMVG::cameras::Pinhole_Intrinsic pinhole_camera;
  SphericalToCubic(spherical_image, cube_images, pinhole_camera);

  if (WriteImage("cube_0.png", cube_images[0]) &&
      WriteImage("cube_1.png", cube_images[1]) &&
      WriteImage("cube_2.png", cube_images[2]) &&
      WriteImage("cube_3.png", cube_images[3]) &&
      WriteImage("cube_4.png", cube_images[4]) &&
      WriteImage("cube_5.png", cube_images[5]))
  return EXIT_SUCCESS;

  return EXIT_FAILURE;
}
