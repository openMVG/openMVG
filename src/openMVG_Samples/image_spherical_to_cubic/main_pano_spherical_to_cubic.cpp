// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/spherical/cubic_image_sampler.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>

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

  const double cubic_image_size = spherical_image.Width() / 4; // Size can be customized
  const openMVG::cameras::Pinhole_Intrinsic pinhole_camera =
    spherical::ComputeCubicCameraIntrinsics(cubic_image_size);

  std::vector<image::Image<image::RGBColor>> cube_images(6);
  spherical::SphericalToCubic
  (
    spherical_image,
    pinhole_camera,
    cube_images,
    image::Sampler2d<image::SamplerNearest>()
  );

  if (WriteImage("cube_0.png", cube_images[0]) &&
      WriteImage("cube_1.png", cube_images[1]) &&
      WriteImage("cube_2.png", cube_images[2]) &&
      WriteImage("cube_3.png", cube_images[3]) &&
      WriteImage("cube_4.png", cube_images[4]) &&
      WriteImage("cube_5.png", cube_images[5]))
  return EXIT_SUCCESS;

  return EXIT_FAILURE;
}
