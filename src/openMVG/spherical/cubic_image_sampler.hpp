// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/spherical/image_resampling.hpp"

#include <array>

namespace openMVG
{
namespace spherical
{

/// Compute a rectilinear camera focal from a given angular desired FoV
double FocalFromPinholeHeight
(
  int h,
  double fov_radian = openMVG::D2R(45) // Camera FoV
)
{
  return h / (2 * tan(fov_radian));
}

const static std::array<openMVG::Mat3,6> GetCubicRotations()
{
  using namespace openMVG;
  return {
    RotationAroundY(D2R(0)) ,   // front
    RotationAroundY(D2R(-90)),   // right
    RotationAroundY(D2R(-180)),  // behind
    RotationAroundY(D2R(-270)),  // left
    RotationAroundX(D2R(-90)),   // up
    RotationAroundX(D2R(+90))    // down
  };
}

openMVG::cameras::Pinhole_Intrinsic ComputeCubicCameraIntrinsics
(
  const int cubic_image_size,
  const double fov = D2R(45)
)
{
  const double focal = spherical::FocalFromPinholeHeight(cubic_image_size, fov);
  const double principal_point_xy = cubic_image_size / 2;
  return cameras::Pinhole_Intrinsic(cubic_image_size,
                                     cubic_image_size,
                                     focal,
                                     principal_point_xy,
                                     principal_point_xy);
}

template <typename ImageT, typename SamplerT>
void SphericalToCubic
(
  const ImageT & equirectangular_image,
  const openMVG::cameras::Pinhole_Intrinsic & pinhole_camera,
  std::vector<ImageT> & cube_images,
  const SamplerT sampler = image::Sampler2d<image::SamplerLinear>()
)
{
  const std::array<Mat3, 6> rot_matrix = GetCubicRotations();
  const std::vector<Mat3> rot_matrix_vec(rot_matrix.cbegin(), rot_matrix.cend());
  SphericalToPinholes(
      equirectangular_image,
      pinhole_camera,
      cube_images,
      rot_matrix_vec,
      sampler);
}

} // namespace spherical
} // namespace openMVG
