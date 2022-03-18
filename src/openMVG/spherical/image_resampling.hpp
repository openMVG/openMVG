// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"
#include "openMVG/image/sample.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{
namespace spherical
{

// Backward rendering of a pinhole image for a given rotation in a panorama
template <typename ImageT, typename SamplerT>
ImageT SphericalToPinhole
(
  const ImageT & equirectangular_image,
  const openMVG::cameras::Pinhole_Intrinsic & pinhole_camera,
  const Mat3 & rot_matrix = Mat3::Identity(),
  const SamplerT sampler = image::Sampler2d<image::SamplerLinear>()
)
{
  using namespace openMVG;
  using namespace openMVG::cameras;
  //
  // Initialize a camera model for each image domain
  // - the equirectangular panorama
  const Intrinsic_Spherical sphere_camera(equirectangular_image.Width()  - 1,
                                          equirectangular_image.Height() - 1);

  // Perform backward/inverse rendering:
  // - For each destination pixel in the pinhole image,
  //   compute where to pick the pixel in the panorama image.
  // This is done by using bearing vector computation

  const int renderer_image_size = pinhole_camera.h();
  ImageT pinhole_image(renderer_image_size, renderer_image_size);

  const int image_width = pinhole_image.Width();
  const int image_height = pinhole_image.Height();

  // Use image coordinate in a matrix to use OpenMVG camera bearing vector vectorization
  Mat2X xy_coords(2, static_cast<int>(image_width * image_height));
  for (int y = 0; y < image_height; ++y)
    for (int x = 0; x < image_width; ++x)
      xy_coords.col(x + image_width * y ) << x +.5 , y + .5;

  // Compute bearing vectors
  const Mat3X bearing_vectors = pinhole_camera(xy_coords);

  // Compute rotation bearings
  const Mat3X rotated_bearings = rot_matrix * bearing_vectors;
  // For every pinhole image pixels
  #pragma omp parallel for
  for (int it = 0; it < rotated_bearings.cols(); ++it)
  {
    // Project the bearing vector to the sphere
    const Vec2 sphere_proj = sphere_camera.project(rotated_bearings.col(it));
    // and use the corresponding pixel location in the panorama
    const Vec2 xy = xy_coords.col(it);
    if (equirectangular_image.Contains(sphere_proj.y(), sphere_proj.x()))
    {
      pinhole_image(xy.y(), xy.x()) = sampler(equirectangular_image, sphere_proj.y(), sphere_proj.x());
    }
  }
  return pinhole_image;
}

// Sample pinhole image from a panorama given some camera rotations
template <typename ImageT, typename SamplerT>
void SphericalToPinholes
(
  const ImageT & equirectangular_image,
  const openMVG::cameras::Pinhole_Intrinsic & pinhole_camera,
  std::vector<ImageT> & pinhole_images,
  const std::vector<Mat3> & rotations,
  const SamplerT sampler = image::Sampler2d<image::SamplerLinear>()
)
{
  pinhole_images.resize(6);
  // render each cube faces
  for (int i_rot = 0; i_rot < rotations.size(); ++i_rot)
  {
    pinhole_images[i_rot] = SphericalToPinhole(
        equirectangular_image,
        pinhole_camera,
        rotations[i_rot],
        sampler);
  }
}

} // namespace spherical
} // namespace openMVG
