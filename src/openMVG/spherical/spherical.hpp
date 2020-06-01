// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014, 2015 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SPHERICAL_SPHERICAL_HPP
#define OPENMVG_SPHERICAL_SPHERICAL_HPP

#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"

namespace openMVG {
namespace spherical {

using namespace openMVG::image;

// Spherical coordinates of placement of tangent images for level 0 (20 x 2 row
// major order)
static const std::vector<double> tangentCentersL0 = {
    1.8850,  -0.9184, 0.6283,  -0.9184, -0.6283, -0.9184, -1.8850, -0.9184,
    3.1416,  -0.9184, 0.6283,  -0.1887, 1.8850,  -0.1887, 3.1416,  -0.1887,
    -1.8850, -0.1887, -0.6283, -0.1887, 1.2566,  0.9184,  2.5133,  0.9184,
    -2.5133, 0.9184,  -1.2566, 0.9184,  0.0000,  0.9184,  1.2566,  0.1887,
    2.5133,  0.1887,  -2.5133, 0.1887,  -1.2566, 0.1887,  0.0000,  0.1887};

// Spherical coordinates of placement of tangent images for level 1 (80 x 2 row
// major order)
static const std::vector<double> tangentCentersL1 = {
    1.8850,  -0.9184, 0.6283,  -0.9184, -0.6283, -0.9184, -1.8850, -0.9184,
    3.1416,  -0.9184, 0.6283,  -0.1887, 1.8850,  -0.1887, 3.1416,  -0.1887,
    -1.8850, -0.1887, -0.6283, -0.1887, 1.2566,  0.9184,  2.5133,  0.9184,
    -2.5133, 0.9184,  -1.2566, 0.9184,  0.0000,  0.9184,  1.2566,  0.1887,
    2.5133,  0.1887,  -2.5133, 0.1887,  -1.2566, 0.1887,  0.0000,  0.1887,
    2.2802,  -0.6998, 1.4897,  -0.6998, 1.8850,  -1.2655, 1.0236,  -0.6998,
    0.2331,  -0.6998, 0.6283,  -1.2655, -0.2331, -0.6998, -1.0236, -0.6998,
    -0.6283, -1.2655, -1.4897, -0.6998, -2.2802, -0.6998, -1.8850, -1.2655,
    3.1416,  -1.2655, -2.7464, -0.6998, 2.7464,  -0.6998, 0.9473,  -0.3506,
    0.6283,  0.1583,  0.3093,  -0.3506, 2.2040,  -0.3506, 1.8850,  0.1583,
    1.5660,  -0.3506, -2.8226, -0.3506, 3.1416,  0.1583,  2.8226,  -0.3506,
    -1.5660, -0.3506, -1.8850, 0.1583,  -2.2040, -0.3506, -0.3093, -0.3506,
    -0.6283, 0.1583,  -0.9473, -0.3506, 0.8614,  0.6998,  1.6519,  0.6998,
    1.2566,  1.2655,  2.1180,  0.6998,  2.9085,  0.6998,  2.5133,  1.2655,
    -2.9085, 0.6998,  -2.1180, 0.6998,  -2.5133, 1.2655,  -1.6519, 0.6998,
    -0.8614, 0.6998,  -1.2566, 1.2655,  0.0000,  1.2655,  -0.3952, 0.6998,
    0.3952,  0.6998,  1.2566,  -0.1583, 1.5756,  0.3506,  0.9376,  0.3506,
    2.5133,  -0.1583, 2.8323,  0.3506,  2.1943,  0.3506,  -2.5133, -0.1583,
    -2.1943, 0.3506,  -2.8323, 0.3506,  -1.2566, -0.1583, -0.9376, 0.3506,
    -1.5756, 0.3506,  0.0000,  -0.1583, 0.3190,  0.3506,  -0.3190, 0.3506};

// Angular resolution of icosahedron levels
static const double angResL0 = 1.1071487665176392; // radians
static const double angResL1 = 0.5891666412353516; // radians

// Floating point negative modulus operation
template <typename T>
const T fnegmod(const T lval, const T rval) {
  return fmod(fmod(lval, rval) + rval, rval);
}

// Rescale value to a new range
template <typename T>
T renormalize(const T value, const T old_min, const T old_max, const T new_min, const T new_max){
  return (new_max - new_min) * (value - old_min) / (old_max - old_min) +
         new_min;
}

void CreateSamplingMap(const size_t base_level, const size_t resolution_level,
                       const size_t rect_h, const size_t rect_w,
                       std::vector<double>& sampling_map) {
  // Tangent image dimension is 2^(s-b)
  const size_t tangent_dim = 1 << (resolution_level - base_level);

  // Allocate space for all the tangent image sample maps
  const size_t num_tangent_imgs =
      20 * std::pow(4L, static_cast<long>(base_level));
  sampling_map.resize(num_tangent_imgs * tangent_dim * tangent_dim * 2);

  // Resolution of the icosahedron
  double sample_resolution;
  switch(base_level)
  {
    case 0:
        sample_resolution = angResL0 / tangent_dim;
        break;
    case 1:
    default:
      sample_resolution = angResL1 / tangent_dim;
      break;
  }

  // Tangent image locations
  std::vector<double> tangent_centers;
  switch(base_level)
  {
    case 0:
        tangent_centers = tangentCentersL0;
        break;
    case 1:
    default:
        tangent_centers = tangentCentersL1;
        break;
  }

  // Go through each tangent image
  for (size_t i = 0; i < num_tangent_imgs; i++){
    // Tangent point for each tangent image in spherical coordinates
    const double lon = tangent_centers[2  * i];
    const double lat = tangent_centers[2 * i + 1];

    // j indexes rows in the tangent image
    for (size_t j = 0; j < tangent_dim; j++){
      const double y =
          (static_cast<double>(j) - static_cast<double>(tangent_dim) / 2.0) *
              sample_resolution +
          sample_resolution / 2.0;
      // k indexes columns in the tangent image
      for (size_t k = 0; k < tangent_dim; k++){
        const double x =
            (static_cast<double>(k) - static_cast<double>(tangent_dim) / 2.0) *
                sample_resolution +
            sample_resolution / 2.0;

        // Compute the gnomonic projection of each (x,y) onto a plane centered at (lon, lat) [output is in spherical coordinates (radians)]
        const double rho = std::sqrt(x * x + y * y);
        const double nu = std::atan(rho);

        // Index in sample map
        const size_t map_idx =
            i * tangent_dim * tangent_dim * 2 + j * tangent_dim * 2 + k * 2;

        // Compute output longitude (modulo 2*PI)
        sampling_map[map_idx] =
            fnegmod(lon +
                     std::atan2(x * std::sin(nu),
                                rho * std::cos(lat) * std::cos(nu) -
                                    y * std::sin(lat) * std::sin(nu)) +
                     M_PI,
                 2 * M_PI) -
            M_PI;

        // Compute output latitude
        sampling_map[map_idx + 1] =
            std::asin(std::cos(nu) * std::sin(lat) +
                      y * std::sin(nu) * std::cos(lat) / rho);

        // Convert the output to coordinates of an equirectangular image
        sampling_map[map_idx] =
            renormalize(sampling_map[map_idx], -M_PI, M_PI, 0.0, static_cast<double>(rect_w));
        sampling_map[map_idx + 1] =
            renormalize(sampling_map[map_idx + 1], -M_PI / 2.0, M_PI / 2.0, 0.0,
                        static_cast<double>(rect_h - 1));
      }
    }
  }
}

/**
 ** @brief Convert equirectangular image to tangent images
 ** @param src input equirectangular image
 ** @param input spherical level (in terms of icosahedron subdivision)
 ** @param base level
 ** @param out output tangent image set
 **/
template <typename Image>
void Equirectangular2Tangent(const Image& src, const size_t resolution_level,
                             const size_t base_level,
                             std::vector<Image>& out) {
  // TODO: should assert a range of allowed base levels (probably {0,1})

  // Allocate the output vector
  const size_t num_tangent_imgs =
      20 * std::pow(4L, static_cast<long>(base_level));
  out.resize(num_tangent_imgs);

  // Create the sampling maps for the tangent images
  std::vector<double> sampling_map;
  CreateSamplingMap(base_level, resolution_level, src.Height(), src.Width(),
                    sampling_map);

  // Tangent image dimension is 2^(s-b)
  const size_t tangent_dim = 1 << (resolution_level - base_level);

  // Create bilinear sampler
  const Sampler2d<SamplerLinear> sampler;

  // Create each tangent image
  for (size_t i = 0; i < num_tangent_imgs; i++){
    
    // Initialize output image
    out[i] = Image(tangent_dim, tangent_dim);

    // Resample to each tangent image
    for (size_t j = 0; j < tangent_dim; j++) {
      for (size_t k = 0; k < tangent_dim; k++) {
        // Index in sample map
        const size_t map_idx =
            i * tangent_dim * tangent_dim * 2 + j * tangent_dim * 2 + k * 2;

        // Sample from the precomputed map
        out[i](j, k) =
            sampler(src, sampling_map[map_idx + 1], sampling_map[map_idx]);
      }
    }
  }
}
}
}  // namespace openMVG

#endif