// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Marc Eder.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SPHERICAL_SPHERICAL_HPP
#define OPENMVG_SPHERICAL_SPHERICAL_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"
#include "openMVG/spherical/ico_constants.hpp"

namespace openMVG {
namespace spherical {

using namespace openMVG::image;

// Floating point negative modulus operation
template <typename T>
const T fnegmod(const T lval, const T rval) {
  return fmod(fmod(lval, rval) + rval, rval);
}

// Rescale value to a new range
template <typename T>
const T Rescale(const T value, const T old_min, const T old_max,
                const T new_min, const T new_max) {
  return (new_max - new_min) * (value - old_min) / (old_max - old_min) +
         new_min;
}

// Compute tangent image dimension
const size_t TangentImageDimension(const size_t base_level,
                                   const size_t resolution_level) {
  return 1 << (resolution_level - base_level);
}

// Compute number of tangent images at a given base level
const size_t NumTangentImages(const size_t base_level) {
  return 20 * (1 << (2 * base_level));
}

// Sampling resolution of the tangent images
const double VertexAngularResolution(const size_t base_level) {
  switch (base_level) {
    case 0:
      return angResL0;
    case 1:
    default:
      return angResL1;
      // case 2:
      //   return angResL2;
  }
}

// Sampling resolution of the tangent images
const double SampleResolution(const size_t base_level,
                              const size_t resolution_level) {
  return VertexAngularResolution(base_level) /
         TangentImageDimension(base_level, resolution_level);
}

void TangentImageCenters(const size_t base_level,
                         std::vector<double> &centers) {
  switch (base_level) {
    case 0:
      centers = tangentCentersL0;
      break;
    case 1:
    default:
      centers = tangentCentersL1;
      break;
      // case 2:
      // centers = tangentCentersL0;
      // break;
  }
}

// lonlat_in: N x 2 matrix in row-major order
// lonlat_out: N x kh x kw x 2 tensor in row-major order
void GnomonicKernel(const std::vector<double> lonlat_in, const size_t kh,
                    const size_t kw, const double res_lon, const double res_lat,
                    std::vector<double> &lonlat_out) {
  // Allocate memory for output
  lonlat_out.resize(lonlat_in.size() * kh * kw);

#pragma omp parallel for
  // Go through each tangent image
  for (size_t i = 0; i < lonlat_in.size() / 2; i++) {
    // Tangent point for each tangent image in spherical coordinates
    const double lon = lonlat_in[2 * i];
    const double lat = lonlat_in[2 * i + 1];

    // j indexes rows in the tangent image
    for (size_t j = 0; j < kh; j++) {
      double y =
          (static_cast<double>(j) - static_cast<double>(kh) / 2.0) * res_lat;
      if (kh % 2 == 0) {
        y += res_lat / 2.0;
      }
      // k indexes columns in the tangent image
      for (size_t k = 0; k < kw; k++) {
        double x =
            (static_cast<double>(k) - static_cast<double>(kw) / 2.0) * res_lon;
        if (kw % 2 == 0) {
          x += res_lon / 2.0;
        }

        // Index in output vector
        const size_t vec_idx = i * kh * kw * 2 + j * kw * 2 + k * 2;

        // Compute the gnomonic projection of each (x,y) onto a plane centered
        // at (lon, lat) [output is in spherical coordinates (radians)]
        const double rho = std::sqrt(x * x + y * y);
        const double nu = std::atan(rho);

        // Compute output longitude (modulo 2*PI)
        lonlat_out[vec_idx] =
            fnegmod(lon +
                        std::atan2(x * std::sin(nu),
                                   rho * std::cos(lat) * std::cos(nu) -
                                       y * std::sin(lat) * std::sin(nu)) +
                        M_PI,
                    2 * M_PI) -
            M_PI;

        // Compute output latitude
        lonlat_out[vec_idx + 1] =
            std::asin(std::cos(nu) * std::sin(lat) +
                      y * std::sin(nu) * std::cos(lat) / rho);
      }
    }
  }
}

// Returns vector of (N x 4 [TL, TR, BL, BR] x 3) 3D coords
void TangentImageCorners(const size_t base_level, const size_t resolution_level,
                         std::vector<double> &corners) {
  // Get the centers of the tangent images in spherical coordinates
  std::vector<double> centers;
  TangentImageCenters(base_level, centers);

  // Get the tangent image info
  const double ang_res = VertexAngularResolution(base_level);
  const double tangent_dim =
      static_cast<double>(TangentImageDimension(base_level, resolution_level));
  const size_t num_tangent_imgs = NumTangentImages(base_level);

  // Compute the corners via the gnomonic projection
  const double d = ang_res * (tangent_dim - 1) / tangent_dim;
  GnomonicKernel(centers, 2, 2, d, d, corners);
}

void CreateEquirectangularToTangentImagesSampleMap(
    const size_t base_level, const size_t resolution_level, const size_t rect_h,
    const size_t rect_w, std::vector<double> &sampling_map) {
  // Tangent image dimension is 2^(s-b)
  const size_t tangent_dim =
      TangentImageDimension(base_level, resolution_level);

  // Allocate space for all the tangent image sample maps
  const size_t num_tangent_imgs = NumTangentImages(base_level);
  sampling_map.resize(num_tangent_imgs * tangent_dim * tangent_dim * 2);

  // Sampling resolution of the tangent images
  const double sample_resolution =
      SampleResolution(base_level, resolution_level);

  // Tangent image locations
  std::vector<double> tangent_centers;
  TangentImageCenters(base_level, tangent_centers);
  GnomonicKernel(tangent_centers, tangent_dim, tangent_dim, sample_resolution,
                 sample_resolution, sampling_map);

#pragma omp parallel for
  // Go through each tangent image
  for (size_t i = 0; i < num_tangent_imgs; i++) {
    // j indexes rows in the tangent image
    for (size_t j = 0; j < tangent_dim; j++) {
      // k indexes columns in the tangent image
      for (size_t k = 0; k < tangent_dim; k++) {
        // Index in sample map
        const size_t map_idx =
            i * tangent_dim * tangent_dim * 2 + j * tangent_dim * 2 + k * 2;

        // Convert the output to coordinates of an equirectangular image
        sampling_map[map_idx] = Rescale(sampling_map[map_idx], -M_PI, M_PI, 0.0,
                                        static_cast<double>(rect_w - 1));
        sampling_map[map_idx + 1] =
            Rescale(sampling_map[map_idx + 1], -M_PI / 2.0, M_PI / 2.0, 0.0,
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
void CreateTangentImages(const Image &src, const size_t resolution_level,
                         const size_t base_level, std::vector<Image> &out) {
  // TODO: should assert a range of allowed base levels (probably {0,1})

  // Allocate the output vector
  const size_t num_tangent_imgs =
      20 * std::pow(4L, static_cast<long>(base_level));
  out.resize(num_tangent_imgs);

  // Create the sampling maps for the tangent images
  std::vector<double> sampling_map;
  CreateEquirectangularToTangentImagesSampleMap(
      base_level, resolution_level, src.Height(), src.Width(), sampling_map);

  // Tangent image dimension is 2^(s-b)
  const size_t tangent_dim = 1 << (resolution_level - base_level);

  // Create bilinear samplerback
  const Sampler2d<SamplerLinear> sampler;

  // Create each tangent image
#pragma omp parallel for
  for (size_t i = 0; i < num_tangent_imgs; i++) {
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

void MapTangentUVToEquirectangular(const double u, const double v,
                                   const size_t base_level,
                                   const size_t resolution_level,
                                   const size_t rect_h, const size_t rect_w,
                                   double &x, double &y) {
  // Compute UV location as 3D coordinates

  // Convert 3D coordinate to spherical coordinate

  // Rescale spherical coordinate to equirectangular image
}
}  // namespace spherical
}  // namespace openMVG

#endif