// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Marc Eder.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SPHERICAL_TANGENT_IMAGES_HPP
#define OPENMVG_SPHERICAL_TANGENT_IMAGES_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"
#include "openMVG/features/regions.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace spherical {

// Floating point negative modulus operation
template <typename T>
inline const T NegMod(const T lval, const T rval) {
  return fmod(fmod(lval, rval) + rval, rval);
}

// Rescale value to a new range
template <typename T>
inline const T Rescale(const T value, const T old_min, const T old_max,
                       const T new_min, const T new_max) {
  return (new_max - new_min) * (value - old_min) / (old_max - old_min) +
         new_min;
}

// Converts spherical coordinates to 3D coordinates
inline const Vec3 ConvertSphericalTo3D(const Vec2 &lonlat) {
  const double x = std::cos(lonlat[1]) * std::cos(lonlat[0]);
  const double y = -std::sin(lonlat[1]);
  const double z = std::cos(lonlat[1]) * std::sin(lonlat[0]);
  return Vec3(x, y, z);
}

// Converts 3D coordinates to spherical coordinates
inline const Vec2 Convert3DToSpherical(const Vec3 &xyz) {
  const double lon = std::atan2(xyz[2], xyz[0]);
  const double lat =
      std::atan2(-xyz[1], std::sqrt(xyz[0] * xyz[0] + xyz[2] * xyz[2]));
  return Vec2(lon, lat);
}

class TangentImages {
 private:
  int base_level;   /* base icosahedron level */
  int sphere_level; /* icosahedron level corresponding to equirectangular input
                     */
  int rect_h;       /* height of the input equirectangular image */
  int rect_w;       /* width of the input equirectangular image */
  int dim;          /* dimension of the tangent images */
  int num;          /* number of tangent images */
  std::vector<Vec3> corners; /* 3D coords of the corners of the computed
                                tangent images */

  /***************************************************************************
  These constants are provided here to avoid two alternative solutions:

    (1) Re-implementing Loop subdivision from scratch for the icosahedron
    (2) Adding a new dependency on CGAL or another library that provides a
  Loop subdivision implementation

  Based on the results in Eder et al., "Tangent Images for Mitigating
  Spherical Distortion," CVPR 2020, we only really need the coordinates on the
  icosahedron for levels 0, 1, and 2, which makes the alternatives fairly
  cumbersome. Instead, this file explicitly provides the necessary data. As
  this is for internal operations of the class, these are private constants.
  ***************************************************************************/

  /* These are the average angle in radians between vertices in a L<b>
   * icosahedron. Only up to base 2 is supported. */
  static const double kVertexAngularResolutionL0;
  static const double kVertexAngularResolutionL1;
  static const double kVertexAngularResolutionL2;

  /* These are the spherical coordinates of the centers of each tangent image
   * for a L<b> icosahedron. Only up to base 2 is supported. */
  static const std::vector<double> kTangentImageCentersL0;
  static const std::vector<double> kTangentImageCentersL1;
  static const std::vector<double> kTangentImageCentersL2;

  static const std::vector<double> kIcosahedronVerticesL0;
  static const std::vector<double> kIcosahedronVerticesL1;
  static const std::vector<double> kIcosahedronVerticesL2;

  /*
    Compute the corners of the tangent images in 3D coordinates for subsequent
    use converting UV coordinates on a tangent plane to XY coordinates in the
    equirectangular image, storing the result in the class variable <corners>
   */
  void ComputeTangentImageCorners();

  /*
    Computes the gnomonic projection of a window with dimensions (<kh> x <kw>)
    centered at <lonlat_in>, which is a flattened N x 2 matrix with the centers
    of N tangent images in spherical coordinates. The output is a flattened N x
    kh x kw x 2 tensor in row-major order.
  */
  void InverseGnomonicKernel(const std::vector<double> lonlat_in, const int kh,
                             const int kw, const double res_lon,
                             const double res_lat,
                             std::vector<double> &lonlat_out) const;

  void CreateEquirectangularToTangentImagesSampleMap(
      std::vector<double> &sampling_map) const;

  /*
    A selector function that grabs the correct hard-coded constant for the
    stored base level
  */
  const double GetVertexAngularResolution() const;
  /*
    A selector function that grabs the correct hard-coded centers (in spherical
    coordinates) of the icosahedron faces at the stored base level. These are
    the centers of projection for creating the tangent images.
  */
  void GetTangentImageCenters(std::vector<double> &centers) const;

  /*
    A selector function to get the hard-coded vertices of the icosahedron at the
    stored base level. These are used for determining whether features are
    in-bounds when converting from tangent images to equirectangular
  */
  void GetIcosahedronVertices(std::vector<double> &triangles) const;

  /*
    This function converts (u, v) coordinates on a tangent image given by <tangent_image_idx> to the corresponding 3D point on the tangent plane to the sphere. This is used in conjunction with RayIntersectsTriangle to determine the valid region of a tangent image.
  */
  const Vec3 TangentUVTo3D(const size_t tangent_image_idx, const double u,
                           const double v) const;


  /*
    This function computes the 3D barycentric coordinates of point <pt> on triangle ABC in order to check if the point falls within the triangle. This is used when checking if a feature is within the valid region of the tangent image when converting back to equirectangular coordinates.
  */
  const bool RayIntersectsTriangle(const Vec3 &pt, const Vec3 &A, const Vec3 &B,
                                   const Vec3 &C) const;

  /*
    Checks if the ray to <pt> on the tangent image given by <tangent_image_idx>
    is within a valid region to be transferred back to the equirectangular
    image. This function is used to avoid converting features detected on the
    padded regions in the tangent image.
  */
  const bool IsValidRegion(const size_t tangent_image_idx,
                           const Vec3 &pt) const;

  /*
    Compute the number of tangent images at this base level
  */
  void ComputeNum();

  /*
    Compute the dimension of tangent images at this base level and spherical
    level
  */
  void ComputeDim();

  /*
    Get barycenters of faces in terms of 3D coordinates
  */
  void ComputeFaceCenters(std::vector<double> &centers) const;

 public:
  /* Constructor */
  TangentImages(const int base_level, const int sphere_level, const int rect_h,
                const int rect_w);

  /*
    Create the tangent images by resampling from the equirectangular image
  */
  template <typename ImageT>
  void CreateTangentImages(const ImageT &rect_img,
                           std::vector<ImageT> &tangent_images);

  /*
    Accepts a vector of Regions that describe the features detected on the set
    of tangent images. This function converts the coordinates of each feature
    from UV coordinates on a tangent image to XY coordinates on an
    equirectangular spherical image. These converted features are returned as a
    new Regions object
   */
  const std::unique_ptr<features::Regions>
  ConvertTangentImageFeaturesToEquirectangular(
      const std::vector<std::unique_ptr<features::Regions>> &regions) const;

  /*
    Returns the FOV of the tangent images in degrees
  */
  const double FOV() const;

  /*
    Returns the base level
  */
  inline const int BaseLevel() const { return this->base_level; }

  /*
    Returns the sphere level
  */
  inline const int SphereLevel() const { return this->sphere_level; }

  /*
    Returns the number of tangent images at this base level
  */
  inline const int Num() const { return this->num; }

  /*
    Returns the dimension of the tangent images in pixels
  */
  inline const int Dim() const { return this->dim; }
};

/* TEMPLATE FUNCTION DEFINITIONS */

/*
  Function to create a set of tangent images from a given equirectangular image
*/
template <typename ImageT>
void TangentImages::CreateTangentImages(const ImageT &rect_img,
                                        std::vector<ImageT> &tangent_images) {
  // Allocate the output vector
  const size_t num_tangent_imgs = this->num;
  tangent_images.resize(num_tangent_imgs);

  // Create the sampling maps for the tangent images
  std::vector<double> sampling_map;
  this->CreateEquirectangularToTangentImagesSampleMap(sampling_map);

  // Tangent image dimension is 2^(s-b)
  const size_t tangent_dim = this->dim;

  // Create bilinear samplerback
  const image::Sampler2d<image::SamplerLinear> sampler;

  // Create each tangent image
#pragma omp parallel for
  for (size_t i = 0; i < num_tangent_imgs; i++) {
    // Initialize output image
    tangent_images[i] = ImageT(tangent_dim, tangent_dim);

    // Resample to each tangent image
    for (size_t j = 0; j < tangent_dim; j++) {
      for (size_t k = 0; k < tangent_dim; k++) {
        // Index in sample map
        const size_t map_idx =
            i * tangent_dim * tangent_dim * 2 + j * tangent_dim * 2 + k * 2;

        // Sample from the precomputed map
        tangent_images[i](j, k) =
            sampler(rect_img, sampling_map[map_idx + 1], sampling_map[map_idx]);
      }
    }
  }
}

}  // namespace spherical
}  // namespace openMVG

#endif