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
#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace spherical {

/*
The spherical to Cartesian conversions in this file transform coordinates
according to this coordinate system:

-Y   +Z
 |   /
 |  /
 | /
 |/
 --------- +X

This 3D rectangular coordinate system is consistently used within computer
vision, and aligns with the 2D image coordinate system that places the origin at
the top-left of an image.

Some clarifying identities:
* (lon, lat) == (0, 0) along the +X axis
* (lon, lat) == (-pi, 0) along the -X axis
* (lon, lat) == (0, p/2) along the +Y axis
* (lon, lat) == (0, -p/2) along the -Y axis
* (lon, lat) == (pi/2, 0) along the +Z axis
* (lon, lat) == (3*pi/2, 0) along the -Z axis

Additionally, we assume a (H, W) equirectangular image relates to spherical
coordinates as:

           0 <--> -pi            W-1 <--> pi*(W-1/W)
              --------------------------
 0 <--> -pi/2 |                        |
              |                        |
              |                        |
H-1 <--> pi/2 |                        |
              --------------------------

*/

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

/*
Computes the forward gnomonic projection of a given point (lon, lat) on the
sphere to a plane tangent to (center_lon, center_lat). Output is in
coordinates on the plane. Equations courtesy of
https://mathworld.wolfram.com/GnomonicProjection.html
*/
inline const Vec2 ForwardGnomonicProjection(const Vec2 &lonlat,
                                            const Vec2 &center_lonlat) {
  const double cos_c = std::sin(center_lonlat[1]) * std::sin(lonlat[1]) +
                       std::cos(center_lonlat[1]) * std::cos(lonlat[1]) *
                           std::cos(lonlat[0] - center_lonlat[0]);
  const double x =
      std::cos(lonlat[1]) * std::sin(lonlat[0] - center_lonlat[0]) / cos_c;
  const double y = (std::cos(center_lonlat[1]) * std::sin(lonlat[1]) -
                    std::sin(center_lonlat[1]) * std::cos(lonlat[1]) *
                        std::cos(lonlat[0] - center_lonlat[0])) /
                   cos_c;

  return Vec2(x, y);
}

/*
  Computes the inverse gnomonic projection of a given point (x, y) on a plane
  tangent to (center_lon, center_lat). Output is in spherical coordinates
  (radians). Equations courtesy of
  https://mathworld.wolfram.com/GnomonicProjection.html
*/
inline const Vec2 InverseGnomonicProjection(const Vec2 &xy,
                                            const Vec2 &center_lonlat) {
  // Compute the inverse gnomonic projection of each (x,y) on a plane
  // centered at (lon, lat) [output is in spherical coordinates (radians)]
  const double rho = std::sqrt(xy[0] * xy[0] + xy[1] * xy[1]);
  const double nu = std::atan(rho);

  // Output longitude (module 2*PI)
  const double lon =
      NegMod(center_lonlat[0] +
                 std::atan2(
                     xy[0] * std::sin(nu),
                     rho * std::cos(center_lonlat[1]) * std::cos(nu) -
                         xy[1] * std::sin(center_lonlat[1]) * std::sin(nu)) +
                 M_PI,
             2 * M_PI) -
      M_PI;

  // Output latitude
  const double lat =
      std::asin(std::cos(nu) * std::sin(center_lonlat[1]) +
                xy[1] * std::sin(nu) * std::cos(center_lonlat[1]) / rho);
  return Vec2(lon, lat);
}

// Checks if a 2D point falls within a triangle
inline bool PointInTriangle2D(const Vec2 &pt, const Vec2 &v1, const Vec2 &v2,
                              const Vec2 &v3) {
  // Lambda to check which side of the triangle edge this point falls on
  auto sign = [](Vec2 pt, Vec2 v0, Vec2 v1) {
    return (pt[0] - v1[0]) * (v0[1] - v1[1]) -
           (v0[0] - v1[0]) * (pt[1] - v1[1]);
  };

  // Check each triangle edge
  const double d1 = sign(pt, v1, v2);
  const double d2 = sign(pt, v2, v3);
  const double d3 = sign(pt, v3, v1);

  // Returns true if all signs are consistent
  return !(((d1 < 0) || (d2 < 0) || (d3 < 0)) &&
           ((d1 > 0) || (d2 > 0) || (d3 > 0)));
}

// Checks if a origin-centered ray in 3D intersects a triangle in 3D using the
// spherical triangle test. Assumes face normal points outward
inline bool RayIntersectSphericalTriangle3D(const Vec3 &ray, const Vec3 &v1,
                                            const Vec3 &v2, const Vec3 &v3) {
  // Normals for each edge plane where the origin is a point on every plane
  const auto n1 = v1.cross(v2);
  const auto n2 = v2.cross(v3);
  const auto n3 = v3.cross(v1);

  // Check which side of the edge plan our ray falls on, denoted by sign
  const double d1 = ray.dot(n1);
  const double d2 = ray.dot(n2);
  const double d3 = ray.dot(n3);

  // Also check that the ray is in the same direction as the face normal (avoids
  // negative-ray ambiguity)
  const auto d_dir = ray.dot((v3 - v2).cross(v1 - v2));

  // Inside if all signs are consistent
  return (((d1 > 0) && (d2 > 0) && (d3 > 0)) ||
          ((d1 < 0) && (d2 < 0) && (d3 < 0))) &&
         (d_dir > 0);
}

// Tangent images class
class TangentImages {
 private:
  int base_level;   /* base icosahedron level */
  int sphere_level; /* icosahedron level corresponding to equirectangular input
                     */
  int rect_h;       /* height of the input equirectangular image */
  int rect_w;       /* width of the input equirectangular image */
  int dim;          /* dimension of the tangent images */
  int num;          /* number of tangent images */

  /***************************************************************************
  These constants are provided here to avoid two alternative solutions:

    (1) Re-implementing Loop subdivision from scratch for the icosahedron
    (2) Adding a new dependency on CGAL or another library that provides a
        Loop subdivision implementation

  Based on the results in Eder et al., "Tangent Images for Mitigating
  Spherical Distortion," CVPR 2020, we only really need the coordinates on the
  icosahedron for levels 0, 1, and 2, which makes the alternatives fairly
  cumbersome. Instead, this implementation hard-codes the necessary data. As
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

  /* These are the spherical coordinates of the vertices of each face of the
   * icosahedron. This is a "triangle soup" storage. */
  static const std::vector<double> kIcosahedronVerticesL0;
  static const std::vector<double> kIcosahedronVerticesL1;
  static const std::vector<double> kIcosahedronVerticesL2;

  /*
    Computes the inverse gnomonic projection of a window with dimensions (<kh>,
    <kw>) centered at <lonlat_in>, which is a flattened N x 2 matrix with the
    centers of N tangent images in spherical coordinates. The output is a
    flattened N x kh x kw x 2 tensor in row-major order.
  */
  void InverseGnomonicKernel(const std::vector<double> &lonlat_in, const int kh,
                             const int kw, const double res_lon,
                             const double res_lat,
                             std::vector<double> &lonlat_out) const;

  /*
    Creates a sampling map that describes, for each pixel on each tangent image,
    where to sample from the input equirectangular image.
  */
  void CreateEquirectangularToTangentImagesSampleMap(
      std::vector<double> &sampling_map) const;

  /*
    A selector function that grabs the correct hard-coded constant for the
    stored base level
  */
  double GetVertexAngularResolution() const;
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
  void GetIcosahedronVertices(std::vector<double> &triangles,
                              const int base_level) const;

  /*
    Returns the pixel coordinates (u, v) corresponding to the vertices of the
    icosahedral face associated with the <tangent_image_idx>-th tangent image
  */
  void ProjectFaceOntoTangentImage(const size_t tangent_image_idx, Vec2 &v0_uv,
                                   Vec2 &v1_uv, Vec2 &v2_uv) const;

  /*
    This function converts (u, v) coordinates on a tangent image given by
    <tangent_image_idx> to the corresponding spherical coordinate.
  */
  const Vec2 TangentUVToSpherical(const size_t tangent_image_idx,
                                  const Vec2 &uv) const;

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
  void CreateTangentImages(
      const ImageT &rect_img, std::vector<ImageT> &tangent_images,
      std::vector<image::Image<unsigned char>> *mask = nullptr) const;
  /*
    Create an equirectangular image by resampling from  the tangent images
  */
  template <typename ImageT>
  void ConvertTangentImagesToEquirectangular(
      const std::vector<ImageT> &tangent_images, ImageT &rect_img) const;

  /*
    This function converts (u, v) coordinates on a tangent image given by
    <tangent_image_idx> to the corresponding pixel in the equirectangular image.
  */
  const Vec2 TangentUVToEquirectangular(const size_t tangent_image_idx,
                                        const Vec2 &uv) const;

  /*
    Given the associated tangent image index (obtainable by the
    GetTangentImageIndex() function), this function converts (x, y) coordinates
    on the equirectangular image to tangent image coordinates (u, v).
  */
  const Vec2 EquirectangularToTangentUV(const size_t tangent_image_idx,
                                        const Vec2 &xy) const;
  /*
    Given a spherical coordinate, computes the index of the intersection face
    (and therefor the tangent image index). This functions leverages the
    subdivision structure of the icosahedron to do this look up in <= 20 + 4 * b
    checks.
  */
  size_t GetTangentImageIndex(const Vec2 &lonlat) const;

  /*
    Given a spherical coordinate, converts it to a pixel coordinate on the
    equirectangular image
  */
  const Vec2 ConvertSphericalToEquirectangular(const Vec2 &lonlat) const;

  /*
    Given a pixel coordinate on the equirectangular image, converts it to a
     spherical coordinate
  */
  const Vec2 ConvertEquirectangularToSpherical(const Vec2 &xy) const;

  /*
    This function takes a grayscale equirectangular image as input, computes
    sparse features on the tangent image representation defined by the calling
    object, and then converts the detected features' coordinates back to pixels
    on the equirectangular image. Optionally an equirectangular mask can be
    passed in as well where non-zero values denote where to compute features.
  */
  std::unique_ptr<features::Regions> ComputeFeaturesOnTangentImages(
      features::Image_describer &image_describer,
      const image::Image<unsigned char> &imageGray,
      image::Image<unsigned char> *feature_mask = nullptr) const;

  /*
    Returns the FOV of the tangent images in degrees
  */
  double FOV() const;

  /*
    Returns the base level
  */
  inline int BaseLevel() const { return this->base_level; }

  /*
    Returns the sphere level
  */
  inline int SphereLevel() const { return this->sphere_level; }

  /*
    Returns the number of tangent images at this base level
  */
  inline int Num() const { return this->num; }

  /*
    Returns the dimension of the tangent images in pixels
  */
  inline int Dim() const { return this->dim; }
};

/* TEMPLATE FUNCTION DEFINITIONS */

/*
  Function to create a set of tangent images from a given equirectangular image
*/
template <typename ImageT>
void TangentImages::CreateTangentImages(
    const ImageT &rect_img, std::vector<ImageT> &tangent_images,
    std::vector<image::Image<unsigned char>> *mask) const {
  // Allocate the output vector
  tangent_images.resize(this->num);
  if (mask) {
    mask->resize(this->num);
  }

  // Create the sampling maps for the tangent images
  std::vector<double> sampling_map;
  this->CreateEquirectangularToTangentImagesSampleMap(sampling_map);

  // Create bilinear sampler
  const image::Sampler2d<image::SamplerLinear> sampler;

  // For creating a face mask for each tangent image
  std::vector<double> face_vertices;
  this->GetIcosahedronVertices(face_vertices);
  std::vector<double> tangent_centers;
  this->GetTangentImageCenters(tangent_centers);

// Create each tangent image
#pragma omp parallel for
  for (size_t i = 0; i < this->num; i++) {
    // Initialize output image
    tangent_images[i] = ImageT(this->dim, this->dim);

    // Do some pre-computation for getting the mask for this tangent image if
    // desired
    Vec2 v0_uv;
    Vec2 v1_uv;
    Vec2 v2_uv;
    if (mask) {
      // Initialize mask image
      (*mask)[i] = image::Image<unsigned char>(this->dim, this->dim);

      // Compute gnomonic projection of face vertices
      this->ProjectFaceOntoTangentImage(i, v0_uv, v1_uv, v2_uv);
    }

    // Resample to each tangent image
    for (size_t j = 0; j < this->dim; j++) {
      for (size_t k = 0; k < this->dim; k++) {
        // Index in sample map
        const size_t map_idx =
            i * this->dim * this->dim * 2 + j * this->dim * 2 + k * 2;

        // Sample from the precomputed map
        tangent_images[i](j, k) =
            sampler(rect_img, sampling_map[map_idx + 1], sampling_map[map_idx]);

        // If we also want to generate masks of the unique regions on each
        // tangent image (i.e. the area that falls within the face of the
        // associated icosahedron)
        if (mask) {
          // Check if the sampling point falls within the associated triangular
          // face by doing a point-in-triangle test on the tangent image.

          // First get the sampling point in spherical coordinates
          const Vec2 lonlat = this->ConvertEquirectangularToSpherical(
              Vec2(sampling_map[map_idx], sampling_map[map_idx + 1]));

          // Apply the forward gnomonic projection onto this tangent image
          const Vec2 uv = ForwardGnomonicProjection(
              lonlat, Vec2(tangent_centers[i * 2], tangent_centers[i * 2 + 1]));

          // If the point falls within the projected face, then make the mask
          // non-zero
          if (PointInTriangle2D(uv, v0_uv, v1_uv, v2_uv)) {
            (*mask)[i](j, k) = 255;  // 255 chosen for visualizability
          }
        }
      }
    }
  }
}

/*
  Function to recreate an equirectangular image from a set of tangent images
*/
template <typename ImageT>
void TangentImages::ConvertTangentImagesToEquirectangular(
    const std::vector<ImageT> &tangent_images, ImageT &rect_img) const {
  // Initialize output image
  rect_img = ImageT(this->rect_w, this->rect_h);

  // Create bilinear sampler
  const image::Sampler2d<image::SamplerLinear> sampler;

// Iterate over all equirectangular image pixels
#pragma omp parallel for
  for (size_t i = 0; i < this->rect_h; i++) {
    for (size_t j = 0; j < this->rect_w; j++) {
      // Equirectangular pixel coordinates
      const auto xy = Vec2(j, i);

      // Convert the pixel coordinate to spherical coordinates
      const auto lonlat = this->ConvertEquirectangularToSpherical(xy);

      // Get the tangent image index to sample from
      const size_t tangent_image_idx = this->GetTangentImageIndex(lonlat);

      // Get the coordinates on the tangent image to sample from
      const auto uv = this->EquirectangularToTangentUV(tangent_image_idx, xy);

      // Sample from the specified tangent image
      rect_img(i, j) = sampler(tangent_images[tangent_image_idx], uv[1], uv[0]);
    }
  }
}
}  // namespace spherical
}  // namespace openMVG

#endif