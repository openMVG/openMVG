// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Marc Eder.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/spherical/tangent_images.hpp"
#include <cmath>

namespace openMVG {
namespace spherical {

/* Constructor */
TangentImages::TangentImages(const int base_level, const int sphere_level,
                             const int rect_h, const int rect_w)
    : base_level(base_level),
      sphere_level(sphere_level),
      rect_h(rect_h),
      rect_w(rect_w) {
  // Ensure that the base level is not outside the supported range
  if (this->base_level < 0 || this->base_level > 2) {
    throw std::invalid_argument("Base level must be in {0, 1, 2}");
  }

  // Compute the number and dim
  this->ComputeNum();
  this->ComputeDim();
}

size_t TangentImages::GetTangentImageIndex(const Vec2 &lonlat) const {
  // Uses the nested structure of the icosahedron to efficiently search for
  // intersecting faces

  // Get the vertices of the level 0 icosahedron as a triangle soup
  std::vector<double> face_vertices;
  this->GetIcosahedronVertices(face_vertices, 0);

  // Find the first face intersection exhaustively
  size_t face_idx = -1;

  for (size_t i = 0; i < 20; i++) {
    // Get the vertices of this face in 3D coordinates
    const Vec3 v0 = ConvertSphericalTo3D(
        Vec2(face_vertices[6 * i], face_vertices[6 * i + 1]));
    const Vec3 v1 = ConvertSphericalTo3D(
        Vec2(face_vertices[6 * i + 2], face_vertices[6 * i + 3]));
    const Vec3 v2 = ConvertSphericalTo3D(
        Vec2(face_vertices[6 * i + 4], face_vertices[6 * i + 5]));
    // std::cout << "ray" << std::endl;
    // std::cout << lonlat.transpose() << std::endl;
    // std::cout << ConvertSphericalTo3D(lonlat).transpose() << std::endl;
    // std::cout << "corners" << std::endl;
    // std::cout << v0.transpose() << std::endl;
    // std::cout << v1.transpose() << std::endl;
    // std::cout << v2.transpose() << std::endl;
    if (RayIntersectSphericalTriangle3D(ConvertSphericalTo3D(lonlat), v0, v1,
                                        v2)) {
      face_idx = i;
      // If we are at a base level zero representation, just return, as our job
      // is done
      if (this->base_level == 0) {
        return face_idx;
      }
    }
  }

  // This should never happen as the icosahedron is water tight, so a spherical
  // ray will intersection with a face. Nevertheless....
  if (face_idx == -1) {
    throw std::runtime_error(
        "Something went wrong with face intersection search!");
  }

  // If this tangent image representation uses a base level above 0, we only
  // need to search 4 faces for each subsequent base level
  for (size_t i = 1; i <= this->base_level; i++) {
    // Get the vertices at this base level
    this->GetIcosahedronVertices(face_vertices, i);

    // At each subsequent base level, we can compute the range of possible
    // indices for the face intersections using this relation below. Because of
    // the way the hard-coded icosahedron is subdivided, each subsequent level
    // turns <face_idx> into 4 new faces numbered {<face_idx>, 20*4^(l-1) +
    // 3*<face_idx> + k}, where k={0,1,2}. Because it's guaranteed to be one of
    // these four faces, we just check the newly-numbered ones. If it's not one
    // of them, then it still <face_idx>.
    size_t start_face_idx = 20 * std::pow(4, i - 1) + 3 * face_idx;

    // Check each of the 3 subdivided faces for intersections. If it's none of
    // those 3, then it's the central face (which keeps the index from the
    // previous subdivision level)
    for (size_t j = 0; j < 3; j++) {
      // Current face index
      const size_t cur_face_idx = start_face_idx + j;

      // Get the vertices of this face in 3D coordinates
      const Vec3 v0 =
          ConvertSphericalTo3D(Vec2(face_vertices[6 * cur_face_idx],
                                    face_vertices[6 * cur_face_idx + 1]));
      const Vec3 v1 =
          ConvertSphericalTo3D(Vec2(face_vertices[6 * cur_face_idx + 2],
                                    face_vertices[6 * cur_face_idx + 3]));
      const Vec3 v2 =
          ConvertSphericalTo3D(Vec2(face_vertices[6 * cur_face_idx + 4],
                                    face_vertices[6 * cur_face_idx + 5]));

      // Check intersection at this face
      if (RayIntersectSphericalTriangle3D(ConvertSphericalTo3D(lonlat), v0, v1,
                                          v2)) {
        face_idx = cur_face_idx;

        // Early stopping option
        if (this->base_level == i) {
          return face_idx;
        }
      }
    }
    // If we've checked the 3 other faces, then the intersecting face is still
    // labeled as <face_idx>.
  }

  return face_idx;
}

/*
  Given the associated tangent image index (obtainable by the
  GetTangentImageIndex() function), this function converts (x, y) coordinates on
  the equirectangular image to tangent image coordinates (u, v).
*/
const Vec2 TangentImages::EquirectangularToTangentUV(
    const size_t tangent_image_idx, const Vec2 &xy) const {
  // Get the centers of the tangent images
  std::vector<double> tangent_centers;
  this->GetTangentImageCenters(tangent_centers);

  // Sampling resolution of the tangent images
  const double sample_resolution =
      this->GetVertexAngularResolution() / this->dim;

  // Get equirectangular coordinates in lon/lat
  const auto lonlat = this->ConvertEquirectangularToSpherical(xy);

  // Center of the tangent image (i.e. center of spherical projection)
  const double center_lon = tangent_centers[2 * tangent_image_idx];
  const double center_lat = tangent_centers[2 * tangent_image_idx + 1];

  // Project the coordinate onto the plane
  const auto uv_n =
      ForwardGnomonicProjection(lonlat, Vec2(center_lon, center_lat));

  // Convert the normalized projected coordinates to indices of the tangent
  // image pixel grid
  const double u = (uv_n[0] - sample_resolution / 2.0) / sample_resolution +
                   static_cast<double>(this->dim) / 2.0;
  const double v = (uv_n[1] - sample_resolution / 2.0) / sample_resolution +
                   static_cast<double>(this->dim) / 2.0;

  return Vec2(u, v);
}

void TangentImages::InverseGnomonicKernel(
    const std::vector<double> &lonlat_in, const int kh, const int kw,
    const double res_lon, const double res_lat,
    std::vector<double> &lonlat_out) const {
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
      const double y =
          (static_cast<double>(j) - static_cast<double>(kh) / 2.0) * res_lat +
          res_lat / 2.0;
      // k indexes columns in the tangent image
      for (size_t k = 0; k < kw; k++) {
        const double x =
            (static_cast<double>(k) - static_cast<double>(kw) / 2.0) * res_lon +
            res_lon / 2.0;

        // Index in output vector
        const size_t vec_idx = i * kh * kw * 2 + j * kw * 2 + k * 2;

        // Compute the inverse gnomonic projection of each (x,y) on a plane
        // centered at (lon, lat)
        const Vec2 out_sph_coords =
            InverseGnomonicProjection(Vec2(x, y), Vec2(lon, lat));

        // std::cout << "index: " << i << ", (" << j << ", " << k << ")"
        //           << std::endl;
        // std::cout << "out_sph_coords: " << out_sph_coords << std::endl;

        // Compute output longitude (modulo 2*PI)
        lonlat_out[vec_idx] = out_sph_coords[0];

        // Compute output latitude
        lonlat_out[vec_idx + 1] = out_sph_coords[1];
      }
    }
  }
}

/*
  Creates a resampling map to create tangent images
*/
void TangentImages::CreateEquirectangularToTangentImagesSampleMap(
    std::vector<double> &sampling_map) const {
  // Allocate space for all the tangent image sample maps
  sampling_map.resize(static_cast<size_t>(this->num) *
                      static_cast<size_t>(this->dim) *
                      static_cast<size_t>(this->dim) * 2);

  // Sampling resolution of the tangent images
  const double sample_resolution =
      this->GetVertexAngularResolution() / this->dim;

  // Get the centers of the tangent images
  std::vector<double> tangent_centers;
  this->GetTangentImageCenters(tangent_centers);

  // Compute the sample map according to the inverse gnomonic projection
  this->InverseGnomonicKernel(tangent_centers, this->dim, this->dim,
                              sample_resolution, sample_resolution,
                              sampling_map);

  // For each tangent image, convert the spherical sampling coordinates to the
  // coordinates of an equirectangular image
#pragma omp parallel for
  // Go through each tangent image
  for (size_t i = 0; i < this->num; i++) {
    // j indexes rows in the tangent image
    for (size_t j = 0; j < this->dim; j++) {
      // k indexes columns in the tangent image
      for (size_t k = 0; k < this->dim; k++) {
        // Index in sample map
        const size_t map_idx =
            i * this->dim * this->dim * 2 + j * this->dim * 2 + k * 2;

        // Sampling location in spherical coordinates
        double lon = sampling_map[map_idx];
        double lat = sampling_map[map_idx + 1];

        // Convert the output to coordinates of an equirectangular image
        const Vec2 xy = this->ConvertSphericalToEquirectangular(Vec2(lon, lat));

        // Store the output pixel coordinates on the equirectangular image
        sampling_map[map_idx] = xy[0];
        sampling_map[map_idx + 1] = xy[1];
      }
    }
  }
}

const Vec2 TangentImages::ConvertSphericalToEquirectangular(
    const Vec2 &lonlat) const {
  const double x = Rescale(lonlat[0], -M_PI, M_PI, 0.0,
                           static_cast<double>(this->rect_w - 1));
  const double y = Rescale(lonlat[1], -M_PI / 2.0, M_PI / 2.0, 0.0,
                           static_cast<double>(this->rect_h - 1));
  return Vec2(x, y);
}

const Vec2 TangentImages::ConvertEquirectangularToSpherical(
    const Vec2 &xy) const {
  const double lon =
      Rescale(xy[0], 0.0, static_cast<double>(this->rect_w - 1), -M_PI, M_PI);
  const double lat = Rescale(xy[1], 0.0, static_cast<double>(this->rect_h - 1),
                             -M_PI / 2.0, M_PI / 2.0);
  return Vec2(lon, lat);
}

void TangentImages::ProjectFaceOntoTangentImage(const size_t tangent_image_idx,
                                                Vec2 &v0_uv, Vec2 &v1_uv,
                                                Vec2 &v2_uv) const {
  // Get the vertices of the icosahedron
  std::vector<double> face_vertices;
  this->GetIcosahedronVertices(face_vertices);

  // Get the centers of each tangent image
  std::vector<double> tangent_centers;
  this->GetTangentImageCenters(tangent_centers);

  // Compute gnomonic projection of face vertices
  const Vec2 v0 = Vec2(face_vertices[6 * tangent_image_idx],
                       face_vertices[6 * tangent_image_idx + 1]);
  const Vec2 v1 = Vec2(face_vertices[6 * tangent_image_idx + 2],
                       face_vertices[6 * tangent_image_idx + 3]);
  const Vec2 v2 = Vec2(face_vertices[6 * tangent_image_idx + 4],
                       face_vertices[6 * tangent_image_idx + 5]);
  v0_uv = ForwardGnomonicProjection(
      v0, Vec2(tangent_centers[tangent_image_idx * 2],
               tangent_centers[tangent_image_idx * 2 + 1]));
  v1_uv = ForwardGnomonicProjection(
      v1, Vec2(tangent_centers[tangent_image_idx * 2],
               tangent_centers[tangent_image_idx * 2 + 1]));
  v2_uv = ForwardGnomonicProjection(
      v2, Vec2(tangent_centers[tangent_image_idx * 2],
               tangent_centers[tangent_image_idx * 2 + 1]));
}

/*
  Number of tangent images is 20(4^b), where b is the base level
*/
void TangentImages::ComputeNum() {
  this->num = 20 * (1 << (2 * this->base_level));
}

/*
  Tangent image dimension is 2 ^ (s-b), where s is sphere level and b is
  base level. Note, tangent images are square.
*/
void TangentImages::ComputeDim() {
  this->dim = 1 << (this->sphere_level - this->base_level);
}

/*
  Computes the field of view of these tangent images as an
  average over all of them as there is slight variation due to the singular
  points of the icosahedron after subdivision
*/
double TangentImages::FOV() const {
  // Accumulator variable for mean computation
  double accum = 0.0;

  // Go through each set of N tangent images
  for (size_t i = 0; i < this->num; i++) {
    Vec3 tl =
        ConvertSphericalTo3D(this->TangentUVToSpherical(i, Vec2(0.0, 0.0)));
    Vec3 tr = ConvertSphericalTo3D(
        this->TangentUVToSpherical(i, Vec2(0.0, this->dim)));

    // Normalize the rays to vectors from the origin
    tl.normalize();
    tr.normalize();

    // Compute the angle between the vectors to get the FOV. Because tangent
    // images are square, we only do this in one direction
    accum += std::acos(tl.dot(tr));
  }

  // Computer average in radians
  accum /= this->num;

  // Return FOV in degrees
  return 180 * accum / M_PI;
}

const Vec2 TangentImages::TangentUVToSpherical(const size_t tangent_image_idx,
                                               const Vec2 &uv) const {
  // Get the reference to the tangent image centers
  std::vector<double> tangent_centers;
  this->GetTangentImageCenters(tangent_centers);

  // Angular resolution between pixels in a tangent image
  const double sample_resolution =
      this->GetVertexAngularResolution() / this->dim;

  // Convert UV coordinates to angular distance and shift them so the origin is
  // at the center pixel
  const double un =
      (uv[0] - static_cast<double>(this->dim) / 2.0) * sample_resolution +
      sample_resolution / 2.0;
  const double vn =
      (uv[1] - static_cast<double>(this->dim) / 2.0) * sample_resolution +
      sample_resolution / 2.0;

  // Compute the inverse gnomonic projection and convert the output to
  // spherical coordinates
  return InverseGnomonicProjection(
      Vec2(un, vn), Vec2(tangent_centers[tangent_image_idx * 2],
                         tangent_centers[tangent_image_idx * 2 + 1]));
}

const Vec2 TangentImages::TangentUVToEquirectangular(
    const size_t tangent_image_idx, const Vec2 &uv) const {
  // Simply convert to spherical and then to equirectangular
  return this->ConvertSphericalToEquirectangular(
      this->TangentUVToSpherical(tangent_image_idx, uv));
}

// Average angle between vertices in a <base_level> icosahedron (radians)
double TangentImages::GetVertexAngularResolution() const {
  switch (this->base_level) {
    case 0:
      return kVertexAngularResolutionL0;
    case 1:
    default:
      return kVertexAngularResolutionL1;
    case 2:
      return kVertexAngularResolutionL2;
  }
}

// Assigns the provided vector to the vector representing the N x 2 matrix of
// tangent image center (i.e. tangent points to the sphere) in spherical
// coordinates (lon, lat)
void TangentImages::GetTangentImageCenters(std::vector<double> &centers) const {
  switch (this->base_level) {
    case 0:
      centers = kTangentImageCentersL0;
      break;
    case 1:
    default:
      centers = kTangentImageCentersL1;
      break;
    case 2:
      centers = kTangentImageCentersL2;
      break;
  }
}

// Assigns the provided vector to the vector representing the N x 3 x 2
// (lon, lat) tensor of tangent image vertices as a "triangle soup",
// automatically using the base level of this instances's representation.
void TangentImages::GetIcosahedronVertices(
    std::vector<double> &triangles) const {
  this->GetIcosahedronVertices(triangles, this->base_level);
}

// Assigns the provided vector to the vector representing the N x 3 x 2
// (lon, lat) tensor of tangent image vertices as a "triangle soup"
void TangentImages::GetIcosahedronVertices(std::vector<double> &triangles,
                                           const int base_level) const {
  switch (base_level) {
    case 0:
      triangles = kIcosahedronVerticesL0;
      break;
    case 1:
    default:
      triangles = kIcosahedronVerticesL1;
      break;
    case 2:
      triangles = kIcosahedronVerticesL2;
      break;
  }
}

std::unique_ptr<features::Regions>
TangentImages::ComputeFeaturesOnTangentImages(
    features::Image_describer &image_describer,
    const image::Image<unsigned char> &imageGray,
    image::Image<unsigned char> *feature_mask) const {
  // Create the tangent images and the valid region mask
  std::vector<image::Image<unsigned char>> t_images;
  std::vector<image::Image<unsigned char>> t_mask;
  this->CreateTangentImages(imageGray, t_images, &t_mask);

  // If a mask is provided, convert that to tangent images and AND it with the
  // valid region mask
  if (feature_mask) {
    std::vector<image::Image<unsigned char>> t_feature_mask;
    std::vector<image::Image<unsigned char>> t_feature_mask_mask;
    this->CreateTangentImages(*feature_mask, t_feature_mask,
                              &t_feature_mask_mask);
#pragma omp parallel for
    for (size_t i = 0; i < t_images.size(); i++) {
      for (size_t j = 0; j < t_images[i].Height(); j++) {
        for (size_t k = 0; k < t_images[i].Width(); k++) {
          if (t_mask[i](j, k) > 0 && t_feature_mask_mask[i](j, k) > 0) {
            t_mask[i](j, k) = 255;
          } else {
            t_mask[i](j, k) = 0;
          }
        }
      }
    }
  }

  // At this point we have out tangent images to detect and describe on as
  // <t_images> and our feature mask is incorporated into our valid region mask
  // in <t_mask>. So next, we compute the descriptors on each tangent
  // image

  // Image_describer returns a unique_ptr, but we want to keep track of these
  // over all tangent images, so we create a vector of raw pointers
  std::vector<std::unique_ptr<features::Regions>> t_regions;
  for (size_t i = 0; i < t_images.size(); i++) {
    // Compute descriptors only in the valid regions
    auto regions = image_describer.Describe(t_images[i], &(t_mask[i]));
    t_regions.push_back(std::move(regions));
  }

  // Output regions container should be the same type as the tangent image ones
  std::unique_ptr<features::Regions> regions_out =
      std::unique_ptr<features::Regions>(t_regions[0]->EmptyClone());

  // Go through all regions and convert them to equirectangular
  for (size_t i = 0; i < t_regions.size(); i++) {
    // Go through each feature in this region, copy it to the output regions
    // container and update its coordinates
    for (size_t j = 0; j < t_regions[i]->RegionCount(); j++) {
      // Grab the location of the features
      const auto coords = t_regions[i]->GetRegionPosition(j);

      // Copy the feature over to the end of the output regions container
      t_regions[i]->CopyRegion(j, regions_out.get());

      // Convert back to equirectangular coordinate
      const auto rect_coords = this->TangentUVToEquirectangular(i, coords);

      // Update the coordinates in the output regions container depending on the
      // region type
      if (dynamic_cast<features::SIFT_Regions *>(regions_out.get())) {
        // SIFT_Regions
        auto sift_regions =
            dynamic_cast<features::SIFT_Regions *>(regions_out.get());
        sift_regions->Features().back().coords() = rect_coords.cast<float>();
      } else if (dynamic_cast<features::AKAZE_Float_Regions *>(
                     regions_out.get())) {
        // AKAZE_Float_Regions
        auto akaze_float_regions =
            dynamic_cast<features::AKAZE_Float_Regions *>(regions_out.get());
        akaze_float_regions->Features().back().coords() =
            rect_coords.cast<float>();
      } else if (dynamic_cast<features::AKAZE_Liop_Regions *>(
                     regions_out.get())) {
        // AKAZE_Liop_Regions
        auto akaze_liop_regions =
            dynamic_cast<features::AKAZE_Liop_Regions *>(regions_out.get());
        akaze_liop_regions->Features().back().coords() =
            rect_coords.cast<float>();
      } else if (dynamic_cast<features::AKAZE_Binary_Regions *>(
                     regions_out.get())) {
        // AKAZE_Binary_Regions
        auto akaze_binary_regions =
            dynamic_cast<features::AKAZE_Binary_Regions *>(regions_out.get());
        akaze_binary_regions->Features().back().coords() =
            rect_coords.cast<float>();
      } else {
        // Throw runtime error if none of the above
        throw std::runtime_error(
            "Invalid Region type! Tangent images support SIFT_Regions, "
            "AKAZE_Binary_Regions, AKAZE_Liop_Regions, and "
            "AKAZE_Float_Regions.");
      }
    }
  }

  // Return the output regions
  return regions_out;
}

/***************************/
/* Constant initialization */
/***************************/
const double TangentImages::kVertexAngularResolutionL0 = 2.2142975330352783;
const double TangentImages::kVertexAngularResolutionL1 = 1.1071487665176392;
const double TangentImages::kVertexAngularResolutionL2 = 0.5891666412353516;

// Spherical coordinates of level 0 tangent images (20 x 2 [lon, lat] row
// major order)
const std::vector<double> TangentImages::kTangentImageCentersL0 = {
    1.8850,  -0.9184, 0.6283,  -0.9184, -0.6283, -0.9184, -1.8850, -0.9184,
    3.1416,  -0.9184, 0.6283,  -0.1887, 1.8850,  -0.1887, 3.1416,  -0.1887,
    -1.8850, -0.1887, -0.6283, -0.1887, 1.2566,  0.9184,  2.5133,  0.9184,
    -2.5133, 0.9184,  -1.2566, 0.9184,  0.0000,  0.9184,  1.2566,  0.1887,
    2.5133,  0.1887,  -2.5133, 0.1887,  -1.2566, 0.1887,  0.0000,  0.1887};

// Spherical coordinates of level 0 tangent images (80 x 2 [lon, lat] row
// major order)
const std::vector<double> TangentImages::kTangentImageCentersL1 = {
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

// Spherical coordinates of level 0 tangent images (320 x 2 [lon, lat] row
// major order)
const std::vector<double> TangentImages::kTangentImageCentersL2 = {
    1.8850e+00,  -9.1844e-01, 6.2832e-01,  -9.1844e-01, -6.2832e-01,
    -9.1844e-01, -1.8850e+00, -9.1844e-01, -3.1416e+00, -9.1844e-01,
    6.2832e-01,  -1.8871e-01, 1.8850e+00,  -1.8871e-01, -3.1416e+00,
    -1.8871e-01, -1.8850e+00, -1.8871e-01, -6.2832e-01, -1.8871e-01,
    1.2566e+00,  9.1844e-01,  2.5133e+00,  9.1844e-01,  -2.5133e+00,
    9.1844e-01,  -1.2566e+00, 9.1844e-01,  -2.2098e-10, 9.1844e-01,
    1.2566e+00,  1.8871e-01,  2.5133e+00,  1.8871e-01,  -2.5133e+00,
    1.8871e-01,  -1.2566e+00, 1.8871e-01,  1.3657e-10,  1.8871e-01,
    2.2870e+00,  -6.9424e-01, 1.4829e+00,  -6.9424e-01, 1.8850e+00,
    -1.2731e+00, 1.0304e+00,  -6.9424e-01, 2.2628e-01,  -6.9424e-01,
    6.2832e-01,  -1.2731e+00, -2.2628e-01, -6.9424e-01, -1.0304e+00,
    -6.9424e-01, -6.2832e-01, -1.2731e+00, -1.4829e+00, -6.9424e-01,
    -2.2870e+00, -6.9424e-01, -1.8850e+00, -1.2731e+00, 3.1416e+00,
    -1.2731e+00, -2.7396e+00, -6.9424e-01, 2.7396e+00,  -6.9424e-01,
    9.5466e-01,  -3.5380e-01, 6.2832e-01,  1.6593e-01,  3.0197e-01,
    -3.5380e-01, 2.2113e+00,  -3.5380e-01, 1.8850e+00,  1.6593e-01,
    1.5586e+00,  -3.5380e-01, -2.8152e+00, -3.5380e-01, 3.1416e+00,
    1.6593e-01,  2.8152e+00,  -3.5380e-01, -1.5586e+00, -3.5380e-01,
    -1.8850e+00, 1.6593e-01,  -2.2113e+00, -3.5380e-01, -3.0197e-01,
    -3.5380e-01, -6.2832e-01, 1.6593e-01,  -9.5466e-01, -3.5380e-01,
    8.5460e-01,  6.9424e-01,  1.6587e+00,  6.9424e-01,  1.2566e+00,
    1.2731e+00,  2.1112e+00,  6.9424e-01,  2.9153e+00,  6.9424e-01,
    2.5133e+00,  1.2731e+00,  -2.9153e+00, 6.9424e-01,  -2.1112e+00,
    6.9424e-01,  -2.5133e+00, 1.2731e+00,  -1.6587e+00, 6.9424e-01,
    -8.5460e-01, 6.9424e-01,  -1.2566e+00, 1.2731e+00,  0.0000e+00,
    1.2731e+00,  -4.0204e-01, 6.9424e-01,  4.0204e-01,  6.9424e-01,
    1.2566e+00,  -1.6593e-01, 1.5830e+00,  3.5380e-01,  9.3029e-01,
    3.5380e-01,  2.5133e+00,  -1.6593e-01, 2.8396e+00,  3.5380e-01,
    2.1869e+00,  3.5380e-01,  -2.5133e+00, -1.6593e-01, -2.1869e+00,
    3.5380e-01,  -2.8396e+00, 3.5380e-01,  -1.2566e+00, -1.6593e-01,
    -9.3029e-01, 3.5380e-01,  -1.5830e+00, 3.5380e-01,  0.0000e+00,
    -1.6593e-01, 3.2634e-01,  3.5380e-01,  -3.2634e-01, 3.5380e-01,
    1.8850e+00,  -7.3639e-01, 1.5948e+00,  -9.9082e-01, 2.1751e+00,
    -9.9082e-01, 6.2832e-01,  -7.3639e-01, 3.3815e-01,  -9.9082e-01,
    9.1849e-01,  -9.9082e-01, -6.2832e-01, -7.3639e-01, -9.1849e-01,
    -9.9082e-01, -3.3815e-01, -9.9082e-01, -1.8850e+00, -7.3639e-01,
    -2.1751e+00, -9.9082e-01, -1.5948e+00, -9.9082e-01, -2.8514e+00,
    -9.9082e-01, -3.1416e+00, -7.3639e-01, 2.8514e+00,  -9.9082e-01,
    7.8649e-01,  -9.5722e-02, 4.7015e-01,  -9.5722e-02, 6.2832e-01,
    -3.7076e-01, 2.0431e+00,  -9.5722e-02, 1.7268e+00,  -9.5722e-02,
    1.8850e+00,  -3.7076e-01, -2.9834e+00, -9.5722e-02, 2.9834e+00,
    -9.5722e-02, 3.1416e+00,  -3.7076e-01, -1.7268e+00, -9.5722e-02,
    -2.0431e+00, -9.5722e-02, -1.8850e+00, -3.7076e-01, -4.7015e-01,
    -9.5722e-02, -7.8649e-01, -9.5722e-02, -6.2832e-01, -3.7076e-01,
    1.2566e+00,  7.3639e-01,  1.5468e+00,  9.9082e-01,  9.6647e-01,
    9.9082e-01,  2.5133e+00,  7.3639e-01,  2.8034e+00,  9.9082e-01,
    2.2231e+00,  9.9082e-01,  -2.5133e+00, 7.3639e-01,  -2.2231e+00,
    9.9082e-01,  -2.8034e+00, 9.9082e-01,  -1.2566e+00, 7.3639e-01,
    -9.6647e-01, 9.9082e-01,  -1.5468e+00, 9.9082e-01,  -2.9017e-01,
    9.9082e-01,  -1.8101e-10, 7.3639e-01,  2.9017e-01,  9.9082e-01,
    1.4148e+00,  9.5722e-02,  1.2566e+00,  3.7076e-01,  1.0985e+00,
    9.5722e-02,  2.6714e+00,  9.5722e-02,  2.5133e+00,  3.7076e-01,
    2.3551e+00,  9.5722e-02,  -2.3551e+00, 9.5722e-02,  -2.5133e+00,
    3.7076e-01,  -2.6714e+00, 9.5722e-02,  -1.0985e+00, 9.5722e-02,
    -1.2566e+00, 3.7076e-01,  -1.4148e+00, 9.5722e-02,  1.5817e-01,
    9.5722e-02,  1.4389e-10,  3.7076e-01,  -1.5817e-01, 9.5722e-02,
    2.4160e+00,  -5.7429e-01, 2.0740e+00,  -6.3680e-01, 2.3739e+00,
    -8.6006e-01, 1.3539e+00,  -5.7429e-01, 1.3961e+00,  -8.6006e-01,
    1.6959e+00,  -6.3680e-01, 1.8850e+00,  -1.4317e+00, 2.2843e+00,
    -1.1600e+00, 1.4856e+00,  -1.1600e+00, 1.1594e+00,  -5.7429e-01,
    8.1741e-01,  -6.3680e-01, 1.1172e+00,  -8.6006e-01, 9.7241e-02,
    -5.7429e-01, 1.3942e-01,  -8.6006e-01, 4.3923e-01,  -6.3680e-01,
    6.2832e-01,  -1.4317e+00, 1.0276e+00,  -1.1600e+00, 2.2900e-01,
    -1.1600e+00, -9.7241e-02, -5.7429e-01, -4.3923e-01, -6.3680e-01,
    -1.3942e-01, -8.6006e-01, -1.1594e+00, -5.7429e-01, -1.1172e+00,
    -8.6006e-01, -8.1741e-01, -6.3680e-01, -6.2832e-01, -1.4317e+00,
    -2.2900e-01, -1.1600e+00, -1.0276e+00, -1.1600e+00, -1.3539e+00,
    -5.7429e-01, -1.6959e+00, -6.3680e-01, -1.3961e+00, -8.6006e-01,
    -2.4160e+00, -5.7429e-01, -2.3739e+00, -8.6006e-01, -2.0740e+00,
    -6.3680e-01, -1.8850e+00, -1.4317e+00, -1.4856e+00, -1.1600e+00,
    -2.2843e+00, -1.1600e+00, 3.1416e+00,  -1.4317e+00, -2.7423e+00,
    -1.1600e+00, 2.7423e+00,  -1.1600e+00, -2.6105e+00, -5.7429e-01,
    -2.9525e+00, -6.3680e-01, -2.6527e+00, -8.6006e-01, 2.6105e+00,
    -5.7429e-01, 2.6527e+00,  -8.6006e-01, 2.9525e+00,  -6.3680e-01,
    1.1119e+00,  -4.1650e-01, 9.4485e-01,  -1.7714e-01, 7.9745e-01,
    -4.5602e-01, 6.2832e-01,  3.2452e-01,  4.7190e-01,  8.0982e-02,
    7.8474e-01,  8.0982e-02,  1.4472e-01,  -4.1650e-01, 4.5918e-01,
    -4.5602e-01, 3.1179e-01,  -1.7714e-01, 2.3685e+00,  -4.1650e-01,
    2.2015e+00,  -1.7714e-01, 2.0541e+00,  -4.5602e-01, 1.8850e+00,
    3.2452e-01,  1.7285e+00,  8.0982e-02,  2.0414e+00,  8.0982e-02,
    1.4014e+00,  -4.1650e-01, 1.7158e+00,  -4.5602e-01, 1.5684e+00,
    -1.7714e-01, -2.6580e+00, -4.1650e-01, -2.8251e+00, -1.7714e-01,
    -2.9725e+00, -4.5602e-01, 3.1416e+00,  3.2452e-01,  2.9852e+00,
    8.0982e-02,  -2.9852e+00, 8.0982e-02,  2.6580e+00,  -4.1650e-01,
    2.9725e+00,  -4.5602e-01, 2.8251e+00,  -1.7714e-01, -1.4014e+00,
    -4.1650e-01, -1.5684e+00, -1.7714e-01, -1.7158e+00, -4.5602e-01,
    -1.8850e+00, 3.2452e-01,  -2.0414e+00, 8.0982e-02,  -1.7285e+00,
    8.0982e-02,  -2.3685e+00, -4.1650e-01, -2.0541e+00, -4.5602e-01,
    -2.2015e+00, -1.7714e-01, -1.4472e-01, -4.1650e-01, -3.1179e-01,
    -1.7714e-01, -4.5918e-01, -4.5602e-01, -6.2832e-01, 3.2452e-01,
    -7.8474e-01, 8.0982e-02,  -4.7190e-01, 8.0982e-02,  -1.1119e+00,
    -4.1650e-01, -7.9745e-01, -4.5602e-01, -9.4485e-01, -1.7714e-01,
    7.2556e-01,  5.7429e-01,  1.0675e+00,  6.3680e-01,  7.6774e-01,
    8.6006e-01,  1.7877e+00,  5.7429e-01,  1.7455e+00,  8.6006e-01,
    1.4457e+00,  6.3680e-01,  1.2566e+00,  1.4317e+00,  8.5732e-01,
    1.1600e+00,  1.6560e+00,  1.1600e+00,  1.9822e+00,  5.7429e-01,
    2.3242e+00,  6.3680e-01,  2.0244e+00,  8.6006e-01,  3.0444e+00,
    5.7429e-01,  3.0022e+00,  8.6006e-01,  2.7024e+00,  6.3680e-01,
    2.5133e+00,  1.4317e+00,  2.1140e+00,  1.1600e+00,  2.9126e+00,
    1.1600e+00,  -3.0444e+00, 5.7429e-01,  -2.7024e+00, 6.3680e-01,
    -3.0022e+00, 8.6006e-01,  -1.9822e+00, 5.7429e-01,  -2.0244e+00,
    8.6006e-01,  -2.3242e+00, 6.3680e-01,  -2.5133e+00, 1.4317e+00,
    -2.9126e+00, 1.1600e+00,  -2.1140e+00, 1.1600e+00,  -1.7877e+00,
    5.7429e-01,  -1.4457e+00, 6.3680e-01,  -1.7455e+00, 8.6006e-01,
    -7.2556e-01, 5.7429e-01,  -7.6774e-01, 8.6006e-01,  -1.0675e+00,
    6.3680e-01,  -1.2566e+00, 1.4317e+00,  -1.6560e+00, 1.1600e+00,
    -8.5732e-01, 1.1600e+00,  0.0000e+00,  1.4317e+00,  -3.9932e-01,
    1.1600e+00,  3.9932e-01,  1.1600e+00,  -5.3108e-01, 5.7429e-01,
    -1.8909e-01, 6.3680e-01,  -4.8890e-01, 8.6006e-01,  5.3108e-01,
    5.7429e-01,  4.8890e-01,  8.6006e-01,  1.8909e-01,  6.3680e-01,
    1.2566e+00,  -3.2452e-01, 1.4131e+00,  -8.0982e-02, 1.1002e+00,
    -8.0982e-02, 1.7402e+00,  4.1650e-01,  1.4258e+00,  4.5602e-01,
    1.5732e+00,  1.7714e-01,  7.7304e-01,  4.1650e-01,  9.4011e-01,
    1.7714e-01,  1.0875e+00,  4.5602e-01,  2.5133e+00,  -3.2452e-01,
    2.6697e+00,  -8.0982e-02, 2.3569e+00,  -8.0982e-02, 2.9969e+00,
    4.1650e-01,  2.6824e+00,  4.5602e-01,  2.8298e+00,  1.7714e-01,
    2.0297e+00,  4.1650e-01,  2.1967e+00,  1.7714e-01,  2.3441e+00,
    4.5602e-01,  -2.5133e+00, -3.2452e-01, -2.3569e+00, -8.0982e-02,
    -2.6697e+00, -8.0982e-02, -2.0297e+00, 4.1650e-01,  -2.3441e+00,
    4.5602e-01,  -2.1967e+00, 1.7714e-01,  -2.9969e+00, 4.1650e-01,
    -2.8298e+00, 1.7714e-01,  -2.6824e+00, 4.5602e-01,  -1.2566e+00,
    -3.2452e-01, -1.1002e+00, -8.0982e-02, -1.4131e+00, -8.0982e-02,
    -7.7304e-01, 4.1650e-01,  -1.0875e+00, 4.5602e-01,  -9.4011e-01,
    1.7714e-01,  -1.7402e+00, 4.1650e-01,  -1.5732e+00, 1.7714e-01,
    -1.4258e+00, 4.5602e-01,  0.0000e+00,  -3.2452e-01, 1.5642e-01,
    -8.0982e-02, -1.5642e-01, -8.0982e-02, 4.8359e-01,  4.1650e-01,
    1.6913e-01,  4.5602e-01,  3.1653e-01,  1.7714e-01,  -4.8359e-01,
    4.1650e-01,  -3.1653e-01, 1.7714e-01,  -1.6913e-01, 4.5602e-01};

// 20 x 3 x 2 (lon, lat) row major layout of level 0 icosahedron vertices
const std::vector<double> TangentImages::kIcosahedronVerticesL0 = {
    0.0000000,  -1.5707964, 2.5132742,  -0.4636476, 1.2566371,  -0.4636476,

    0.0000000,  -1.5707964, 1.2566371,  -0.4636476, 0.0000000,  -0.4636476,

    0.0000000,  -1.5707964, 0.0000000,  -0.4636476, -1.2566371, -0.4636476,

    0.0000000,  -1.5707964, -1.2566371, -0.4636476, -2.5132742, -0.4636476,

    0.0000000,  -1.5707964, -2.5132742, -0.4636476, 2.5132742,  -0.4636476,

    0.0000000,  -0.4636476, 1.2566371,  -0.4636476, 0.6283185,  0.4636476,

    1.2566371,  -0.4636476, 2.5132742,  -0.4636476, 1.8849555,  0.4636476,

    2.5132742,  -0.4636476, -2.5132742, -0.4636476, 3.1415927,  0.4636476,

    -2.5132742, -0.4636476, -1.2566371, -0.4636476, -1.8849555, 0.4636476,

    -1.2566371, -0.4636476, 0.0000000,  -0.4636476, -0.6283185, 0.4636476,

    0.0000000,  1.5707964,  0.6283185,  0.4636476,  1.8849555,  0.4636476,

    0.0000000,  1.5707964,  1.8849555,  0.4636476,  3.1415927,  0.4636476,

    0.0000000,  1.5707964,  3.1415927,  0.4636476,  -1.8849555, 0.4636476,

    0.0000000,  1.5707964,  -1.8849555, 0.4636476,  -0.6283185, 0.4636476,

    0.0000000,  1.5707964,  -0.6283185, 0.4636476,  0.6283185,  0.4636476,

    1.2566371,  -0.4636476, 1.8849555,  0.4636476,  0.6283185,  0.4636476,

    3.1415927,  0.4636476,  1.8849555,  0.4636476,  2.5132742,  -0.4636476,

    -1.8849555, 0.4636476,  3.1415927,  0.4636476,  -2.5132742, -0.4636476,

    -0.6283185, 0.4636476,  -1.8849555, 0.4636476,  -1.2566371, -0.4636476,

    0.6283185,  0.4636476,  -0.6283185, 0.4636476,  0.0000000,  -0.4636476};

// 80 x 3 x 2 (lon, lat) row major layout of level 1 icosahedron vertices
const std::vector<double> TangentImages::kIcosahedronVerticesL1 = {
    2.5132742,  -1.0172219, 1.8849555,  -0.5535744, 1.2566371,  -1.0172219,

    1.2566371,  -1.0172219, 0.6283185,  -0.5535744, 0.0000000,  -1.0172219,

    0.0000000,  -1.0172219, -0.6283185, -0.5535744, -1.2566371, -1.0172219,

    -1.2566371, -1.0172219, -1.8849555, -0.5535744, -2.5132742, -1.0172219,

    2.5132742,  -1.0172219, -2.5132742, -1.0172219, 3.1415927,  -0.5535743,

    0.6283185,  -0.5535744, 0.9424778,  -0.0000000, 0.3141592,  -0.0000000,

    1.8849555,  -0.5535744, 2.1991148,  -0.0000000, 1.5707964,  -0.0000000,

    3.1415927,  -0.5535743, -2.8274333, -0.0000000, 2.8274333,  -0.0000000,

    -1.8849555, -0.5535744, -1.5707964, -0.0000000, -2.1991148, -0.0000000,

    -0.6283185, -0.5535744, -0.3141592, -0.0000000, -0.9424778, -0.0000000,

    0.6283185,  1.0172219,  1.2566371,  0.5535744,  1.8849555,  1.0172219,

    1.8849555,  1.0172219,  2.5132742,  0.5535744,  3.1415927,  1.0172219,

    3.1415927,  1.0172219,  -2.5132742, 0.5535744,  -1.8849555, 1.0172219,

    -1.8849555, 1.0172219,  -1.2566371, 0.5535744,  -0.6283185, 1.0172219,

    0.6283185,  1.0172219,  -0.6283185, 1.0172219,  0.0000000,  0.5535743,

    0.9424778,  -0.0000000, 1.5707964,  -0.0000000, 1.2566371,  0.5535744,

    2.1991148,  -0.0000000, 2.8274333,  -0.0000000, 2.5132742,  0.5535744,

    -2.8274333, -0.0000000, -2.1991148, -0.0000000, -2.5132742, 0.5535744,

    -1.5707964, -0.0000000, -0.9424778, -0.0000000, -1.2566371, 0.5535744,

    -0.3141592, -0.0000000, 0.3141592,  -0.0000000, 0.0000000,  0.5535743,

    2.5132742,  -1.0172219, 2.5132742,  -0.4636476, 1.8849555,  -0.5535744,

    1.8849555,  -0.5535744, 1.2566371,  -0.4636476, 1.2566371,  -1.0172219,

    1.2566371,  -1.0172219, 3.1415927,  -1.5707964, 2.5132742,  -1.0172219,

    1.2566371,  -1.0172219, 1.2566371,  -0.4636476, 0.6283185,  -0.5535744,

    0.6283185,  -0.5535744, 0.0000000,  -0.4636476, 0.0000000,  -1.0172219,

    0.0000000,  -1.0172219, 3.1415927,  -1.5707964, 1.2566371,  -1.0172219,

    0.0000000,  -1.0172219, 0.0000000,  -0.4636476, -0.6283185, -0.5535744,

    -0.6283185, -0.5535744, -1.2566371, -0.4636476, -1.2566371, -1.0172219,

    -1.2566371, -1.0172219, 3.1415927,  -1.5707964, 0.0000000,  -1.0172219,

    -1.2566371, -1.0172219, -1.2566371, -0.4636476, -1.8849555, -0.5535744,

    -1.8849555, -0.5535744, -2.5132742, -0.4636476, -2.5132742, -1.0172219,

    -2.5132742, -1.0172219, 3.1415927,  -1.5707964, -1.2566371, -1.0172219,

    2.5132742,  -1.0172219, 3.1415927,  -1.5707964, -2.5132742, -1.0172219,

    -2.5132742, -1.0172219, -2.5132742, -0.4636476, 3.1415927,  -0.5535743,

    3.1415927,  -0.5535743, 2.5132742,  -0.4636476, 2.5132742,  -1.0172219,

    0.6283185,  -0.5535744, 1.2566371,  -0.4636476, 0.9424778,  -0.0000000,

    0.9424778,  -0.0000000, 0.6283185,  0.4636476,  0.3141592,  -0.0000000,

    0.3141592,  -0.0000000, 0.0000000,  -0.4636476, 0.6283185,  -0.5535744,

    1.8849555,  -0.5535744, 2.5132742,  -0.4636476, 2.1991148,  -0.0000000,

    2.1991148,  -0.0000000, 1.8849555,  0.4636476,  1.5707964,  -0.0000000,

    1.5707964,  -0.0000000, 1.2566371,  -0.4636476, 1.8849555,  -0.5535744,

    3.1415927,  -0.5535743, -2.5132742, -0.4636476, -2.8274333, -0.0000000,

    -2.8274333, -0.0000000, 3.1415927,  0.4636476,  2.8274333,  -0.0000000,

    2.8274333,  -0.0000000, 2.5132742,  -0.4636476, 3.1415927,  -0.5535743,

    -1.8849555, -0.5535744, -1.2566371, -0.4636476, -1.5707964, -0.0000000,

    -1.5707964, -0.0000000, -1.8849555, 0.4636476,  -2.1991148, -0.0000000,

    -2.1991148, -0.0000000, -2.5132742, -0.4636476, -1.8849555, -0.5535744,

    -0.6283185, -0.5535744, 0.0000000,  -0.4636476, -0.3141592, -0.0000000,

    -0.3141592, -0.0000000, -0.6283185, 0.4636476,  -0.9424778, -0.0000000,

    -0.9424778, -0.0000000, -1.2566371, -0.4636476, -0.6283185, -0.5535744,

    0.6283185,  1.0172219,  0.6283185,  0.4636476,  1.2566371,  0.5535744,

    1.2566371,  0.5535744,  1.8849555,  0.4636476,  1.8849555,  1.0172219,

    1.8849555,  1.0172219,  0.0000000,  1.5707964,  0.6283185,  1.0172219,

    1.8849555,  1.0172219,  1.8849555,  0.4636476,  2.5132742,  0.5535744,

    2.5132742,  0.5535744,  3.1415927,  0.4636476,  3.1415927,  1.0172219,

    3.1415927,  1.0172219,  0.0000000,  1.5707964,  1.8849555,  1.0172219,

    3.1415927,  1.0172219,  3.1415927,  0.4636476,  -2.5132742, 0.5535744,

    -2.5132742, 0.5535744,  -1.8849555, 0.4636476,  -1.8849555, 1.0172219,

    -1.8849555, 1.0172219,  0.0000000,  1.5707964,  3.1415927,  1.0172219,

    -1.8849555, 1.0172219,  -1.8849555, 0.4636476,  -1.2566371, 0.5535744,

    -1.2566371, 0.5535744,  -0.6283185, 0.4636476,  -0.6283185, 1.0172219,

    -0.6283185, 1.0172219,  0.0000000,  1.5707964,  -1.8849555, 1.0172219,

    0.6283185,  1.0172219,  0.0000000,  1.5707964,  -0.6283185, 1.0172219,

    -0.6283185, 1.0172219,  -0.6283185, 0.4636476,  0.0000000,  0.5535743,

    0.0000000,  0.5535743,  0.6283185,  0.4636476,  0.6283185,  1.0172219,

    0.9424778,  -0.0000000, 1.2566371,  -0.4636476, 1.5707964,  -0.0000000,

    1.5707964,  -0.0000000, 1.8849555,  0.4636476,  1.2566371,  0.5535744,

    1.2566371,  0.5535744,  0.6283185,  0.4636476,  0.9424778,  -0.0000000,

    2.1991148,  -0.0000000, 2.5132742,  -0.4636476, 2.8274333,  -0.0000000,

    2.8274333,  -0.0000000, 3.1415927,  0.4636476,  2.5132742,  0.5535744,

    2.5132742,  0.5535744,  1.8849555,  0.4636476,  2.1991148,  -0.0000000,

    -2.8274333, -0.0000000, -2.5132742, -0.4636476, -2.1991148, -0.0000000,

    -2.1991148, -0.0000000, -1.8849555, 0.4636476,  -2.5132742, 0.5535744,

    -2.5132742, 0.5535744,  3.1415927,  0.4636476,  -2.8274333, -0.0000000,

    -1.5707964, -0.0000000, -1.2566371, -0.4636476, -0.9424778, -0.0000000,

    -0.9424778, -0.0000000, -0.6283185, 0.4636476,  -1.2566371, 0.5535744,

    -1.2566371, 0.5535744,  -1.8849555, 0.4636476,  -1.5707964, -0.0000000,

    -0.3141592, -0.0000000, 0.0000000,  -0.4636476, 0.3141592,  -0.0000000,

    0.3141592,  -0.0000000, 0.6283185,  0.4636476,  0.0000000,  0.5535743,

    0.0000000,  0.5535743,  -0.6283185, 0.4636476,  -0.3141592, -0.0000000};

// 320 x 3 x 2 (lon, lat) row major layout of level 2 icosahedron vertices
const std::vector<double> TangentImages::kIcosahedronVerticesL2 = {
    2.1106849,  -0.8159013, 1.6592263,  -0.8159014, 1.8849555,  -1.0964646,

    0.8540478,  -0.8159013, 0.4025893,  -0.8159013, 0.6283185,  -1.0964646,

    -0.4025893, -0.8159013, -0.8540478, -0.8159013, -0.6283185, -1.0964646,

    -1.6592263, -0.8159014, -2.1106849, -0.8159013, -1.8849555, -1.0964646,

    3.1415927,  -1.0964645, -2.9158635, -0.8159013, 2.9158635,  -0.8159013,

    0.7883530,  -0.2750545, 0.6283185,  -0.0106841, 0.4682840,  -0.2750545,

    2.0449901,  -0.2750545, 1.8849555,  -0.0106841, 1.7249211,  -0.2750545,

    -2.9815581, -0.2750545, 3.1415927,  -0.0106841, 2.9815581,  -0.2750545,

    -1.7249211, -0.2750545, -1.8849555, -0.0106841, -2.0449901, -0.2750545,

    -0.4682840, -0.2750545, -0.6283185, -0.0106841, -0.7883530, -0.2750545,

    1.0309078,  0.8159013,  1.4823663,  0.8159014,  1.2566371,  1.0964646,

    2.2875450,  0.8159013,  2.7390034,  0.8159013,  2.5132742,  1.0964646,

    -2.7390034, 0.8159013,  -2.2875450, 0.8159013,  -2.5132742, 1.0964646,

    -1.4823663, 0.8159014,  -1.0309078, 0.8159013,  -1.2566371, 1.0964646,

    0.0000000,  1.0964645,  -0.2257293, 0.8159013,  0.2257293,  0.8159013,

    1.2566371,  0.0106841,  1.4166715,  0.2750545,  1.0966026,  0.2750545,

    2.5132742,  0.0106841,  2.6733086,  0.2750545,  2.3532395,  0.2750545,

    -2.5132742, 0.0106841,  -2.3532395, 0.2750545,  -2.6733086, 0.2750545,

    -1.2566371, 0.0106841,  -1.0966026, 0.2750545,  -1.4166715, 0.2750545,

    0.0000000,  0.0106841,  0.1600345,  0.2750545,  -0.1600345, 0.2750545,

    2.5132742,  -0.7203593, 2.2301655,  -0.5267599, 2.1106849,  -0.8159013,

    1.5397458,  -0.5267599, 1.2566371,  -0.7203593, 1.6592263,  -0.8159014,

    1.2566371,  -1.3140846, 2.5132742,  -1.3140846, 1.8849555,  -1.0964646,

    1.2566371,  -0.7203593, 0.9735283,  -0.5267599, 0.8540478,  -0.8159013,

    0.2831088,  -0.5267599, 0.0000000,  -0.7203593, 0.4025893,  -0.8159013,

    0.0000000,  -1.3140846, 1.2566371,  -1.3140846, 0.6283185,  -1.0964646,

    0.0000000,  -0.7203593, -0.2831088, -0.5267599, -0.4025893, -0.8159013,

    -0.9735283, -0.5267599, -1.2566371, -0.7203593, -0.8540478, -0.8159013,

    -1.2566371, -1.3140846, 0.0000000,  -1.3140846, -0.6283185, -1.0964646,

    -1.2566371, -0.7203593, -1.5397458, -0.5267599, -1.6592263, -0.8159014,

    -2.2301655, -0.5267599, -2.5132742, -0.7203593, -2.1106849, -0.8159013,

    -2.5132742, -1.3140846, -1.2566371, -1.3140846, -1.8849555, -1.0964646,

    2.5132742,  -1.3140846, -2.5132742, -1.3140846, 3.1415927,  -1.0964645,

    -2.5132742, -0.7203593, -2.7963829, -0.5267599, -2.9158635, -0.8159013,

    2.7963829,  -0.5267599, 2.5132742,  -0.7203593, 2.9158635,  -0.8159013,

    0.9735283,  -0.5267599, 1.1019347,  -0.2514759, 0.7883530,  -0.2750545,

    0.7830209,  0.2514759,  0.4736161,  0.2514759,  0.6283185,  -0.0106841,

    0.1547024,  -0.2514759, 0.2831088,  -0.5267599, 0.4682840,  -0.2750545,

    2.2301655,  -0.5267599, 2.3585718,  -0.2514759, 2.0449901,  -0.2750545,

    2.0396581,  0.2514759,  1.7302531,  0.2514759,  1.8849555,  -0.0106841,

    1.4113395,  -0.2514759, 1.5397458,  -0.5267599, 1.7249211,  -0.2750545,

    -2.7963829, -0.5267599, -2.6679766, -0.2514759, -2.9815581, -0.2750545,

    -2.9868903, 0.2514759,  2.9868903,  0.2514759,  3.1415927,  -0.0106841,

    2.6679766,  -0.2514759, 2.7963829,  -0.5267599, 2.9815581,  -0.2750545,

    -1.5397458, -0.5267599, -1.4113395, -0.2514759, -1.7249211, -0.2750545,

    -1.7302531, 0.2514759,  -2.0396581, 0.2514759,  -1.8849555, -0.0106841,

    -2.3585718, -0.2514759, -2.2301655, -0.5267599, -2.0449901, -0.2750545,

    -0.2831088, -0.5267599, -0.1547024, -0.2514759, -0.4682840, -0.2750545,

    -0.4736161, 0.2514759,  -0.7830209, 0.2514759,  -0.6283185, -0.0106841,

    -1.1019347, -0.2514759, -0.9735283, -0.5267599, -0.7883530, -0.2750545,

    0.6283185,  0.7203593,  0.9114273,  0.5267599,  1.0309078,  0.8159013,

    1.6018468,  0.5267599,  1.8849555,  0.7203593,  1.4823663,  0.8159014,

    1.8849555,  1.3140846,  0.6283185,  1.3140846,  1.2566371,  1.0964646,

    1.8849555,  0.7203593,  2.1680644,  0.5267599,  2.2875450,  0.8159013,

    2.8584838,  0.5267599,  3.1415927,  0.7203593,  2.7390034,  0.8159013,

    3.1415927,  1.3140846,  1.8849555,  1.3140846,  2.5132742,  1.0964646,

    3.1415927,  0.7203593,  -2.8584838, 0.5267599,  -2.7390034, 0.8159013,

    -2.1680644, 0.5267599,  -1.8849555, 0.7203593,  -2.2875450, 0.8159013,

    -1.8849555, 1.3140846,  3.1415927,  1.3140846,  -2.5132742, 1.0964646,

    -1.8849555, 0.7203593,  -1.6018468, 0.5267599,  -1.4823663, 0.8159014,

    -0.9114273, 0.5267599,  -0.6283185, 0.7203593,  -1.0309078, 0.8159013,

    -0.6283185, 1.3140846,  -1.8849555, 1.3140846,  -1.2566371, 1.0964646,

    0.6283185,  1.3140846,  -0.6283185, 1.3140846,  0.0000000,  1.0964645,

    -0.6283185, 0.7203593,  -0.3452097, 0.5267599,  -0.2257293, 0.8159013,

    0.3452097,  0.5267599,  0.6283185,  0.7203593,  0.2257293,  0.8159013,

    1.1019347,  -0.2514759, 1.4113395,  -0.2514759, 1.2566371,  0.0106841,

    1.7302531,  0.2514759,  1.6018468,  0.5267599,  1.4166715,  0.2750545,

    0.9114273,  0.5267599,  0.7830209,  0.2514759,  1.0966026,  0.2750545,

    2.3585718,  -0.2514759, 2.6679766,  -0.2514759, 2.5132742,  0.0106841,

    2.9868903,  0.2514759,  2.8584838,  0.5267599,  2.6733086,  0.2750545,

    2.1680644,  0.5267599,  2.0396581,  0.2514759,  2.3532395,  0.2750545,

    -2.6679766, -0.2514759, -2.3585718, -0.2514759, -2.5132742, 0.0106841,

    -2.0396581, 0.2514759,  -2.1680644, 0.5267599,  -2.3532395, 0.2750545,

    -2.8584838, 0.5267599,  -2.9868903, 0.2514759,  -2.6733086, 0.2750545,

    -1.4113395, -0.2514759, -1.1019347, -0.2514759, -1.2566371, 0.0106841,

    -0.7830209, 0.2514759,  -0.9114273, 0.5267599,  -1.0966026, 0.2750545,

    -1.6018468, 0.5267599,  -1.7302531, 0.2514759,  -1.4166715, 0.2750545,

    -0.1547024, -0.2514759, 0.1547024,  -0.2514759, 0.0000000,  0.0106841,

    0.4736161,  0.2514759,  0.3452097,  0.5267599,  0.1600345,  0.2750545,

    -0.3452097, 0.5267599,  -0.4736161, 0.2514759,  -0.1600345, 0.2750545,

    2.1106849,  -0.8159013, 1.8849555,  -0.5535744, 1.6592263,  -0.8159014,

    1.6592263,  -0.8159014, 1.2566371,  -1.0172219, 1.8849555,  -1.0964646,

    1.8849555,  -1.0964646, 2.5132742,  -1.0172219, 2.1106849,  -0.8159013,

    0.8540478,  -0.8159013, 0.6283185,  -0.5535743, 0.4025893,  -0.8159013,

    0.4025893,  -0.8159013, 0.0000000,  -1.0172219, 0.6283185,  -1.0964646,

    0.6283185,  -1.0964646, 1.2566371,  -1.0172219, 0.8540478,  -0.8159013,

    -0.4025893, -0.8159013, -0.6283185, -0.5535743, -0.8540478, -0.8159013,

    -0.8540478, -0.8159013, -1.2566371, -1.0172219, -0.6283185, -1.0964646,

    -0.6283185, -1.0964646, 0.0000000,  -1.0172219, -0.4025893, -0.8159013,

    -1.6592263, -0.8159014, -1.8849555, -0.5535744, -2.1106849, -0.8159013,

    -2.1106849, -0.8159013, -2.5132742, -1.0172219, -1.8849555, -1.0964646,

    -1.8849555, -1.0964646, -1.2566371, -1.0172219, -1.6592263, -0.8159014,

    3.1415927,  -1.0964645, -2.5132742, -1.0172219, -2.9158635, -0.8159013,

    -2.9158635, -0.8159013, -3.1415927, -0.5535743, 2.9158635,  -0.8159013,

    2.9158635,  -0.8159013, 2.5132742,  -1.0172219, 3.1415927,  -1.0964645,

    0.7883530,  -0.2750545, 0.9424778,  -0.0000000, 0.6283185,  -0.0106841,

    0.6283185,  -0.0106841, 0.3141593,  -0.0000000, 0.4682840,  -0.2750545,

    0.4682840,  -0.2750545, 0.6283185,  -0.5535743, 0.7883530,  -0.2750545,

    2.0449901,  -0.2750545, 2.1991148,  -0.0000000, 1.8849555,  -0.0106841,

    1.8849555,  -0.0106841, 1.5707964,  -0.0000000, 1.7249211,  -0.2750545,

    1.7249211,  -0.2750545, 1.8849555,  -0.5535744, 2.0449901,  -0.2750545,

    -2.9815581, -0.2750545, -2.8274333, -0.0000000, 3.1415927,  -0.0106841,

    3.1415927,  -0.0106841, 2.8274333,  -0.0000000, 2.9815581,  -0.2750545,

    2.9815581,  -0.2750545, -3.1415927, -0.5535743, -2.9815581, -0.2750545,

    -1.7249211, -0.2750545, -1.5707964, -0.0000000, -1.8849555, -0.0106841,

    -1.8849555, -0.0106841, -2.1991148, -0.0000000, -2.0449901, -0.2750545,

    -2.0449901, -0.2750545, -1.8849555, -0.5535744, -1.7249211, -0.2750545,

    -0.4682840, -0.2750545, -0.3141593, -0.0000000, -0.6283185, -0.0106841,

    -0.6283185, -0.0106841, -0.9424778, -0.0000000, -0.7883530, -0.2750545,

    -0.7883530, -0.2750545, -0.6283185, -0.5535743, -0.4682840, -0.2750545,

    1.0309078,  0.8159013,  1.2566371,  0.5535744,  1.4823663,  0.8159014,

    1.4823663,  0.8159014,  1.8849555,  1.0172219,  1.2566371,  1.0964646,

    1.2566371,  1.0964646,  0.6283185,  1.0172219,  1.0309078,  0.8159013,

    2.2875450,  0.8159013,  2.5132742,  0.5535743,  2.7390034,  0.8159013,

    2.7390034,  0.8159013,  3.1415927,  1.0172219,  2.5132742,  1.0964646,

    2.5132742,  1.0964646,  1.8849555,  1.0172219,  2.2875450,  0.8159013,

    -2.7390034, 0.8159013,  -2.5132742, 0.5535743,  -2.2875450, 0.8159013,

    -2.2875450, 0.8159013,  -1.8849555, 1.0172219,  -2.5132742, 1.0964646,

    -2.5132742, 1.0964646,  3.1415927,  1.0172219,  -2.7390034, 0.8159013,

    -1.4823663, 0.8159014,  -1.2566371, 0.5535744,  -1.0309078, 0.8159013,

    -1.0309078, 0.8159013,  -0.6283185, 1.0172219,  -1.2566371, 1.0964646,

    -1.2566371, 1.0964646,  -1.8849555, 1.0172219,  -1.4823663, 0.8159014,

    0.0000000,  1.0964645,  -0.6283185, 1.0172219,  -0.2257293, 0.8159013,

    -0.2257293, 0.8159013,  -0.0000000, 0.5535743,  0.2257293,  0.8159013,

    0.2257293,  0.8159013,  0.6283185,  1.0172219,  0.0000000,  1.0964645,

    1.2566371,  0.0106841,  1.5707964,  -0.0000000, 1.4166715,  0.2750545,

    1.4166715,  0.2750545,  1.2566371,  0.5535744,  1.0966026,  0.2750545,

    1.0966026,  0.2750545,  0.9424778,  -0.0000000, 1.2566371,  0.0106841,

    2.5132742,  0.0106841,  2.8274333,  -0.0000000, 2.6733086,  0.2750545,

    2.6733086,  0.2750545,  2.5132742,  0.5535743,  2.3532395,  0.2750545,

    2.3532395,  0.2750545,  2.1991148,  -0.0000000, 2.5132742,  0.0106841,

    -2.5132742, 0.0106841,  -2.1991148, -0.0000000, -2.3532395, 0.2750545,

    -2.3532395, 0.2750545,  -2.5132742, 0.5535743,  -2.6733086, 0.2750545,

    -2.6733086, 0.2750545,  -2.8274333, -0.0000000, -2.5132742, 0.0106841,

    -1.2566371, 0.0106841,  -0.9424778, -0.0000000, -1.0966026, 0.2750545,

    -1.0966026, 0.2750545,  -1.2566371, 0.5535744,  -1.4166715, 0.2750545,

    -1.4166715, 0.2750545,  -1.5707964, -0.0000000, -1.2566371, 0.0106841,

    0.0000000,  0.0106841,  0.3141593,  -0.0000000, 0.1600345,  0.2750545,

    0.1600345,  0.2750545,  -0.0000000, 0.5535743,  -0.1600345, 0.2750545,

    -0.1600345, 0.2750545,  -0.3141593, -0.0000000, 0.0000000,  0.0106841,

    2.5132742,  -0.7203593, 2.5132742,  -0.4636476, 2.2301655,  -0.5267599,

    2.2301655,  -0.5267599, 1.8849555,  -0.5535744, 2.1106849,  -0.8159013,

    2.1106849,  -0.8159013, 2.5132742,  -1.0172219, 2.5132742,  -0.7203593,

    1.5397458,  -0.5267599, 1.2566371,  -0.4636476, 1.2566371,  -0.7203593,

    1.2566371,  -0.7203593, 1.2566371,  -1.0172219, 1.6592263,  -0.8159014,

    1.6592263,  -0.8159014, 1.8849555,  -0.5535744, 1.5397458,  -0.5267599,

    1.2566371,  -1.3140846, 3.1415927,  -1.5707964, 2.5132742,  -1.3140846,

    2.5132742,  -1.3140846, 2.5132742,  -1.0172219, 1.8849555,  -1.0964646,

    1.8849555,  -1.0964646, 1.2566371,  -1.0172219, 1.2566371,  -1.3140846,

    1.2566371,  -0.7203593, 1.2566371,  -0.4636476, 0.9735283,  -0.5267599,

    0.9735283,  -0.5267599, 0.6283185,  -0.5535743, 0.8540478,  -0.8159013,

    0.8540478,  -0.8159013, 1.2566371,  -1.0172219, 1.2566371,  -0.7203593,

    0.2831088,  -0.5267599, 0.0000000,  -0.4636476, 0.0000000,  -0.7203593,

    0.0000000,  -0.7203593, 0.0000000,  -1.0172219, 0.4025893,  -0.8159013,

    0.4025893,  -0.8159013, 0.6283185,  -0.5535743, 0.2831088,  -0.5267599,

    0.0000000,  -1.3140846, 3.1415927,  -1.5707964, 1.2566371,  -1.3140846,

    1.2566371,  -1.3140846, 1.2566371,  -1.0172219, 0.6283185,  -1.0964646,

    0.6283185,  -1.0964646, 0.0000000,  -1.0172219, 0.0000000,  -1.3140846,

    0.0000000,  -0.7203593, 0.0000000,  -0.4636476, -0.2831088, -0.5267599,

    -0.2831088, -0.5267599, -0.6283185, -0.5535743, -0.4025893, -0.8159013,

    -0.4025893, -0.8159013, 0.0000000,  -1.0172219, 0.0000000,  -0.7203593,

    -0.9735283, -0.5267599, -1.2566371, -0.4636476, -1.2566371, -0.7203593,

    -1.2566371, -0.7203593, -1.2566371, -1.0172219, -0.8540478, -0.8159013,

    -0.8540478, -0.8159013, -0.6283185, -0.5535743, -0.9735283, -0.5267599,

    -1.2566371, -1.3140846, 3.1415927,  -1.5707964, 0.0000000,  -1.3140846,

    0.0000000,  -1.3140846, 0.0000000,  -1.0172219, -0.6283185, -1.0964646,

    -0.6283185, -1.0964646, -1.2566371, -1.0172219, -1.2566371, -1.3140846,

    -1.2566371, -0.7203593, -1.2566371, -0.4636476, -1.5397458, -0.5267599,

    -1.5397458, -0.5267599, -1.8849555, -0.5535744, -1.6592263, -0.8159014,

    -1.6592263, -0.8159014, -1.2566371, -1.0172219, -1.2566371, -0.7203593,

    -2.2301655, -0.5267599, -2.5132742, -0.4636476, -2.5132742, -0.7203593,

    -2.5132742, -0.7203593, -2.5132742, -1.0172219, -2.1106849, -0.8159013,

    -2.1106849, -0.8159013, -1.8849555, -0.5535744, -2.2301655, -0.5267599,

    -2.5132742, -1.3140846, 3.1415927,  -1.5707964, -1.2566371, -1.3140846,

    -1.2566371, -1.3140846, -1.2566371, -1.0172219, -1.8849555, -1.0964646,

    -1.8849555, -1.0964646, -2.5132742, -1.0172219, -2.5132742, -1.3140846,

    2.5132742,  -1.3140846, 3.1415927,  -1.5707964, -2.5132742, -1.3140846,

    -2.5132742, -1.3140846, -2.5132742, -1.0172219, 3.1415927,  -1.0964645,

    3.1415927,  -1.0964645, 2.5132742,  -1.0172219, 2.5132742,  -1.3140846,

    -2.5132742, -0.7203593, -2.5132742, -0.4636476, -2.7963829, -0.5267599,

    -2.7963829, -0.5267599, -3.1415927, -0.5535743, -2.9158635, -0.8159013,

    -2.9158635, -0.8159013, -2.5132742, -1.0172219, -2.5132742, -0.7203593,

    2.7963829,  -0.5267599, 2.5132742,  -0.4636476, 2.5132742,  -0.7203593,

    2.5132742,  -0.7203593, 2.5132742,  -1.0172219, 2.9158635,  -0.8159013,

    2.9158635,  -0.8159013, -3.1415927, -0.5535743, 2.7963829,  -0.5267599,

    0.9735283,  -0.5267599, 1.2566371,  -0.4636476, 1.1019347,  -0.2514759,

    1.1019347,  -0.2514759, 0.9424778,  -0.0000000, 0.7883530,  -0.2750545,

    0.7883530,  -0.2750545, 0.6283185,  -0.5535743, 0.9735283,  -0.5267599,

    0.7830209,  0.2514759,  0.6283185,  0.4636476,  0.4736161,  0.2514759,

    0.4736161,  0.2514759,  0.3141593,  -0.0000000, 0.6283185,  -0.0106841,

    0.6283185,  -0.0106841, 0.9424778,  -0.0000000, 0.7830209,  0.2514759,

    0.1547024,  -0.2514759, 0.0000000,  -0.4636476, 0.2831088,  -0.5267599,

    0.2831088,  -0.5267599, 0.6283185,  -0.5535743, 0.4682840,  -0.2750545,

    0.4682840,  -0.2750545, 0.3141593,  -0.0000000, 0.1547024,  -0.2514759,

    2.2301655,  -0.5267599, 2.5132742,  -0.4636476, 2.3585718,  -0.2514759,

    2.3585718,  -0.2514759, 2.1991148,  -0.0000000, 2.0449901,  -0.2750545,

    2.0449901,  -0.2750545, 1.8849555,  -0.5535744, 2.2301655,  -0.5267599,

    2.0396581,  0.2514759,  1.8849555,  0.4636476,  1.7302531,  0.2514759,

    1.7302531,  0.2514759,  1.5707964,  -0.0000000, 1.8849555,  -0.0106841,

    1.8849555,  -0.0106841, 2.1991148,  -0.0000000, 2.0396581,  0.2514759,

    1.4113395,  -0.2514759, 1.2566371,  -0.4636476, 1.5397458,  -0.5267599,

    1.5397458,  -0.5267599, 1.8849555,  -0.5535744, 1.7249211,  -0.2750545,

    1.7249211,  -0.2750545, 1.5707964,  -0.0000000, 1.4113395,  -0.2514759,

    -2.7963829, -0.5267599, -2.5132742, -0.4636476, -2.6679766, -0.2514759,

    -2.6679766, -0.2514759, -2.8274333, -0.0000000, -2.9815581, -0.2750545,

    -2.9815581, -0.2750545, -3.1415927, -0.5535743, -2.7963829, -0.5267599,

    -2.9868903, 0.2514759,  3.1415927,  0.4636476,  2.9868903,  0.2514759,

    2.9868903,  0.2514759,  2.8274333,  -0.0000000, 3.1415927,  -0.0106841,

    3.1415927,  -0.0106841, -2.8274333, -0.0000000, -2.9868903, 0.2514759,

    2.6679766,  -0.2514759, 2.5132742,  -0.4636476, 2.7963829,  -0.5267599,

    2.7963829,  -0.5267599, -3.1415927, -0.5535743, 2.9815581,  -0.2750545,

    2.9815581,  -0.2750545, 2.8274333,  -0.0000000, 2.6679766,  -0.2514759,

    -1.5397458, -0.5267599, -1.2566371, -0.4636476, -1.4113395, -0.2514759,

    -1.4113395, -0.2514759, -1.5707964, -0.0000000, -1.7249211, -0.2750545,

    -1.7249211, -0.2750545, -1.8849555, -0.5535744, -1.5397458, -0.5267599,

    -1.7302531, 0.2514759,  -1.8849555, 0.4636476,  -2.0396581, 0.2514759,

    -2.0396581, 0.2514759,  -2.1991148, -0.0000000, -1.8849555, -0.0106841,

    -1.8849555, -0.0106841, -1.5707964, -0.0000000, -1.7302531, 0.2514759,

    -2.3585718, -0.2514759, -2.5132742, -0.4636476, -2.2301655, -0.5267599,

    -2.2301655, -0.5267599, -1.8849555, -0.5535744, -2.0449901, -0.2750545,

    -2.0449901, -0.2750545, -2.1991148, -0.0000000, -2.3585718, -0.2514759,

    -0.2831088, -0.5267599, 0.0000000,  -0.4636476, -0.1547024, -0.2514759,

    -0.1547024, -0.2514759, -0.3141593, -0.0000000, -0.4682840, -0.2750545,

    -0.4682840, -0.2750545, -0.6283185, -0.5535743, -0.2831088, -0.5267599,

    -0.4736161, 0.2514759,  -0.6283185, 0.4636476,  -0.7830209, 0.2514759,

    -0.7830209, 0.2514759,  -0.9424778, -0.0000000, -0.6283185, -0.0106841,

    -0.6283185, -0.0106841, -0.3141593, -0.0000000, -0.4736161, 0.2514759,

    -1.1019347, -0.2514759, -1.2566371, -0.4636476, -0.9735283, -0.5267599,

    -0.9735283, -0.5267599, -0.6283185, -0.5535743, -0.7883530, -0.2750545,

    -0.7883530, -0.2750545, -0.9424778, -0.0000000, -1.1019347, -0.2514759,

    0.6283185,  0.7203593,  0.6283185,  0.4636476,  0.9114273,  0.5267599,

    0.9114273,  0.5267599,  1.2566371,  0.5535744,  1.0309078,  0.8159013,

    1.0309078,  0.8159013,  0.6283185,  1.0172219,  0.6283185,  0.7203593,

    1.6018468,  0.5267599,  1.8849555,  0.4636476,  1.8849555,  0.7203593,

    1.8849555,  0.7203593,  1.8849555,  1.0172219,  1.4823663,  0.8159014,

    1.4823663,  0.8159014,  1.2566371,  0.5535744,  1.6018468,  0.5267599,

    1.8849555,  1.3140846,  0.0000000,  1.5707964,  0.6283185,  1.3140846,

    0.6283185,  1.3140846,  0.6283185,  1.0172219,  1.2566371,  1.0964646,

    1.2566371,  1.0964646,  1.8849555,  1.0172219,  1.8849555,  1.3140846,

    1.8849555,  0.7203593,  1.8849555,  0.4636476,  2.1680644,  0.5267599,

    2.1680644,  0.5267599,  2.5132742,  0.5535743,  2.2875450,  0.8159013,

    2.2875450,  0.8159013,  1.8849555,  1.0172219,  1.8849555,  0.7203593,

    2.8584838,  0.5267599,  3.1415927,  0.4636476,  3.1415927,  0.7203593,

    3.1415927,  0.7203593,  3.1415927,  1.0172219,  2.7390034,  0.8159013,

    2.7390034,  0.8159013,  2.5132742,  0.5535743,  2.8584838,  0.5267599,

    3.1415927,  1.3140846,  0.0000000,  1.5707964,  1.8849555,  1.3140846,

    1.8849555,  1.3140846,  1.8849555,  1.0172219,  2.5132742,  1.0964646,

    2.5132742,  1.0964646,  3.1415927,  1.0172219,  3.1415927,  1.3140846,

    3.1415927,  0.7203593,  3.1415927,  0.4636476,  -2.8584838, 0.5267599,

    -2.8584838, 0.5267599,  -2.5132742, 0.5535743,  -2.7390034, 0.8159013,

    -2.7390034, 0.8159013,  3.1415927,  1.0172219,  3.1415927,  0.7203593,

    -2.1680644, 0.5267599,  -1.8849555, 0.4636476,  -1.8849555, 0.7203593,

    -1.8849555, 0.7203593,  -1.8849555, 1.0172219,  -2.2875450, 0.8159013,

    -2.2875450, 0.8159013,  -2.5132742, 0.5535743,  -2.1680644, 0.5267599,

    -1.8849555, 1.3140846,  0.0000000,  1.5707964,  3.1415927,  1.3140846,

    3.1415927,  1.3140846,  3.1415927,  1.0172219,  -2.5132742, 1.0964646,

    -2.5132742, 1.0964646,  -1.8849555, 1.0172219,  -1.8849555, 1.3140846,

    -1.8849555, 0.7203593,  -1.8849555, 0.4636476,  -1.6018468, 0.5267599,

    -1.6018468, 0.5267599,  -1.2566371, 0.5535744,  -1.4823663, 0.8159014,

    -1.4823663, 0.8159014,  -1.8849555, 1.0172219,  -1.8849555, 0.7203593,

    -0.9114273, 0.5267599,  -0.6283185, 0.4636476,  -0.6283185, 0.7203593,

    -0.6283185, 0.7203593,  -0.6283185, 1.0172219,  -1.0309078, 0.8159013,

    -1.0309078, 0.8159013,  -1.2566371, 0.5535744,  -0.9114273, 0.5267599,

    -0.6283185, 1.3140846,  0.0000000,  1.5707964,  -1.8849555, 1.3140846,

    -1.8849555, 1.3140846,  -1.8849555, 1.0172219,  -1.2566371, 1.0964646,

    -1.2566371, 1.0964646,  -0.6283185, 1.0172219,  -0.6283185, 1.3140846,

    0.6283185,  1.3140846,  0.0000000,  1.5707964,  -0.6283185, 1.3140846,

    -0.6283185, 1.3140846,  -0.6283185, 1.0172219,  0.0000000,  1.0964645,

    0.0000000,  1.0964645,  0.6283185,  1.0172219,  0.6283185,  1.3140846,

    -0.6283185, 0.7203593,  -0.6283185, 0.4636476,  -0.3452097, 0.5267599,

    -0.3452097, 0.5267599,  -0.0000000, 0.5535743,  -0.2257293, 0.8159013,

    -0.2257293, 0.8159013,  -0.6283185, 1.0172219,  -0.6283185, 0.7203593,

    0.3452097,  0.5267599,  0.6283185,  0.4636476,  0.6283185,  0.7203593,

    0.6283185,  0.7203593,  0.6283185,  1.0172219,  0.2257293,  0.8159013,

    0.2257293,  0.8159013,  -0.0000000, 0.5535743,  0.3452097,  0.5267599,

    1.1019347,  -0.2514759, 1.2566371,  -0.4636476, 1.4113395,  -0.2514759,

    1.4113395,  -0.2514759, 1.5707964,  -0.0000000, 1.2566371,  0.0106841,

    1.2566371,  0.0106841,  0.9424778,  -0.0000000, 1.1019347,  -0.2514759,

    1.7302531,  0.2514759,  1.8849555,  0.4636476,  1.6018468,  0.5267599,

    1.6018468,  0.5267599,  1.2566371,  0.5535744,  1.4166715,  0.2750545,

    1.4166715,  0.2750545,  1.5707964,  -0.0000000, 1.7302531,  0.2514759,

    0.9114273,  0.5267599,  0.6283185,  0.4636476,  0.7830209,  0.2514759,

    0.7830209,  0.2514759,  0.9424778,  -0.0000000, 1.0966026,  0.2750545,

    1.0966026,  0.2750545,  1.2566371,  0.5535744,  0.9114273,  0.5267599,

    2.3585718,  -0.2514759, 2.5132742,  -0.4636476, 2.6679766,  -0.2514759,

    2.6679766,  -0.2514759, 2.8274333,  -0.0000000, 2.5132742,  0.0106841,

    2.5132742,  0.0106841,  2.1991148,  -0.0000000, 2.3585718,  -0.2514759,

    2.9868903,  0.2514759,  3.1415927,  0.4636476,  2.8584838,  0.5267599,

    2.8584838,  0.5267599,  2.5132742,  0.5535743,  2.6733086,  0.2750545,

    2.6733086,  0.2750545,  2.8274333,  -0.0000000, 2.9868903,  0.2514759,

    2.1680644,  0.5267599,  1.8849555,  0.4636476,  2.0396581,  0.2514759,

    2.0396581,  0.2514759,  2.1991148,  -0.0000000, 2.3532395,  0.2750545,

    2.3532395,  0.2750545,  2.5132742,  0.5535743,  2.1680644,  0.5267599,

    -2.6679766, -0.2514759, -2.5132742, -0.4636476, -2.3585718, -0.2514759,

    -2.3585718, -0.2514759, -2.1991148, -0.0000000, -2.5132742, 0.0106841,

    -2.5132742, 0.0106841,  -2.8274333, -0.0000000, -2.6679766, -0.2514759,

    -2.0396581, 0.2514759,  -1.8849555, 0.4636476,  -2.1680644, 0.5267599,

    -2.1680644, 0.5267599,  -2.5132742, 0.5535743,  -2.3532395, 0.2750545,

    -2.3532395, 0.2750545,  -2.1991148, -0.0000000, -2.0396581, 0.2514759,

    -2.8584838, 0.5267599,  3.1415927,  0.4636476,  -2.9868903, 0.2514759,

    -2.9868903, 0.2514759,  -2.8274333, -0.0000000, -2.6733086, 0.2750545,

    -2.6733086, 0.2750545,  -2.5132742, 0.5535743,  -2.8584838, 0.5267599,

    -1.4113395, -0.2514759, -1.2566371, -0.4636476, -1.1019347, -0.2514759,

    -1.1019347, -0.2514759, -0.9424778, -0.0000000, -1.2566371, 0.0106841,

    -1.2566371, 0.0106841,  -1.5707964, -0.0000000, -1.4113395, -0.2514759,

    -0.7830209, 0.2514759,  -0.6283185, 0.4636476,  -0.9114273, 0.5267599,

    -0.9114273, 0.5267599,  -1.2566371, 0.5535744,  -1.0966026, 0.2750545,

    -1.0966026, 0.2750545,  -0.9424778, -0.0000000, -0.7830209, 0.2514759,

    -1.6018468, 0.5267599,  -1.8849555, 0.4636476,  -1.7302531, 0.2514759,

    -1.7302531, 0.2514759,  -1.5707964, -0.0000000, -1.4166715, 0.2750545,

    -1.4166715, 0.2750545,  -1.2566371, 0.5535744,  -1.6018468, 0.5267599,

    -0.1547024, -0.2514759, 0.0000000,  -0.4636476, 0.1547024,  -0.2514759,

    0.1547024,  -0.2514759, 0.3141593,  -0.0000000, 0.0000000,  0.0106841,

    0.0000000,  0.0106841,  -0.3141593, -0.0000000, -0.1547024, -0.2514759,

    0.4736161,  0.2514759,  0.6283185,  0.4636476,  0.3452097,  0.5267599,

    0.3452097,  0.5267599,  -0.0000000, 0.5535743,  0.1600345,  0.2750545,

    0.1600345,  0.2750545,  0.3141593,  -0.0000000, 0.4736161,  0.2514759,

    -0.3452097, 0.5267599,  -0.6283185, 0.4636476,  -0.4736161, 0.2514759,

    -0.4736161, 0.2514759,  -0.3141593, -0.0000000, -0.1600345, 0.2750545,

    -0.1600345, 0.2750545,  -0.0000000, 0.5535743,  -0.3452097, 0.5267599};

}  // namespace spherical
}  // namespace openMVG