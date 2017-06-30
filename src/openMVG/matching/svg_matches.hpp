// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_SVG_MATCHES_HPP
#define OPENMVG_MATCHING_SVG_MATCHES_HPP

#include <string>
#include <utility>
#include <vector>

#include <openMVG/matching/indMatch.hpp>
#include <openMVG/features/feature_container.hpp>

namespace openMVG {
namespace matching {

  /**
   * @brief Return a string containing the svg file content to display
   * a scene containing two images and their feature matches:
   * image are exported side by side, feature depicted by circles and
   * corresponding features are connected by a line.
   *
   * @param[in] left_image_path Left image path. For compactness of the output
   * svg file the image is exported as a link.
   * @param[in] left_image_size Left image size {width,height}.
   * @param[in] left_features Left image features.
   * @param[in] right_image_path Right image path. For compactness of the output
   * svg file the image is exported as a link.
   * @param[in] right_image_size Right image size {width,height}.
   * @param[in] right_features Right image features.
   * @param[in] matches Corresponding features indexes.
   * @param[in] b_vertical_display Display images in a vertical or horizontal
   * canvas.
   * @param[in] feature_circle_radius Radius of the circle used to represent
   * the features.
   * @param[in] stroke_size Stroke size used to display the line between the
   * corresponding features.
   */
std::string Matches2SVGString
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const matching::IndMatches & matches,
  const bool b_vertical_display = true,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
);

/**
 * @brief Saves a svg file containing two images and their feature matches:
 * image are exported side by side, feature depicted by circles and
 * corresponding features are connected by a line.
 *
 * @param[in] left_image_path Left image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] left_image_size Left image size {width,height}.
 * @param[in] left_features Left image features.
 * @param[in] right_image_path Right image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] right_image_size Right image size {width,height}.
 * @param[in] right_features Right image features.
 * @param[in] matches Corresponding features indexes.
 * @param[in] b_vertical_display Display images in a vertical or horizontal
 * canvas.
 * @param[in] svg_filename Path to the output SVG file.
 * @param[in] feature_circle_radius Radius of the circle used to represent
 * the features.
 * @param[in] stroke_size Stroke size used to display the line between the
 * corresponding features.
 */
bool Matches2SVG
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const matching::IndMatches & matches,
  const std::string & svg_filename,
  const bool b_vertical_display = true,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
);

/**
 * @brief Saves a svg file containing two images and their inlier features
 * matches: image are exported side by side, feature depicted by circles and
 * corresponding features (inliers) are connected by a line.
 *
 * @param[in] left_image_path Left image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] left_image_size Left image size {width,height}.
 * @param[in] left_features Left image features.
 * @param[in] right_image_path Right image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] right_image_size Right image size {width,height}.
 * @param[in] right_features Right image features.
 * @param[in] matches Corresponding features indexes.
 * @param[in] inliers Index of the matches array that corresponds to inliers.
 * @param[in] b_vertical_display Display images in a vertical or horizontal
 * canvas.
 * @param[in] svg_filename Path to the output SVG file.
 * @param[in] feature_circle_radius Radius of the circle used to represent
 * the features.
 * @param[in] stroke_size Stroke size used to display the line between the
 * corresponding features.
 */
bool InlierMatches2SVG
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const matching::IndMatches & matches,
  const std::vector<uint32_t> & inliers,
  const std::string & svg_filename,
  const bool b_vertical_display = true,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
);

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_SVG_MATCHES_HPP
