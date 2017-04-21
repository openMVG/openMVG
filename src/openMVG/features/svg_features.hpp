// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_SVG_FEATURES_HPP
#define OPENMVG_FEATURES_SVG_FEATURES_HPP

#include <string>
#include <utility>

#include <openMVG/features/feature.hpp>
#include <openMVG/features/feature_container.hpp>

namespace openMVG {
namespace features {

/**
 * @brief Saves a svg file containing the image and its features
 *
 * @param[in] image_path Left image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] image_size Image size {width,height}.
 * @param[in] features image features.
 * @param[in] svg_filename Path to the output SVG file.
 * @param[in] feature_circle_radius Radius of the circle used to represent
 * the features.
 * @param[in] stroke_size Stroke size used to display the line between the
 * corresponding features.
 */
bool Features2SVG
(
  const std::string & image_path,
  const std::pair<size_t,size_t> & image_size,
  const features::PointFeatures & features,
  const std::string & svg_filename,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
);

/**
 * @brief Saves a svg file containing the image and its scale invariant features
 *
 * @param[in] image_path Left image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] image_size Image size {width,height}.
 * @param[in] features image features.
 * @param[in] svg_filename Path to the output SVG file.
 * @param[in] feature_circle_radius Radius of the circle used to represent
 * the features.
 * @param[in] stroke_size Stroke size used to display the line between the
 * corresponding features.
 */
bool Features2SVG
(
  const std::string & image_path,
  const std::pair<size_t,size_t> & image_size,
  const features::SIOPointFeatures & features,
  const std::string & svg_filename,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
);

/**
 * @brief Saves a svg file containing two images and their features
 *
 * @param[in] left_image_path Left image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] left_image_size Left image size {width,height}.
 * @param[in] left_features Left image features.
 * @param[in] right_image_path Right image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] right_image_size Right image size {width,height}.
 * @param[in] right_features Right image features.
 * @param[in] b_vertical_display Display images in a vertical or horizontal
 * canvas.
 * @param[in] svg_filename Path to the output SVG file.
 * @param[in] feature_circle_radius Radius of the circle used to represent
 * the features.
 * @param[in] stroke_size Stroke size used to display the line between the
 * corresponding features.
 */
bool Features2SVG
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const std::string & svg_filename,
  const bool b_vertical_display = true,
  const double feature_circle_radius = 4.0,
  const double stroke_size = 2.0
);

/**
 * @brief Saves a svg file containing two images and their features
 *
 * @param[in] left_image_path Left image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] left_image_size Left image size {width,height}.
 * @param[in] left_features Left image features.
 * @param[in] right_image_path Right image path. For compactness of the output
 * svg file the image is exported as a link.
 * @param[in] right_image_size Right image size {width,height}.
 * @param[in] right_features Right image features.
 * @param[in] b_vertical_display Display images in a vertical or horizontal
 * canvas.
 * @param[in] svg_filename Path to the output SVG file.
 * @param[in] stroke_size Stroke size used to display the line between the
 * corresponding features.
 */
bool Features2SVG
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::SIOPointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::SIOPointFeatures & right_features,
  const std::string & svg_filename,
  const bool b_vertical_display = true,
  const double stroke_size = 2.0
);

}  // namespace features
}  // namespace openMVG

#endif // OPENMVG_FEATURES_SVG_FEATURES_HPP
