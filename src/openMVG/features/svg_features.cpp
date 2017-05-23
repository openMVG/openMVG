// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/features/svg_features.hpp>
#include "third_party/vectorGraphics/svgDrawer.hpp"

namespace openMVG {
namespace features {

bool Features2SVG
(
  const std::string & image_path,
  const std::pair<size_t,size_t> & image_size,
  const features::PointFeatures & features,
  const std::string & svg_filename,
  const double feature_circle_radius,
  const double stroke_size
)
{
  svg::svgDrawer svgStream(image_size.first, image_size.second);

  // Draw image
  svgStream.drawImage(image_path, image_size.first, image_size.second);

  // Draw features
  for (const features::PointFeature & feat_it : features) {
    // Draw the feature (circle)
    svgStream.drawCircle(
      feat_it.x(), feat_it.y(), feature_circle_radius,
      svg::svgStyle().stroke("yellow", stroke_size));
  }

  // Save the SVG file
  std::ofstream svgFile( svg_filename.c_str() );
  if (svgFile.is_open())
  {
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
    return true;
  }
  return false;
}

bool Features2SVG
(
  const std::string & image_path,
  const std::pair<size_t,size_t> & image_size,
  const features::SIOPointFeatures & features,
  const std::string & svg_filename,
  const double stroke_size
)
{
  svg::svgDrawer svgStream(image_size.first, image_size.second);

  // Draw image
  svgStream.drawImage(image_path, image_size.first, image_size.second);

  // Draw features
  for (const features::SIOPointFeature & feat_it : features) {
    // Draw the feature (circle)
    svgStream.drawCircle(
      feat_it.x(), feat_it.y(), feat_it.scale(),
      svg::svgStyle().stroke("yellow", stroke_size));
  }

  // Save the SVG file
  std::ofstream svgFile( svg_filename.c_str() );
  if (svgFile.is_open())
  {
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
    return true;
  }
  return false;
}

bool Features2SVG
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::PointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::PointFeatures & right_features,
  const std::string & svg_filename,
  const bool b_vertical_display,
  const double feature_circle_radius,
  const double stroke_size
)
{
  const size_t svg_w =
    b_vertical_display ?
    left_image_size.first :
    left_image_size.first + right_image_size.first;
  const size_t svg_h =
    b_vertical_display ?
    left_image_size.second + right_image_size.second :
    left_image_size.second;
  const size_t svg_offset_x =
    b_vertical_display ?
    0 :
    left_image_size.first;
  const size_t svg_offset_y =
    b_vertical_display ?
    left_image_size.second :
    0;

  svg::svgDrawer svgStream(svg_w, svg_h);

  // Draw image side by side
  svgStream.drawImage(left_image_path, left_image_size.first, left_image_size.second);
  svgStream.drawImage(right_image_path, right_image_size.first, right_image_size.second,
    b_vertical_display ? 0 : left_image_size.first,
    b_vertical_display ? left_image_size.second : 0);

  // Display feature circles
  for (const features::PointFeature & feat_it : left_features) {
    svgStream.drawCircle(
      feat_it.x(), feat_it.y(), feature_circle_radius,
      svg::svgStyle().stroke("yellow", stroke_size));
  }
  for (const features::PointFeature & feat_it : right_features) {
    svgStream.drawCircle(
      feat_it.x() + svg_offset_x, feat_it.y() + svg_offset_y, feature_circle_radius,
      svg::svgStyle().stroke("yellow", stroke_size));
  }

  // Save the SVG file
  std::ofstream svgFile( svg_filename.c_str() );
  if (svgFile.is_open())
  {
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
    return true;
  }
  return false;
}

bool Features2SVG
(
  const std::string & left_image_path,
  const std::pair<size_t,size_t> & left_image_size,
  const features::SIOPointFeatures & left_features,
  const std::string & right_image_path,
  const std::pair<size_t,size_t> & right_image_size,
  const features::SIOPointFeatures & right_features,
  const std::string & svg_filename,
  const bool b_vertical_display,
  const double stroke_size
)
{
  const size_t svg_w =
    b_vertical_display ?
    left_image_size.first :
    left_image_size.first + right_image_size.first;
  const size_t svg_h =
    b_vertical_display ?
    left_image_size.second + right_image_size.second :
    left_image_size.second;
  const size_t svg_offset_x =
    b_vertical_display ?
    0 :
    left_image_size.first;
  const size_t svg_offset_y =
    b_vertical_display ?
    left_image_size.second :
    0;

  svg::svgDrawer svgStream(svg_w, svg_h);

  // Draw image side by side
  svgStream.drawImage(left_image_path, left_image_size.first, left_image_size.second);
  svgStream.drawImage(right_image_path, right_image_size.first, right_image_size.second,
    b_vertical_display ? 0 : left_image_size.first,
    b_vertical_display ? left_image_size.second : 0);

  // Display feature circles
  for (const features::SIOPointFeature & feat_it : left_features) {
    svgStream.drawCircle(
      feat_it.x(), feat_it.y(), feat_it.scale(),
      svg::svgStyle().stroke("yellow", stroke_size));
  }
  for (const features::SIOPointFeature & feat_it : right_features) {
    svgStream.drawCircle(
      feat_it.x() + svg_offset_x, feat_it.y() + svg_offset_y, feat_it.scale(),
      svg::svgStyle().stroke("yellow", stroke_size));
  }

  // Save the SVG file
  std::ofstream svgFile( svg_filename.c_str() );
  if (svgFile.is_open())
  {
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
    return true;
  }
  return false;
}


}  // namespace features
}  // namespace openMVG
