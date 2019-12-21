// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2011, 2012, 2013, 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SVG_HISTOGRAM_HPP
#define SVG_HISTOGRAM_HPP

#include "svgDrawer.hpp"

#include <fstream>
#include <vector>

namespace svg {

/// Helper function to draw a SVG histogram
/// ____
/// |  |   ___ |
/// |  |__|  | |
/// |  |  |  | |
/// -----------|

template<typename T>
static std::string stringifier(const T & t)
{
  std::ostringstream os;
  os << t;
  return os.str();
}

template<typename T>
bool drawHistogram(const std::vector<T> & vec_value,
  const std::pair<float, float> & range,
  const std::string & sFilename,
  const float W, const float H)
{
  if (vec_value.empty()) {
    return false;
  }
  //-- Max value
  const T maxi = *max_element(vec_value.cbegin(), vec_value.cend());
  const size_t n = vec_value.size();

  const float scaleFactorY = H / static_cast<float>(maxi);
  const float scaleFactorX = W / static_cast<float>(n);

  svgDrawer svg_stream;

  for (size_t i = 0; i < vec_value.size(); ++i)
  {
    const size_t dist = i;
    const T val = vec_value[i];
    std::ostringstream os;
    os << '(' << range.first + dist/float(n) * (range.second-range.first) << ',' << val << ')';
    const svgAttributes style = svgAttributes().fill("blue").stroke("black", 1.0).tooltip(os.str());
    svg_stream << drawRectangle(
      scaleFactorX * dist, H-val * scaleFactorY,
      scaleFactorX, val * scaleFactorY,
      style);
    //_________
    //|       |_________
    //|       ||       |
    //|       ||       |
    //|       ||       |
    //0    sFactorX  2*sFactorX
  }
  const svgAttributes axis_style = svgAttributes().stroke("black", 1.0f);
  // Draw X Axis
  svg_stream << drawText(.05f*W, 1.2f*H, .1f*H, stringifier(range.first), "black");
  svg_stream << drawText(   W, 1.2f*H, .1f*H, stringifier(range.second), "black");
  svg_stream << drawLine(0, 1.1f*H, W, 1.1f*H, axis_style);
  // Draw Y Axis
  svg_stream << drawText(1.2f*W, .1f*H, .1f*H, stringifier(maxi), "black");
  svg_stream << drawText(1.2f*W, H, .1f*H, "0", "black");
  svg_stream << drawLine(1.1f*W, 0, 1.1f*W, H, axis_style);

  std::ofstream svg_file_stream( sFilename.c_str());
  if (svg_file_stream)
  {
    svg_file_stream << svg_stream.closeSvgFile().str();
    return true;
  }
  return false;
}

} // namespace svg

#endif // SVG_HISTOGRAM_HPP
