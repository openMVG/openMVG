// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SVG_ATTRIBUTES_HPP
#define SVG_ATTRIBUTES_HPP

#include <sstream>
#include <string>

namespace svg {

//
// Helper class to handle some draw SVG attributes.
// see https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute
//

/// Class to handle SVG some draw attributes
/// Note: Attributes define how elements should be rendered.
class svgAttributes
{
  std::string fill_color_ = "";
  std::string stroke_color_ = "black";
  std::string tooltip_ = "";
  float stroke_width_ = 1.0f;
  float opacity_ = 1.0f; // Must be in the range [0; 1]
public:

  svgAttributes() = default;

  // Configure fill color
  svgAttributes & fill(const std::string & sFillCol)
  { fill_color_ = sFillCol; return * this;}

  // Configure stroke color and width
  svgAttributes & stroke(const std::string & sStrokeCol, const float strokeWitdh = 1.f)
  { stroke_color_ = sStrokeCol;  stroke_width_ = strokeWitdh; return * this;}

  // Configure with no stroke
  svgAttributes & noStroke()
  { stroke_color_ = "";  stroke_width_ = 0.f; return * this;}

  // Configure tooltip
  svgAttributes & tooltip(const std::string & sTooltip)
  { tooltip_ = sTooltip; return * this;}

  svgAttributes & opacity(const float & opacity)
  { opacity_ = opacity; return *this; }

  const std::string getSvgStream() const {
    std::ostringstream os;

    if (!stroke_color_.empty())
      os << " stroke=\"" << stroke_color_ << "\" stroke-width=\"" << stroke_width_ << "\"";
    if (fill_color_.empty())
      os << " fill=\"none\"";
    else
      os << " fill=\"" << fill_color_ << "\"";
    if (opacity_ > 0)
      os << " opacity=\"" << opacity_ << "\"";
    if (!tooltip_.empty())
      os << " tooltip=\"enable\">"
      << "<title>" << tooltip_ << "</title>";

    return os.str();
  }

  bool bTooltip() const { return !tooltip_.empty();}
};

} // namespace svg

#endif // SVG_ATTRIBUTES_HPP
