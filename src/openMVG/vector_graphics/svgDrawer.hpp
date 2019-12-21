// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2011, 2012, 2013, 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SVG_DRAWER_HPP
#define SVG_DRAWER_HPP

#include "svgAttributes.hpp"
#include "svgElement.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace svg {

//
// Helper class to stream SVG element to a stream/string
// See https://developer.mozilla.org/en-US/docs/Web/SVG
//

/// Class to handle concatenation of SVG elements and export of the XML content to a file.
/// Classic usage:
/// svgDrawer svg_surface;
/// svg_surface << drawCircle(1, 1, 0.5);
//  Then export "svgDrawer.closeSvgFile().str()" to a .svg file.
class svgDrawer{
public:
  ///Constructor
  svgDrawer(const size_t width = 0, const size_t height = 0)
  {
    svg_stream
      << "<?xml version=\"1.0\" standalone=\"yes\"?>\n"
      << "<!-- SVG graphic -->" << "\n"
      << "<svg xmlns='http://www.w3.org/2000/svg'"
      << " xmlns:xlink='http://www.w3.org/1999/xlink'" << "\n";

    if (width > 0 && height > 0)
      svg_stream
        <<"width=\"" << width << "px\" height=\""<< height << "px\""
        << " preserveAspectRatio=\"xMinYMin meet\""
        << " viewBox=\"0 0 " << width << ' ' << height <<"\"";

    svg_stream <<" version=\"1.1\">\n";
  }

  svgDrawer &operator << (const std::string & svg_element)
  {
    svg_stream << svg_element;
    return *this;
  }

  ///Close the svg tag.
  std::ostringstream & closeSvgFile()
  {
    svg_stream <<"</svg>";
    return svg_stream;
  }
private:
  std::ostringstream svg_stream;
};

} // namespace svg

#endif // SVG_DRAWER_HPP
