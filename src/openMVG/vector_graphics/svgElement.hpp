// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SVG_ELEMENT_HPP
#define SVG_ELEMENT_HPP

#include <iterator>
#include <sstream>
#include <string>

namespace svg {

//
// Collection of functions to add ability to describe SVG element
// see https://developer.mozilla.org/en-US/docs/Web/SVG/Element
//

///Circle draw -> x,y position and radius
static std::string drawCircle(float cx, float cy, float r, const svgAttributes & style)
{
  std::ostringstream os;
  os  << "<circle cx=\"" << cx << "\"" << " cy=\"" << cy << "\""
    << " r=\"" << r << "\""
    << style.getSvgStream() + (style.bTooltip() ? "</circle>" : "/>\n");
  return os.str();
}

///Line draw -> start and end point
static std::string drawLine(float ax, float ay, float bx, float by, const svgAttributes & style)
{
  std::ostringstream os;
  os <<
    "<polyline points=\""<< ax << "," << ay << "," << bx << "," << by <<"\""
    << style.getSvgStream() +  (style.bTooltip() ? "</polyline>" : "/>\n");
  return os.str();
}

///Reference to an image (path must be relative to the svg file)
static std::string drawImage(const std::string & simagePath, int W, int H,
  int posx = 0, int posy = 0, float opacity =1.f)
{
  std::ostringstream os;
  os <<
    "<image x=\""<< posx << "\"" << " y=\""<< posy << "\""
    << " width=\""<< W << "px\"" << " height=\""<< H << "px\""
    << " opacity=\""<< opacity << "\""
    << " xlink:href=\"" << simagePath << "\"" << "/>\n";
  return os.str();
}

///Circle draw -> x,y position and width and height
static std::string drawRectangle(float cx, float cy, float W, float H,
  const svgAttributes & style)
{
  std::ostringstream os;
  os << "<rect x=\"" << cx << "\""
     << " y=\"" << cy << "\""
     << " width=\"" << W << "\""
     << " height=\"" << H << "\""
     << style.getSvgStream() + (style.bTooltip() ? "</rect>" : "/>\n");
  return os.str();
}

///Square draw -> x,y position and size
static std::string drawSquare(float cx, float cy, float W,
  const svgAttributes & style)
{
  return drawRectangle(cx, cy, W, W, style);
}

///Text display -> x,y position, font size
static std::string drawText(float cx, float cy, float fontSize = 1.0f,
  const std::string & stext = "", const std::string & scol = "")
{
  std::ostringstream os;
  os << "<text" << " x=\"" << cx << "\"" << " y=\"" << cy << "\""
    << " font-size=\"" << fontSize << "\"";
  if (!scol.empty())
    os << " fill=\"" << scol << "\"";

  os << ">" << stext << "</text>\n";
  return os.str();
}

// Draw a polyline for point data saved in non contiguous format
// X and Y data are saved in two different container
template< typename DataInputIteratorX, typename DataInputIteratorY>
static std::string drawPolyline(DataInputIteratorX xStart, DataInputIteratorX xEnd,
   DataInputIteratorY yStart, DataInputIteratorY yEnd,
  const svgAttributes & style)
{
  std::ostringstream os;
  os << "<polyline points=\"";

  DataInputIteratorY itery = yStart;
  for (DataInputIteratorX iterx = xStart;
    iterx != xEnd;
    std::advance(iterx, 1), std::advance(itery, 1))
  {
    os << *iterx << ',' << *itery << ' ';
  }
  os << "\"" << style.getSvgStream() + (style.bTooltip() ? "</polyline>" : "/>\n");
  return os.str();
}

// Draw a polyline for point data saved in contiguous format (x0,y0,x1,y1,...)
template< typename DataInputIteratorXY>
static std::string drawPolyline(DataInputIteratorXY points,
  const svgAttributes & style)
{
  std::ostringstream os;
  os << "<polyline points=\"";

  std::copy(points.cbegin(), points.cend(),
    std::ostream_iterator<typename DataInputIteratorXY::value_type>(os, ","));

  os << "\""
    << style.getSvgStream() + (style.bTooltip() ? "</polyline>" : "/>\n");
  return os.str();
}

} // namespace svg

#endif // SVG_ELEMENT_HPP
