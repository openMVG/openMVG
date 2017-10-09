
/**
 * @file svgDrawer.h
 * @brief Simple svg drawing from C++
 * @author Pierre MOULON
 *
 * Copyright (c) 2011, 2012, 2013 Pierre MOULON
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef THE_FRENCH_LEAF_SVG_DRAWER_H
#define THE_FRENCH_LEAF_SVG_DRAWER_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace svg {

/// Basic svg style
class svgStyle
{
  std::string _sFillCol, _sStrokeCol, _sToolTip;
  float _strokeW;
public:

  svgStyle():_sFillCol(""), _sStrokeCol("black"), _sToolTip(""),
    _strokeW(1.0f) {}

  // Configure fill color
  svgStyle & fill(const std::string & sFillCol)
  { _sFillCol = sFillCol; return * this;}

  // Configure stroke color and width
  svgStyle & stroke(const std::string & sStrokeCol, float strokeWitdh = 1.f)
  { _sStrokeCol = sStrokeCol;  _strokeW = strokeWitdh; return * this;}

  // Configure with no stroke
  svgStyle & noStroke()
  { _sStrokeCol = "";  _strokeW = 0.f; return * this;}

  // Configure tooltip
  svgStyle & tooltip(const std::string & sTooltip)
  { _sToolTip = sTooltip; return * this;}

  const std::string getSvgStream() const {
    std::ostringstream os;

    if (!_sStrokeCol.empty())
      os << " stroke=\"" << _sStrokeCol << "\" stroke-width=\"" << _strokeW << "\"";
    if (_sFillCol.empty())
      os << " fill=\"none\"";
    else
      os << " fill=\"" << _sFillCol << "\"";
    if (!_sToolTip.empty())
      os << " tooltip=\"enable\">"
      << "<title>" << _sToolTip << "</title>";

    return os.str();
  }

  bool bTooltip() const { return !_sToolTip.empty();}
};


/// Basic class to handle simple SVG draw.
/// You can draw line, square, rectangle, text and image (xlink)
class svgDrawer{
public:
  ///Constructor
  svgDrawer(size_t W = 0, size_t H = 0)
  {
    svgStream << "<?xml version=\"1.0\" standalone=\"yes\"?>\n";

    svgStream << "<!-- SVG graphic -->" << std::endl
      << "<svg xmlns='http://www.w3.org/2000/svg'"
      << " xmlns:xlink='http://www.w3.org/1999/xlink'" << "\n";

    if (W > 0 && H > 0)
      svgStream <<"width=\"" << W << "px\" height=\""<< H << "px\""
      << " preserveAspectRatio=\"xMinYMin meet\""
      << " viewBox=\"0 0 " << W << ' ' << H <<"\"";

    svgStream <<" version=\"1.1\">" << std::endl;
  }
  ///Circle draw -> x,y position and radius
  void drawCircle(float cx, float cy, float r,
    const svgStyle & style)
  {
    svgStream << "<circle cx=\"" << cx << "\"" << " cy=\"" << cy << "\""
      << " r=\"" << r << "\""
      << style.getSvgStream() + (style.bTooltip() ? "</circle>" : "/>\n");  }

  ///Line draw -> start and end point
  void drawLine(float ax, float ay, float bx, float by,
    const svgStyle & style)
  {
    svgStream <<
      "<polyline points=\""<< ax << "," << ay << "," << bx << "," << by <<"\""
      << style.getSvgStream() +  (style.bTooltip() ? "</polyline>" : "/>\n");
  }

  ///Reference to an image (path must be relative to the svg file)
  void drawImage(const std::string & simagePath, int W, int H,
    int posx = 0, int posy = 0, float opacity =1.f)
  {
    svgStream <<
      "<image x=\""<< posx << "\"" << " y=\""<< posy << "\""
      << " width=\""<< W << "px\"" << " height=\""<< H << "px\""
      << " opacity=\""<< opacity << "\""
      << " xlink:href=\"" << simagePath << "\"" << "/>\n";
  }

  ///Square draw -> x,y position and size
  void drawSquare(float cx, float cy, float W,
    const svgStyle & style)
  {
    drawRectangle(cx, cy, W, W, style);
  }

  ///Circle draw -> x,y position and width and height
  void drawRectangle(float cx, float cy, float W, float H,
    const svgStyle & style)
  {
    svgStream << "<rect x=\"" << cx << "\""
      << " y=\"" << cy << "\""
      << " width=\"" << W << "\""
      << " height=\"" << H << "\""
      << style.getSvgStream() + (style.bTooltip() ? "</rect>" : "/>\n");
  }

  ///Text display -> x,y position, font size
  void drawText(float cx, float cy, float fontSize = 1.0f, const std::string & stext ="",
    const std::string & scol = "")
  {
    svgStream << "<text" << " x=\"" << cx << "\"" << " y=\"" << cy << "\""
      << " font-size=\"" << fontSize << "\"";
    if (!scol.empty())
      svgStream << " fill=\"" << scol << "\"";

    svgStream << ">" << stext << "</text>\n";
  }

  template< typename DataInputIteratorX, typename DataInputIteratorY>
  void drawPolyline(DataInputIteratorX xStart, DataInputIteratorX xEnd,
     DataInputIteratorY yStart, DataInputIteratorY yEnd,
    const svgStyle & style)
  {
    svgStream << "<polyline points=\"";

    DataInputIteratorY itery = yStart;
    for (DataInputIteratorX iterx = xStart;
      iterx != xEnd;
      std::advance(iterx, 1), std::advance(itery, 1))
    {
      svgStream << *iterx << ',' << *itery << ' ';
    }
    svgStream << "\""
      << style.getSvgStream() + (style.bTooltip() ? "</polyline>" : "/>\n");
  }

  ///Close the svg tag.
  std::ostringstream & closeSvgFile()
  {
    svgStream <<"</svg>";
    return svgStream;
  }
private:
  std::ostringstream svgStream;
};

/// Helper to draw a SVG histogram
/// ____
/// |  |   ___ |
/// |  |__|  | |
/// |  |  |  | |
/// -----------|
struct svgHisto
{
  template<typename T>
  static std::string stringifier(const T & t)
  {
    std::ostringstream os;
    os << t;
    return os.str();
  }
  template<typename T>
  void draw(const std::vector<T> & vec_value,
    const std::pair<float, float> & range,
    const std::string & sFilename,
    const float W, const float H)
  {
    if (vec_value.empty())  {
      return;
    }
    //-- Max value
    const T maxi = *max_element(vec_value.begin(), vec_value.end());
    const size_t n = vec_value.size();

    const float scaleFactorY = H / static_cast<float>(maxi);
    const float scaleFactorX = W / static_cast<float>(n);

    svgDrawer svgStream;

    for (size_t i = 0; i < vec_value.size(); ++i)
    {
      const size_t dist = i;
      const T val = vec_value[i];
      std::ostringstream os;
      os << '(' << range.first + dist/float(n) * (range.second-range.first) << ',' << val << ')';
      svgStyle style = svgStyle().fill("blue").stroke("black", 1.0).tooltip(os.str());
      svgStream.drawRectangle(
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
    svgStyle styleAxis = svgStyle().stroke("black", 1.0f);
    // Draw X Axis
    svgStream.drawText(.05f*W, 1.2f*H, .1f*H, stringifier(range.first), "black");
    svgStream.drawText(   W, 1.2f*H, .1f*H, stringifier(range.second), "black");
    svgStream.drawLine(0, 1.1f*H, W, 1.1f*H, styleAxis);
    // Draw Y Axis
    svgStream.drawText(1.2f*W, .1f*H, .1f*H, stringifier(maxi), "black");
    svgStream.drawText(1.2f*W, H, .1f*H, "0", "black");
    svgStream.drawLine(1.1f*W, 0, 1.1f*W, H, styleAxis);

    std::ofstream svgFileStream( sFilename.c_str());
    svgFileStream << svgStream.closeSvgFile().str();
  }
};

} // namespace svg

#endif // THE_FRENCH_LEAF_SVG_DRAWER_H
