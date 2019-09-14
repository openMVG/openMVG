// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPHICS_COLOR_GRADIENT_HPP
#define OPENMVG_GRAPHICS_COLOR_GRADIENT_HPP

#include <algorithm>
#include <vector>

namespace openMVG
{
namespace graphics
{


// Port made from http://www.andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients
// - Use c++11 style
// - Add multiple possible color gradient initialization
//
// --------------------
// How to use the code:
// --------------------
// 1. Create the object and configure the color gradient
// --------------------
// Color_Gradient heatMapGradient(Color_Gradient::k2BlueRedHeatMap());
// // Re-Initialize to use a gradient based on 5 colors
// //heatMapGradient = Color_Gradient(Color_Gradient::k5ColorHeatMap());
// --------------------
// 2. Get the color corresponding to the % range you want
// --------------------
// float r,g,b;
// const float interpolated_ratio = .5f;
// heatMapGradient.getColor(interpolated_ratio, r, g, b);
// // [r,g,b] are in the range [0, 1]
//
class Color_Gradient
{
public:
  struct Color_Point  // Internal class used to store colors at different points in the gradient.
  {
    float r, g, b;   // Red, green and blue values of our color.
    float val;       // Position of our color along the gradient (between 0 and 1).
    Color_Point(float red, float green, float blue, float value)
      : r(red), g(green), b(blue), val(value) {}
  };

  using Color_Points = std::vector<Color_Point>;

  //-- Places a 5 color heatmap gradient into the "color" vector:
  static const Color_Points k5ColorHeatMap()
  {
    return
    {
      {0, 0, 1,   0.0f},      // Blue
      {0, 1, 1,   0.25f},     // Cyan
      {0, 1, 0,   0.5f},      // Green
      {1, 1, 0,   0.75f},     // Yellow
      {1, 0, 0,   1.0f}       // Red
    };
  }

  //-- Places a 2 color heatmap gradient (Blue to Red):
  static const Color_Points k2BlueRedHeatMap()
  {
    return
    {
      {0, 0, 1,   0.0f}, // Blue
      {1, 0, 0,   1.0f}  // Red
    };
  }

private:
  Color_Points color_map;      // An array of color points in ascending value

public:
  //-- Default constructor:
  explicit Color_Gradient(const Color_Points & rhs_color = k5ColorHeatMap())
    : color_map(rhs_color)  {}


  //-- Inputs a (value) between 0 and 1 and outputs the (red), (green) and (blue)
  //-- values representing that position in the gradient.
  void getColor
  (
    const float value,
    float & red,
    float & green,
    float & blue
  ) const
  {
    if (color_map.empty())
      return;

    for (int i = 0; i < color_map.size(); ++i)
    {
      const Color_Point &currC = color_map[i];
      if (value < currC.val)
      {
        const Color_Point &prevC  = color_map[ std::max(0, i-1) ];
        const float valueDiff    = (prevC.val - currC.val);
        const float fractBetween = (valueDiff==0) ? 0 : (value - currC.val) / valueDiff;
        red   = (prevC.r - currC.r) * fractBetween + currC.r;
        green = (prevC.g - currC.g) * fractBetween + currC.g;
        blue  = (prevC.b - currC.b) * fractBetween + currC.b;
        return;
      }
    }
    red   = color_map.back().r;
    green = color_map.back().g;
    blue  = color_map.back().b;
    return;
  }
};

} // namespace graphics
} // namespace openMVG

#endif // OPENMVG_GRAPHICS_COLOR_GRADIENT_HPP
