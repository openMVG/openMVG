// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graphics/color_gradient.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/image_io.hpp"
#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::graphics;

TEST(ColorGradient, BoundCheck) {
  const Color_Gradient heatMapGradient(Color_Gradient::k2BlueRedHeatMap());

  float r, g, b;

  heatMapGradient.getColor(0.f, r, g, b);
  EXPECT_EQ(0.f, r);
  EXPECT_EQ(0.f, g);
  EXPECT_EQ(1.f, b);

  heatMapGradient.getColor(1.f, r, g, b);
  EXPECT_EQ(1.f, r);
  EXPECT_EQ(0.f, g);
  EXPECT_EQ(0.f, b);

  heatMapGradient.getColor(.5f, r, g, b);
  EXPECT_EQ(.5f, r);
  EXPECT_EQ(0.f, g);
  EXPECT_EQ(.5f, b);
}

TEST(ColorGradient, GradientDemo) {

  const int w = 100, h = 10;
  Image<RGBColor> image(w, h);
  image.fill(BLACK);


  {
    const Color_Gradient heatMapGradient(Color_Gradient::k2BlueRedHeatMap());
    float r, g, b;

    for (int i = 0; i < w; ++i)
    {
      heatMapGradient.getColor(static_cast<float>(i) / w, r, g, b);
      image.col(i).fill(RGBColor(r * 255, g * 255, b * 255));
    }
    EXPECT_TRUE(WriteImage("color_gradient_2Colors.png", image));
  }

  {
    const Color_Gradient heatMapGradient;
    float r, g, b;

    for (int i = 0; i < w; ++i)
    {
      heatMapGradient.getColor(static_cast<float>(i) / w, r, g, b);
      image.col(i).fill(RGBColor(r * 255, g * 255, b * 255));
    }
    EXPECT_TRUE(WriteImage("color_gradient_5Colors.png", image));
  }

}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
