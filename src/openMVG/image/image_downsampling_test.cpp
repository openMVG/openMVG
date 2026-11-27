// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_resampling.hpp"
#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::image;

TEST(Downsampling, AveragesTwoByTwoBlocks) {
  for (int width : {4, 5}) {
    Image<float> input(width, 3);
    for (int y = 0; y < input.Height(); ++y)
      for (int x = 0; x < input.Width(); ++x)
        input(y, x) = 10 * y + x;
    Image<float> output;
    ImageHalfSample(input, output);
    EXPECT_EQ(width / 2, output.Width());
    EXPECT_EQ(1, output.Height());
    EXPECT_NEAR(5.5f, output(0, 0), 1e-6);
    EXPECT_NEAR(7.5f, output(0, 1), 1e-6);
  }
}

TEST(Downsampling, AveragesByteAndColorPixels) {
  Image<unsigned char> gray(2, 2);
  gray << 0, 40, 80, 120;
  Image<unsigned char> gray_output;
  ImageHalfSample(gray, gray_output);
  EXPECT_EQ(60, int(gray_output(0, 0)));

  Image<RGBColor> color(2, 2);
  for (int i = 0; i < 4; ++i)
    color(i / 2, i % 2) = RGBColor(10 * i, 40 + 10 * i, 80 + 10 * i);
  Image<RGBColor> color_output;
  ImageHalfSample(color, color_output);
  EXPECT_EQ(15, int(color_output(0, 0).r()));
  EXPECT_EQ(55, int(color_output(0, 0).g()));
  EXPECT_EQ(95, int(color_output(0, 0).b()));
}

TEST(Downsampling, PreservesCheckerboardMeanAcrossLevels) {
  Image<float> input(4, 4);
  for (int y = 0; y < 4; ++y)
    for (int x = 0; x < 4; ++x)
      input(y, x) = (x + y) % 2;
  Image<float> half, quarter;
  ImageHalfSample(input, half);
  for (int y = 0; y < 2; ++y)
    for (int x = 0; x < 2; ++x)
      EXPECT_NEAR(0.5f, half(y, x), 1e-6);
  ImageHalfSample(half, quarter);
  EXPECT_NEAR(0.5f, quarter(0, 0), 1e-6);
}

TEST(Downsampling, HandlesImagesSmallerThanABlock) {
  for (const Vec2 & size : {Vec2(0, 0), Vec2(1, 1), Vec2(1, 4), Vec2(4, 1)}) {
    const Image<float> input(size(0), size(1));
    Image<float> output;
    ImageHalfSample(input, output);
    EXPECT_EQ(input.Width() / 2, output.Width());
    EXPECT_EQ(input.Height() / 2, output.Height());
    EXPECT_EQ(0, output.size());
  }
}

TEST(Downsampling, KeepsDecimationDistinct) {
  Image<float> input(4, 4);
  for (int y = 0; y < 4; ++y)
    for (int x = 0; x < 4; ++x)
      input(y, x) = 10 * y + x;
  Image<float> output;
  ImageDecimate(input, output);
  EXPECT_EQ(0.0f, output(0, 0));
  EXPECT_EQ(2.0f, output(0, 1));
  EXPECT_EQ(20.0f, output(1, 0));
  EXPECT_EQ(22.0f, output(1, 1));
}

int main() {
  TestResult result;
  return TestRegistry::runAllTests(result);
}
