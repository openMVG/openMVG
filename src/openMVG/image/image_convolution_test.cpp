// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_convolution.hpp"
#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::image;

TEST(Convolution, ReflectsBothEndsOfRamp) {
  Image<float> input(5, 5);
  for (int y = 0; y < input.Height(); ++y)
    for (int x = 0; x < input.Width(); ++x)
      input(y, x) = x;
  Image<float> output;
  const Vec3 mean(0.25, 0.5, 0.25), identity(0, 1, 0);
  ImageSeparableConvolution(input, mean, identity, output);
  EXPECT_NEAR(0.5f, output(2, 0), 1e-6);
  EXPECT_NEAR(3.5f, output(2, 4), 1e-6);
  EXPECT_NEAR(2.0f, output(2, 2), 1e-6);
}

TEST(Convolution, ReflectedDerivativeVanishesAtBothEnds) {
  Image<float> input(5, 5);
  for (int y = 0; y < input.Height(); ++y)
    for (int x = 0; x < input.Width(); ++x)
      input(y, x) = x;
  Image<float> output;
  const Vec3 derivative(-0.5, 0, 0.5), identity(0, 1, 0);
  ImageSeparableConvolution(input, derivative, identity, output);
  EXPECT_NEAR(0.0f, output(2, 0), 1e-6);
  EXPECT_NEAR(0.0f, output(2, 4), 1e-6);
  EXPECT_NEAR(1.0f, output(2, 2), 1e-6);
}

TEST(Convolution, MatchesReflectedReference) {
  const auto reflect = [](int index, int size) {
    if (index < 0) return -index;
    if (index >= size) return 2 * size - 2 - index;
    return index;
  };
  Image<float> input(9, 7);
  for (int y = 0; y < input.Height(); ++y)
    for (int x = 0; x < input.Width(); ++x)
      input(y, x) = (x * x + 3 * y) % 17;
  for (int size : {3, 5}) {
    Vec kernel_x(size), kernel_y(size);
    for (int i = 0; i < size; ++i) {
      kernel_x(i) = double(i + 1) / (size * (size + 1) / 2);
      kernel_y(i) = double(size - i) / (size * (size + 1) / 2);
    }
    Image<float> output;
    ImageSeparableConvolution(input, kernel_x, kernel_y, output);
    for (int y = 0; y < input.Height(); ++y) {
      for (int x = 0; x < input.Width(); ++x) {
        double expected = 0;
        for (int ky = 0; ky < size; ++ky)
          for (int kx = 0; kx < size; ++kx)
            expected += kernel_y(ky) * kernel_x(kx) *
              input(reflect(y + ky - size / 2, input.Height()),
                    reflect(x + kx - size / 2, input.Width()));
        EXPECT_NEAR(expected, output(y, x), 2e-6);
      }
    }
  }
}

TEST(Convolution, PreservesConstantImage) {
  Image<float> input(5, 5, true, 7.0f), output;
  const Vec3 kernel(0.25, 0.5, 0.25);
  ImageSeparableConvolution(input, kernel, kernel, output);
  EXPECT_MATRIX_NEAR(input, output, 1e-6);
}

int main() {
  TestResult result;
  return TestRegistry::runAllTests(result);
}
