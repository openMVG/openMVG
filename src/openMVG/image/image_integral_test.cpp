// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image_integral.hpp"

#include "testing/testing.h"

#include <random>

using namespace openMVG;
using namespace openMVG::image;

TEST(integral_image, empty) {
  Image<uint8_t> image;

  Image<uint8_t> integral_image;
  IntegralImage(image, &integral_image);
  EXPECT_EQ(0, integral_image.Width());
  EXPECT_EQ(0, integral_image.Height());
}

TEST(integral_image, sum) {

  using Matuchar = Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>;

  // Assert an image with 0 produce a zero integral_image
  {
    const Matuchar mat = Matuchar::Zero(3, 7);
    const Image<uint8_t> image_zeros(mat);
    Image<uint8_t>  integral_image;
    IntegralImage(image_zeros, &integral_image);
    EXPECT_EQ(0, integral_image(2, 6));
  }

  // Define an image from a matrix with some non zero value
  Matuchar mat(3, 7);
  mat.row(0) = Eigen::ArrayXf::LinSpaced(7, 0, 7).cast<uint8_t>();
  mat.row(1) = 2 + Eigen::ArrayXf::LinSpaced(7, 0, 7).cast<uint8_t>();
  mat.row(2) = 4 + Eigen::ArrayXf::LinSpaced(7, 0, 7).cast<uint8_t>();

  const Image<uint8_t> image(mat);
  // Since the image is tiny, we don't risk any overflow
  {
    Image<uint8_t> integral_image;
    IntegralImage(image, &integral_image);
    EXPECT_EQ(mat.sum(), integral_image(2, 6));
    EXPECT_TRUE( 0 != integral_image(2, 6));
  }
  // Try with another templated matrix output
  {
    Image<int> integral_image;
    IntegralImage(image, &integral_image);
    EXPECT_EQ(mat.sum(), integral_image(2, 6));
    EXPECT_TRUE( 0 != integral_image(2, 6));
  }
  // Try with another templated matrix output
  {
    Image<uint32_t> integral_image;
    IntegralImage(image, &integral_image);
    EXPECT_EQ(mat.sum(), integral_image(2, 6));
    EXPECT_TRUE( 0 != integral_image(2, 6));
  }
}

// Assert that BlockSum is working for every kernel size
TEST(integral_image, BlockSum) {

  // Define an image from a matrix with some non zero value
  std::random_device rd;
  std::mt19937 gen(std::mt19937::default_seed);

  // Feed the image with some random pixel values
  Image<int> image(100, 200);
  std::uniform_int_distribution<> pix(0, 127);
  for (int i = 0; i < 100; ++i)
  for (int j = 0; j < 200; ++j)
    image(j, i) = pix(gen);

  // Note:
  // we are using a int matrix since else .block(i,j,p,q).sum() will overflow

  Image<uint32_t> integral_image;
  IntegralImage(image, &integral_image);

  static const int max_kernel_size = 10;

  for (int i = 0; i < 60; ++i)
  {
    // Kernel size generator
    std::uniform_int_distribution<> kernel_width(1, max_kernel_size);
    const int block_size = kernel_width(gen);

    // Compute the kernel position boundary for a random point position
    std::uniform_int_distribution<> dis_x(max_kernel_size + 1, 100 - max_kernel_size - 1);
    std::uniform_int_distribution<> dis_y(max_kernel_size + 1, 200 - max_kernel_size - 1);
    const Eigen::Vector2i left_corner = {dis_x(gen), dis_y(gen)};
    const Eigen::Vector2i right_corner = left_corner + Eigen::Vector2i(block_size, block_size);

    EXPECT_EQ(image.GetMat().block(left_corner.y(),
                                   left_corner.x(),
                                   block_size + 1,
                                   block_size + 1).sum(),
              BlockSum(left_corner, right_corner, integral_image));
  }
}

TEST(integral_image, BorderAndSinglePixelRectangles) {
  Image<int> image(5, 4);
  for (int y = 0; y < image.Height(); ++y)
    for (int x = 0; x < image.Width(); ++x)
      image(y, x) = 10 * y + x + 1;
  Image<uint32_t> integral;
  IntegralImage(image, &integral);
  for (int top = 0; top < image.Height(); ++top)
    for (int left = 0; left < image.Width(); ++left)
      for (int bottom = top; bottom < image.Height(); ++bottom)
        for (int right = left; right < image.Width(); ++right) {
          int expected = 0;
          for (int y = top; y <= bottom; ++y)
            for (int x = left; x <= right; ++x)
              expected += image(y, x);
          EXPECT_EQ(expected, BlockSum(Eigen::Vector2i(left, top),
                                      Eigen::Vector2i(right, bottom), integral));
        }
}

TEST(integral_image, AvoidsOverflowWhenStoredSumsFit) {
  const Image<int> image(5, 5, true, 80000000);
  Image<int> integral;
  IntegralImage(image, &integral);
  EXPECT_EQ(2000000000, integral(4, 4));
  EXPECT_EQ(320000000, BlockSum(Eigen::Vector2i(3, 3),
                              Eigen::Vector2i(4, 4), integral));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
