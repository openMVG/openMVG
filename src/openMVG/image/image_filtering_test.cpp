// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_filtering.hpp"

#include "testing/testing.h"

#include <iostream>

using namespace openMVG;
using namespace openMVG::image;
using namespace std;

TEST(Image, Convolution)
{
  Image<unsigned char> in(250,250);
  for (int i = 10; i < 250-10; i++)
    for (int j = 10; j < 250-10; j++)
    {
      in(j,i) = rand()%127+127;
    }

  EXPECT_TRUE(in(5,5) == 0);
  EXPECT_TRUE(in(249-5,249-5) == 0);

  // filter
  Image<unsigned char> outFiltered(250,250);
  ImageGaussianFilter( in, 6.0, outFiltered);

  EXPECT_TRUE(WriteImage("in.png", in));
  EXPECT_TRUE(WriteImage("outfilter.png", outFiltered));

  // Check that gaussian filtering have smooth at border of the white random square
  EXPECT_TRUE(outFiltered(5,5)>0);
  EXPECT_TRUE(outFiltered(250-5,250-5)>0);
}

TEST(Image, Convolution_MeanBoxFilter_Separable)
{
  Image<unsigned char> in(40,40);
  in.block(10,10,20,20).fill(255.f);
  Vec3 meanBoxFilterKernel(1.f/3.f, 1.f/3.f, 1.f/3.f);
  Image<unsigned char> out;
  ImageSeparableConvolution( in , meanBoxFilterKernel, meanBoxFilterKernel , out);
}

TEST(Image, Convolution_MeanBoxFilter)
{
  Image<unsigned char> in(40,40);
  in.block(10,10,20,20).fill(255.f);
  Mat3 meanBoxFilterKernel;
  meanBoxFilterKernel.fill(1.f/9.f);
  Image<unsigned char> out;
  ImageConvolution(in, meanBoxFilterKernel, out);
}

TEST(Image, Convolution_Scharr_X_Y)
{
  Image<float> in(40,40);
  in.block(10,10,20,20).fill(255.f);

  Image<float> outFiltered(40,40);

  ImageScaledScharrXDerivative( in, outFiltered, 1);

  // X dir
  EXPECT_EQ(127.5f, outFiltered(20,10));
  EXPECT_EQ(-127.5f, outFiltered(20,30));
  // Y dir
  EXPECT_EQ(0.f, outFiltered(10,20));
  EXPECT_EQ(0.f, outFiltered(30,20));
  // Check it exist a vertical black band
  EXPECT_EQ(0.f, outFiltered.block(0,10+3,40,20-2*3).array().abs().sum());

  EXPECT_TRUE(WriteImage("in_Scharr.png", Image<unsigned char>(in.cast<unsigned char>())));
  EXPECT_TRUE(WriteImage("out_ScharrX.png", Image<unsigned char>(outFiltered.cast<unsigned char>())));

  outFiltered.fill(0.0f);
  ImageScaledScharrYDerivative( in, outFiltered, 1);

  // X dir
  EXPECT_EQ(0.f, outFiltered(20,10));
  EXPECT_EQ(0.f, outFiltered(20,30));
  // Y dir
  EXPECT_EQ(127.5f, outFiltered(10,20));
  EXPECT_EQ(-127.5f, outFiltered(30,20));
  // Check it exist a horizontal black band
  EXPECT_EQ(0.f, outFiltered.block(10+3,0,20-2*3,40).array().abs().sum());
  EXPECT_TRUE(WriteImage("out_ScharrY.png", Image<unsigned char>(outFiltered.cast<unsigned char>())));
}

TEST(Image, Convolution_Sobel_X_Y)
{
  Image<float> in(40,40);
  in.block(10,10,20,20).fill(255.f);

  Image<float> outFiltered(40,40);

  ImageSobelXDerivative( in, outFiltered);

  // X dir
  EXPECT_EQ(127.5f, outFiltered(20,10));
  EXPECT_EQ(-127.5f, outFiltered(20,30));
  // Y dir
  EXPECT_EQ(0.f, outFiltered(10,20));
  EXPECT_EQ(0.f, outFiltered(30,20));
  // Check it exist a vertical black band
  EXPECT_EQ(0.f, outFiltered.block(0,10+3,40,20-2*3).array().abs().sum());

  EXPECT_TRUE(WriteImage("in_Scharr.png", Image<unsigned char>(in.cast<unsigned char>())));
  EXPECT_TRUE(WriteImage("out_SobelX.png", Image<unsigned char>(outFiltered.cast<unsigned char>())));

  outFiltered.fill(0.0f);
  ImageSobelYDerivative( in, outFiltered);

  // X dir
  EXPECT_EQ(0.f, outFiltered(20,10));
  EXPECT_EQ(0.f, outFiltered(20,30));
  // Y dir
  EXPECT_EQ(127.5f, outFiltered(10,20));
  EXPECT_EQ(-127.5f, outFiltered(30,20));
  // Check it exist a horizontal black band
  EXPECT_EQ(0.f, outFiltered.block(10+3,0,20-2*3,40).array().abs().sum());

  EXPECT_TRUE(WriteImage("out_SobelY.png", Image<unsigned char>(outFiltered.cast<unsigned char>())));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
