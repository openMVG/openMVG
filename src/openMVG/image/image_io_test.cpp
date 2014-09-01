// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <cstdio>
#include <iostream>
#include <string>

#include "openMVG/image/image.hpp"
#include "testing/testing.h"

using namespace openMVG;
using std::string;

TEST(ReadJpg, Jpg_Color) {
  Image<RGBColor> image;
  string jpg_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_color.jpg";
  EXPECT_TRUE(ReadImage(jpg_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(3, image.Depth());
  EXPECT_EQ(image(0,0), RGBColor(255, 125, 11));
  EXPECT_EQ(image(0,1), RGBColor( 20, 127, 255));
}

TEST(ReadJpg, Jpg_Monochrome) {
  Image<unsigned char> image;
  string jpg_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_monochrome.jpg";
  EXPECT_TRUE(ReadImage(jpg_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}

TEST(ReadPng, Png_Color) {
  Image<RGBAColor> image;
  string png_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_color.png";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));
  // Depth is 4 (RGBA by default)
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(4, image.Depth());
  EXPECT_EQ(image(0,0), RGBAColor(255, 125, 10, 255));
  EXPECT_EQ(image(0,1), RGBAColor( 20, 127, 255,255));
}

TEST(ReadPng, Png_Monochrome) {
  Image<unsigned char> image;
  string png_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_monochrome.png";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}

TEST(GetFormat, filenames) {
  EXPECT_EQ(GetFormat("something.jpg"), openMVG::Jpg);
  EXPECT_EQ(GetFormat("something.png"), openMVG::Png);
  EXPECT_EQ(GetFormat("something.pnm"), openMVG::Pnm);
  EXPECT_EQ(GetFormat("something.tif"), openMVG::Tiff);
  EXPECT_EQ(GetFormat("/some/thing.JpG"), openMVG::Jpg);
  EXPECT_EQ(GetFormat("/some/thing.pNG"), openMVG::Png);
  EXPECT_EQ(GetFormat("some/thing.PNm"), openMVG::Pnm);
  EXPECT_EQ(GetFormat("some/thing.TIf"), openMVG::Tiff);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.JPG"), openMVG::Jpg);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.PNG"), openMVG::Png);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.PNM"), openMVG::Pnm);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.TIF"), openMVG::Tiff);
}

TEST(ImageIOTest, Png_Out) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string out_filename = ("test_write_png.png");
  EXPECT_TRUE(WriteImage(out_filename.c_str(), image));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}

TEST(ImageIOTest, Png_Out_Color) {
  Image<RGBColor> image(1,2);
  image(0,0) = RGBColor(255,127,0);
  image(1,0) = RGBColor(0,127,255);
  string out_filename = ("test_write_png_color.png");
  EXPECT_TRUE(WriteImage(out_filename.c_str(), image));

  Image<RGBColor> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}

TEST(ImageIOTest, InvalidFiles) {
  Image<unsigned char> image;
  string filename = string(THIS_SOURCE_DIR) + "/donotexist.jpg";
  EXPECT_FALSE(ReadImage(filename.c_str(), &image));
  EXPECT_FALSE(ReadImage("hopefully_unexisting_file", &image));
  remove(filename.c_str());
}

TEST(ImageIOTest, Jpg) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string filename = ("test_write_jpg.jpg");
  EXPECT_TRUE(WriteJpg(filename.c_str(), image, 100));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(filename.c_str());
}

TEST(ReadPnm, Pgm) {
  Image<unsigned char> image;
  string pgm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels.pgm";
  EXPECT_TRUE(ReadImage(pgm_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}

TEST(ReadPnm, PgmComments) {
  Image<unsigned char> image;
  string pgm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_gray.pgm";
  EXPECT_TRUE(ReadImage(pgm_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}


TEST(ImageIOTest, Pgm) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string out_filename = "test_write_pnm.pgm";
  EXPECT_TRUE(WriteImage(out_filename.c_str(),image));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}

TEST(ReadPnm, Ppm) {
  Image<RGBColor> image;
  string ppm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels.ppm";
  EXPECT_TRUE(ReadImage(ppm_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(3, image.Depth());
  EXPECT_EQ(image(0,0), RGBColor( (unsigned char)255));
  EXPECT_EQ(image(0,1), RGBColor( (unsigned char)0));
}

TEST(ImageIOTest, Ppm) {
  Image<RGBColor> image(1,2);
  image(0,0) = RGBColor((unsigned char)255);
  image(1,0) = RGBColor((unsigned char)0);
  string out_filename = "test_write_pnm.ppm";
  EXPECT_TRUE(WriteImage(out_filename.c_str(), image));

  Image<RGBColor> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}

TEST(ImageIOTest, Tiff_Gray) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string filename = ("test_write_tiff.tif");
  EXPECT_TRUE(WriteImage(filename.c_str(), image));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(filename.c_str());
}

TEST(ImageIOTest, Tiff_RGB) {
  Image<RGBColor> image(1,2);
  image(0,0) = RGBColor((unsigned char)255);
  image(1,0) = RGBColor((unsigned char)0);
  string filename = ("test_write_tiff.tif");
  EXPECT_TRUE(WriteImage(filename.c_str(), image));

  Image<RGBColor> read_image;
  EXPECT_TRUE(ReadImage(filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(filename.c_str());
}

TEST(ImageIOTest, Tiff_RGBA) {
  Image<RGBAColor> image(1,2);
  image(0,0) = RGBAColor(255, 125, 10, 255);
  image(1,0) = RGBAColor(2, 3, 4, 255);
  string filename = ("test_write_tiff.tif");
  EXPECT_TRUE(WriteImage(filename.c_str(), image));

  Image<RGBAColor> read_image;
  EXPECT_TRUE(ReadImage(filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(filename.c_str());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
