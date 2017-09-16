// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013-2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "exif_IO_EasyExif.hpp"

#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <memory>

using namespace std;
using namespace openMVG;
using namespace openMVG::exif;

const std::string sImg =
  stlplus::folder_part(
  stlplus::folder_part(
  stlplus::folder_up(std::string(THIS_SOURCE_DIR))))
    + "/openMVG_Samples/imageData/Exif_Test/100_7100.JPG";

TEST(Matching, Exif_IO_easyexif_ReadData_invalidFile)
{
  std::unique_ptr<Exif_IO> exif_io ( new Exif_IO_EasyExif( "tmp.jpg" ) );

  EXPECT_FALSE( exif_io->doesHaveExifInfo());
}

TEST(Matching, Exif_IO_easyexif_ReadData)
{
  std::unique_ptr<Exif_IO> exif_io ( new Exif_IO_EasyExif( sImg ) );

  EXPECT_TRUE( exif_io->doesHaveExifInfo());

  EXPECT_EQ( "EASTMAN KODAK COMPANY", exif_io->getBrand());
  EXPECT_EQ( "KODAK Z612 ZOOM DIGITAL CAMERA", exif_io->getModel());

  EXPECT_EQ( 2832, exif_io->getWidth());
  EXPECT_EQ( 2128, exif_io->getHeight());
  EXPECT_NEAR( 5.85, exif_io->getFocal(), 1e-2);

  EXPECT_EQ( "", exif_io->getLensModel());

  double val;
  EXPECT_FALSE(exif_io->GPSLatitude(&val));
  EXPECT_FALSE(exif_io->GPSLongitude(&val));
  EXPECT_FALSE(exif_io->GPSAltitude(&val));
}

TEST(Matching, Exif_IO_easyexif_Read_GPS_Data)
{
  const std::string sImg_gps = std::string(THIS_SOURCE_DIR) + "/image_data/gps_tag.jpg";
  std::unique_ptr<Exif_IO> exif_io ( new Exif_IO_EasyExif( sImg_gps ) );
  double val;

  EXPECT_TRUE(exif_io->GPSLatitude(&val));
  EXPECT_NEAR(47.5129, val, 1e-4);

  EXPECT_TRUE(exif_io->GPSLongitude(&val));
  EXPECT_NEAR(2.1513, val, 1e-4);

  EXPECT_TRUE(exif_io->GPSAltitude(&val));
  EXPECT_NEAR(120, val, 1e-0);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

