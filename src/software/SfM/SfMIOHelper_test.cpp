// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/SfM/SfMIOHelper.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::SfMIO;

TEST(SfMIOHelper, EmptyFile) {
  std::ostringstream os;
  os.str("");

  const std::string sListsFile = stlplus::create_filespec("./","lists","txt");
  std::ofstream file(sListsFile.c_str());
  file << os.str();
  file.close();

  // Read data from the lists.txt file
  std::vector<CameraInfo> vec_camImageNames;
  std::vector<IntrinsicCameraInfo> vec_intrinsicGroups;

  EXPECT_FALSE(
    openMVG::SfMIO::loadImageList( vec_camImageNames, vec_intrinsicGroups, sListsFile));

  EXPECT_EQ(0, vec_intrinsicGroups.size());
}

TEST(SfMIOHelper, UniqueIntrinsicGroup) {

  std::ostringstream os;
  os   //ImaName;W;H;FocalPix;KMatrix
    << "0.jpg;2592;1936;2052.91;0;1278.59;0;2052.91;958.71;0;0;1";

  const std::string sListsFile = stlplus::create_filespec("./","lists","txt");
  std::ofstream file(sListsFile.c_str());
  file << os.str();
  file.close();

  // Read data from the lists.txt file
  std::vector<CameraInfo> vec_camImageNames;
  std::vector<IntrinsicCameraInfo> vec_intrinsicGroups;

  EXPECT_TRUE(
    openMVG::SfMIO::loadImageList( vec_camImageNames, vec_intrinsicGroups, sListsFile));

  EXPECT_EQ(1, vec_intrinsicGroups.size());
}

TEST(SfMIOHelper, SameCameraDifferentFocal) {

  std::ostringstream os;
  os   //ImaName;W;H;FocalPix;KMatrix
    << "DSC00402.JPG;4912;3264;3344;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA" <<'\n'
    << "DSC00403.JPG;4912;3264;6644;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA";

  const std::string sListsFile = stlplus::create_filespec("./","lists","txt");
  std::ofstream file(sListsFile.c_str());
  file << os.str();
  file.close();

  // Read data from the lists.txt file
  std::vector<CameraInfo> vec_camImageNames;
  std::vector<IntrinsicCameraInfo> vec_intrinsicGroups;

  EXPECT_TRUE(
    openMVG::SfMIO::loadImageList( vec_camImageNames, vec_intrinsicGroups, sListsFile));

  EXPECT_EQ(2, vec_intrinsicGroups.size());
}

TEST(SfMIOHelper, ManyCameraDifferentFocal) {

  std::ostringstream os;
  os   //ImaName;W;H;FocalPix;CamMaker;CamName
    << "DSC00402.JPG;4912;3264;3344.34;SONY;NEX-3N" <<'\n'
    << "100_7100.JPG;2832;2128;2881.25;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA" <<'\n'
    << "DSC00403.JPG;4912;3264;3344.34;SONY;NEX-3N" <<'\n' // same group as DSC00402
    << "100_7101.JPG;2832;2128;4881.25;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA" <<'\n'
    << "IMG_3266.JPG;5472;3648;4377.6;Canon;Canon EOS 70D" << '\n'
    << "100_7102.JPG;2832;2128;6881.25;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA" <<'\n'
    << "IMG_3267.JPG;5472;3648;6677.6;Canon;Canon EOS 70D" << '\n'
    << "IMG_3210.JPG;5616;3744;13260;Canon;Canon EOS 5D Mark II" << '\n'
    << "IMG_3211.JPG;5616;3744;10260;Canon;Canon EOS 5D Mark II" << '\n'
    << "IMG_3212.JPG;5616;3744;13260;Canon;Canon EOS 5D Mark II" << '\n' // same group as IMG_3210
    << "IMG_3212.JPG;5616;3744;Xylus;Junior"; // not known camera

  const std::string sListsFile = stlplus::create_filespec("./","lists","txt");
  std::ofstream file(sListsFile.c_str());
  file << os.str();
  file.close();

  // Read data from the lists.txt file
  std::vector<CameraInfo> vec_camImageNames;
  std::vector<IntrinsicCameraInfo> vec_intrinsicGroups;

  EXPECT_TRUE(
    openMVG::SfMIO::loadImageList( vec_camImageNames, vec_intrinsicGroups, sListsFile));

  EXPECT_EQ(9, vec_intrinsicGroups.size());
  // Check intrinsic group Ids correctness
  const size_t intrinsicGTIDs [] = {0,1,0,2,3,4,5,6,7,6,8};
  for (size_t i =0; i < 9; ++i) {
    EXPECT_EQ(intrinsicGTIDs[i], vec_camImageNames[i].m_intrinsicId);
  }
}

TEST(SfMIOHelper, KnowAndUnknowCamera) {

  std::ostringstream os;
  os   //ImaName;W;H;FocalPix;CamMaker;CamName
    << "DSC00402.JPG;4912;3264;3344.34;SONY;NEX-3N" <<'\n'
    << "0.jpg;4912;3264;3344.34;0;2456;0;3344.34;1632;0;0;1";

  const std::string sListsFile = stlplus::create_filespec("./","lists","txt");
  std::ofstream file(sListsFile.c_str());
  file << os.str();
  file.close();

  // Read data from the lists.txt file
  std::vector<CameraInfo> vec_camImageNames;
  std::vector<IntrinsicCameraInfo> vec_intrinsicGroups;

  EXPECT_TRUE(
    openMVG::SfMIO::loadImageList( vec_camImageNames, vec_intrinsicGroups, sListsFile));

  EXPECT_EQ(2, vec_intrinsicGroups.size());
}

TEST(SfMIOHelper, ThreeIntrinsicGroup_KMatrix) {

  std::ostringstream os;
  os   //ImaName;W;H;FocalPix;KMatrix
    << "0.jpg;2592;1936;2052.91;0;1278.59;0;2052.91;958.71;0;0;1" << '\n'
    << "1.jpg;2592;1936;2052.91;0;1278.59;0;2052.91;958.71;0;0;1" << '\n'
    << "2.jpg;2592;1936;2059.94;0;1274.91;0;2059.94;967.70;0;0;1" << '\n'
    << "3.jpg;2592;1936;2044.66;0;1253.00;0;2044.66;981.52;0;0;1" << '\n'
    << "4.jpg;2592;1936;2052.91;0;1278.59;0;2052.91;958.71;0;0;1" << '\n'
    << "5.jpg;2592;1936;2059.94;0;1274.91;0;2059.94;967.70;0;0;1" << '\n'
    << "6.jpg;2592;1936;2044.66;0;1253.00;0;2044.66;981.52;0;0;1";

  const std::string sListsFile = stlplus::create_filespec("./","lists","txt");
  std::ofstream file(sListsFile.c_str());
  file << os.str();
  file.close();

  // Read data from the lists.txt file
  std::vector<CameraInfo> vec_camImageNames;
  std::vector<IntrinsicCameraInfo> vec_intrinsicGroups;

  EXPECT_TRUE(
    openMVG::SfMIO::loadImageList( vec_camImageNames, vec_intrinsicGroups, sListsFile));

  EXPECT_EQ(3, vec_intrinsicGroups.size());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
