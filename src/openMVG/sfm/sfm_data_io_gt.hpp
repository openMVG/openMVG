
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_IO_GT_HPP
#define OPENMVG_SFM_DATA_IO_GT_HPP

#include "sfm_data_io.hpp"

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <vector>

namespace openMVG {
namespace sfm {

static bool read_openMVG_Camera(const std::string & camName, cameras::Pinhole_Intrinsic & cam, geometry::Pose3 & pose)
{
  std::vector<double> val;
  if (stlplus::extension_part(camName) == "bin")
  {
    std::ifstream in(camName.c_str(), std::ios::in|std::ios::binary);
    if (!in.is_open())
    {
      std::cerr << "Error: failed to open file '" << camName << "' for reading" << std::endl;
      return false;
    }
    val.resize(12);
    in.read((char*)&val[0],(std::streamsize)12*sizeof(double));
    if (in.fail())
    {
      val.clear();
      return false;
    }
  }
  else
  {
    std::ifstream ifs;
    ifs.open( camName.c_str(), std::ifstream::in);
    if (!ifs.is_open()) {
      std::cerr << "Error: failed to open file '" << camName << "' for reading" << std::endl;
      return false;
    }
    while (ifs.good())
    {
      double valT;
      ifs >> valT;
      if (!ifs.fail())
        val.push_back(valT);
    }
  }

  //Update the camera from file value
  Mat34 P;
  if (stlplus::extension_part(camName) == "bin")
  {
    P << val[0], val[3], val[6], val[9],
      val[1], val[4], val[7], val[10],
      val[2], val[5], val[8], val[11];
  }
  else
  {
    P << val[0], val[1], val[2], val[3],
      val[4], val[5], val[6], val[7],
      val[8], val[9], val[10], val[11];
  }
  Mat3 K, R;
  Vec3 t;
  KRt_From_P(P, &K, &R, &t);
  cam = cameras::Pinhole_Intrinsic(0,0,K);
  // K.transpose() is applied to give [R t] to the constructor instead of P = K [R t]
  pose = geometry::Pose3(K.transpose() * P);
  return true;
}

static bool read_Strecha_Camera(const std::string & camName, cameras::Pinhole_Intrinsic & cam, geometry::Pose3 & pose)
{
  std::ifstream ifs;
  ifs.open( camName.c_str(), std::ifstream::in);
  if (!ifs.is_open()) {
    std::cerr << "Error: failed to open file '" << camName << "' for reading" << std::endl;
    return false;
  }
  std::vector<double> val;
  while (ifs.good() && !ifs.eof())
  {
    double valT;
    ifs >> valT;
    if (!ifs.fail())
      val.push_back(valT);
  }

  if (val.size() == 3*3 +3 +3*3 +3 + 3 || val.size() == 26) //Strecha cam
  {
    Mat3 K, R;
    K << val[0], val[1], val[2],
      val[3], val[4], val[5],
      val[6], val[7], val[8];
    R << val[12], val[13], val[14],
      val[15], val[16], val[17],
      val[18], val[19], val[20];

    Vec3 C (val[21], val[22], val[23]);
    // R need to be transposed
    cam = cameras::Pinhole_Intrinsic(0,0,K);
    cam.setWidth(val[24]);
    cam.setHeight(val[25]);
    pose = geometry::Pose3(R.transpose(), C);
  }
  else
  {
    return false;
  }
  return true;
}

/**
@brief Reads a set of Pinhole Cameras and its poses from a ground truth dataset.
@param[in] sGTPath, the directory where the camera files are located.
@param[out] sfm_data, the SfM_Data structure to put views/poses/intrinsics in.
**/
bool readGt(
  std::string sGTPath,
  SfM_Data & sfm_data
  )
{
  // IF GT_Folder exists, perform evaluation of the quality of rotation estimates
  if (!stlplus::is_folder(sGTPath)) {
    std::cout << std::endl << "There is not valid GT data to read from " << sGTPath << std::endl;
    return false;
  }

  // Switch between case to choose the file reader according to the file types in GT path
  bool (*fcnReadCamPtr)(const std::string &, cameras::Pinhole_Intrinsic &, geometry::Pose3&);
  std::string suffix;
  if (!stlplus::folder_wildcard(sGTPath, "*.bin", false, true).empty())
  {
    std::cout << "\nusing openMVG Camera";
    fcnReadCamPtr = &read_openMVG_Camera;
    suffix = "bin";
  }
  else if (!stlplus::folder_wildcard(sGTPath, "*.camera", false, true).empty())
  {
    std::cout << "\nusing Strechas Camera";
    fcnReadCamPtr = &read_Strecha_Camera;
    suffix = "camera";
  }
  else
  {
    throw std::logic_error(std::string("No camera found in ") + sGTPath);
  }

  std::cout << std::endl << "Read rotation and translation estimates" << std::endl;
  const std::string sImgPath = stlplus::folder_down(stlplus::folder_up(sGTPath), "images");
  std::map< std::string, Mat3 > map_R_gt;
  //Try to read .suffix camera (parse camera names)
  std::vector<std::string> vec_camfilenames =
    stlplus::folder_wildcard(sGTPath, "*."+suffix, false, true);
  std::sort(vec_camfilenames.begin(), vec_camfilenames.end());
  if (!vec_camfilenames.empty())
  {
    IndexT id = 0;
    for (std::vector<std::string>::const_iterator iter = vec_camfilenames.begin();
      iter != vec_camfilenames.end(); ++iter, ++id)
    {
      geometry::Pose3 pose;
      std::shared_ptr<cameras::Pinhole_Intrinsic> pinholeIntrinsic = std::make_shared<cameras::Pinhole_Intrinsic>();
      bool loaded = fcnReadCamPtr(stlplus::create_filespec(sGTPath, *iter), *pinholeIntrinsic.get(), pose);
      if (!loaded)
      {
        std::cout << "Failed to load: " << *iter << std::endl;
        return false;
      }

      const std::string sImgName = stlplus::basename_part(*iter);
      const std::string sImgFile = stlplus::create_filespec(sImgPath, sImgName);

      // Generate UID
      if (!stlplus::file_exists(sImgFile))
      {
        throw std::logic_error("Impossible to generate UID from this file, because it does not exists: "+sImgPath);
      }
      exif::Exif_IO_EasyExif exifReader;
      exifReader.open( sImgFile );
      const size_t uid = computeUID(exifReader, sImgName);

      // Update intrinsics with width and height of image
      sfm_data.views.emplace((IndexT)uid, std::make_shared<View>(stlplus::basename_part(*iter), (IndexT)uid, id, id, pinholeIntrinsic->w(), pinholeIntrinsic->h()));
      sfm_data.poses.emplace(id, pose);
      sfm_data.intrinsics.emplace(id, pinholeIntrinsic);
    }
  }
  return true;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_IO_GT_HPP