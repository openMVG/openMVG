// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_MIDDLEBURY_HPP
#define IO_READ_GT_MIDDLEBURY_HPP

#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"

// The feature of the MiddleBury's Data:
// 1. all the gt information are store in one file
// 2. the image file's name is stored in the gt file
class SfM_Data_GT_Loader_MiddleBury : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
public:
  bool loadGT() override
  {
    // Check all the files under the path
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );

    // Because the MiddleBury's Data store all the data in one file
    // So make sure there is only one file under the gt_dir
    if (gt_files.size()!=1)
    {
      std::cerr << "Error: Maybe give wrong gt_dir!"<<std::endl
        <<"Make sure there in only one file under the gt_dir!"<<std::endl;
      return false;
    }

    // Load the gt_data from the file
    std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_,gt_files[0]).c_str(), std::ifstream::in);
    if (!gt_file)
    {
      std::cerr << "Error: Failed to open file '" << gt_files[0] << "' for reading" << std::endl;
      return false;
    }

    // Image_number
    int image_count = 0;
    gt_file >> image_count;

    cameras_data_.reserve(image_count);

    for(int i=0;i<image_count;i++)
    {
      std::string image_name;
      Mat3 K, R;
      Vec3 t;

      gt_file >> image_name;
      gt_file >> K(0,0) >> K(0,1) >> K(0,2)
              >> K(1,0) >> K(1,1) >> K(1,2)
              >> K(2,0) >> K(2,1) >> K(2,2);
      gt_file >> R(0,0) >> R(0,1) >> R(0,2)
              >> R(1,0) >> R(1,1) >> R(1,2)
              >> R(2,0) >> R(2,1) >> R(2,2);
      gt_file >> t(0,0) >> t(1,0) >> t(2,0);

      Mat34 P;
      P_From_KRt(K, R, t, &P);
      cameras_data_.emplace_back(P);

      // Parse image name
      images_.push_back( image_name );

    }
    gt_file.close();

    return true;
  }

  bool loadImages() override
  {
    return LoadImages(this->image_dir_, this->images_, this->cameras_data_, this->sfm_data_);
  }
};

#endif // IO_READ_GT_MIDDLEBURY_HPP
