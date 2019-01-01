// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_KITTI_HPP
#define IO_READ_GT_KITTI_HPP

#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"

#include <fstream>
#include <iomanip>

// The feature of the Kitti's Data:
// 1. all the gt information are stored in two file:
//  - ID.txt for the camera poses,
//  - calib.txt for the camera intrinsic parameters.
// 2. the gt information's line number is the same with the corresponding image file name
class SfM_Data_GT_Loader_Kitti : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
public:
  bool loadGT() override
  {
    // Check all the files under the path
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );

    // Kitti's Data store all the data in two file
    // ID.txt    -> camera location ([R|C])
    // calib.txt -> camera calibration
    // So make sure there is only one file under the gt_dir
    if (gt_files.size() != 2 || gt_files.empty())
    {
      std::cerr << "Error: Maybe give wrong gt_dir!" << std::endl
        << "Make sure there you have only those 2 files under the gt_dir! (<ID>.txt and calib.txt)" << std::endl;
      return false;
    }

    // Read the camera calibration
    Mat3 calibration_matrix;
    auto calib_file_it = std::find(gt_files.begin(), gt_files.end(), "calib.txt");
    if (calib_file_it != gt_files.cend())
    {
      std::ifstream calib_file( stlplus::create_filespec(this->gt_dir_, "calib.txt"));
      if (!calib_file)
      {
        std::cerr << "Cannot open the calib.txt file" << std::endl;
        return false;
      }
      std::string temp;
      calib_file >> temp;
      Eigen::Matrix<double, 3, 4, Eigen::RowMajor> P;
      for (int i = 0; i < 12; ++i)
        calib_file >> *(P.data() + i);
      calibration_matrix = P.block<3, 3>(0, 0);
    }
    else
    {
      std::cerr << "Cannot find the expected calib.txt file" << std::endl;
      return false;
    }
    gt_files.erase(calib_file_it);

    // Load the gt_data from the file
    std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_, gt_files[0]));
    if (!gt_file)
    {
      std::cerr << "Error: Failed to open file '" << gt_files[0] << "' for reading" << std::endl;
      return false;
    }
    int image_number_count = 0;
    while (gt_file)
    {
      std::string line;
      std::getline(gt_file, line);
      ++image_number_count;
    }
    cameras_data_.reserve(image_number_count);

    gt_file.clear(std::ios::goodbit);
    gt_file.seekg(std::ios::beg);

    int frame_index = 0;
    while (gt_file)
    {
      Eigen::Matrix<double, 3, 4, Eigen::RowMajor> RC;
      for (int i = 0; i < 12; ++i)
      {
        gt_file >> *(RC.data() + i);
      }
      if (!gt_file)
        break;

      const Mat3 R = RC.block<3, 3>(0, 0);
      const Vec3 C = RC.block<3, 1>(0, 3);
      Mat34 P;
      P_From_KRt(calibration_matrix, R, - R * C, &P);
      cameras_data_.emplace_back(P);

      // Parse image name
      std::ostringstream os;
      os << std::setw(6) << std::setfill('0') << frame_index << ".png";
      images_.emplace_back( os.str() );

      frame_index++;
    }
    gt_file.close();

    return true;
  }

  bool loadImages() override
  {
    return LoadImages(this->image_dir_, this->images_, this->cameras_data_, this->sfm_data_);
  }
};


#endif // IO_READ_GT_KITTI_HPP
