// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_BLENDED_MVS_HPP
#define IO_READ_GT_BLENDED_MVS_HPP

#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"

// BlendedMVS [1] Data importer:
// [1] BlendedMVS: A Large-scale Dataset for Generalized Multi-view Stereo Networks},
// Yao, Yao and Luo, Zixin and Li, Shiwei and Zhang, Jingyang and Ren,
//   Yufan and Zhou, Lei and Fang, Tian and Quan, Long,
// Computer Vision and Pattern Recognition (CVPR), 2020
//--------------
// Data storage:
//--------------
// ./cams
//  | id.txt
//  | ...
// ./blended_images
//  | id.jpg
//  | id_masked.jpg
//
// Note: id encoded on 8 digits
// Camera information:
//  extrinsic
//  [P] =   [R    |t]
//          [0 0 0 1]
//
// intrinsic
// [K]
//
// DEPTH_MIN DEPTH_INTERVAL (DEPTH_NUM DEPTH_MAX) // Ignored here


class SfM_Data_GT_Loader_BlendedMVS : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
public:
  bool loadGT() override
  {
    // Check if there is any camera files under the path
    std::vector<std::string> gt_files = stlplus::folder_wildcard(this->gt_dir_, "*_cam.txt");
    std::sort(gt_files.begin(), gt_files.end());

    if (gt_files.empty())
    {
      OPENMVG_LOG_ERROR << "The provided GT directory does not have any *_cam.txt file.";
      return false;
    }

    cameras_data_.reserve(gt_files.size());

    // Load the gt_data from the file
    for ( const auto & iter_gt_file : gt_files)
    {
      std::ifstream camera_data_file( stlplus::create_filespec(this->gt_dir_, iter_gt_file).c_str(), std::ifstream::in);
      if (!camera_data_file)
      {
        OPENMVG_LOG_ERROR << "Error: Failed to open file '" << iter_gt_file << "' for reading";
        continue;
      }
      // Read Projection and intrinsic matrices
      Mat3 R, K;
      Vec3 t;

      std::string line;
      std::getline(camera_data_file, line);
      if (line == "extrinsic")
      {
        std::vector<double> val;
        for (int i=0; i < 16; ++i)
        {
          double valT;
          camera_data_file >> valT;
          val.push_back(valT);
        }
        R << val[0], val[1], val[2],
             val[4], val[5], val[6],
             val[8], val[9], val[10];
        t << val[3], val[7], val[11];
      }

      std::getline(camera_data_file, line);
      std::getline(camera_data_file, line);
      std::getline(camera_data_file, line);
      if (line == "intrinsic")
      {
        std::vector<double> val;
        for (int i=0; i < 9; ++i)
        {
          double valT;
          camera_data_file >> valT;
          val.push_back(valT);
        }
        K << val[0], val[1], val[2],
             val[3], val[4], val[5],
             val[6], val[7], val[8];
      }
      Mat34 P;
      P_From_KRt(K, R, t, &P);
      cameras_data_.emplace_back(P);

      camera_data_file.close();

      const std::string image_name = stlplus::basename_part(iter_gt_file).substr(0, 8) + ".jpg";
      images_.emplace_back( image_name );
    }

    return true;
  }

  bool loadImages() override
  {
    return LoadImages(this->image_dir_, this->images_, this->cameras_data_, this->sfm_data_);
  }
};

#endif // IO_READ_GT_BLENDED_MVS_HPP
