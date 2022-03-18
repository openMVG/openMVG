// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_ETH3D_HPP
#define IO_READ_GT_ETH3D_HPP

#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"

#include <algorithm>

struct Cameras_Data
{
    int id_;
    std::string model_name_;
    int width_;
    int height_;
    std::vector<double> parameter_;
};


// The feature of the ETH3D's Data:
// 1. all the gt information are stored in one file
// 2. the cameras' intrinsic information are stored in camera.txt
class SfM_Data_GT_Loader_ETH_3D : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
  std::map<int, Cameras_Data> camera_datas; // Store all the dataset camera data
public:
  bool loadGT() override
  {
    // Check all the files under the path
    const std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );

    // Because the ETH3D's Data store all the data in many files
    // Make sure there we have the desired file on disk
    if ( !std::count(gt_files.cbegin(), gt_files.cend(), std::string("cameras.txt"))
        || !std::count(gt_files.cbegin(), gt_files.cend(), std::string("images.txt")))
    {
      std::cerr << "Error: Maybe give wrong gt_dir!"<<std::endl
        <<"Make sure there in only two files(images.txt cameras.txt) under the gt_dir!"<<std::endl;
      return false;
    }

    // Read the camera file
    // Fix name "cameras.txt"
    // Reference:https://github.com/ETH3D/format-loader
    std::ifstream camera_data_file( stlplus::create_filespec(this->gt_dir_,"cameras.txt"), std::ifstream::in);
    if (!camera_data_file)
    {
      std::cerr << "Error: Failed to open file '" << stlplus::create_filespec(this->gt_dir_,"cameras.txt") << "' for reading" << std::endl;
      return false;
    }
    while (camera_data_file)
    {
      std::string line;
      std::getline(camera_data_file, line);
      if (line.size() == 0 || line[0] == '#')
      {
        continue;
      }

      Cameras_Data temp_camera;
      std::istringstream line_stream(line);
      line_stream >> temp_camera.id_ >> temp_camera.model_name_
                  >> temp_camera.width_ >> temp_camera.height_;
      if (temp_camera.model_name_ != "PINHOLE")
      {
        std::cerr << "Only the pinhole camera model is supported.\n"
         << "Please consider use the undistorted version of the dataset.\n"
         << "Or contribute to add the: " << temp_camera.model_name_ << " camera model to OpenMVG."<< std::endl;
        return false;
      }
      while (line_stream)
      {
        temp_camera.parameter_.emplace_back();
        line_stream >> temp_camera.parameter_.back();
      }
      camera_datas.insert({temp_camera.id_,temp_camera});
    }
    camera_data_file.close();

    // Load the gt_data from the file
    std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_, "images.txt"), std::ifstream::in);
    if (!gt_file)
    {
      return false;
    }
    int image_number_count = 0;
    while (gt_file)
    {
      std::string line;
      std::getline(gt_file, line);
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      image_number_count++;
      std::getline(gt_file, line);
    }
    cameras_data_.reserve(image_number_count);

    gt_file.clear(std::ios::goodbit);
    gt_file.seekg(std::ios::beg);

    while (gt_file)
    {
      std::string line;
      std::getline(gt_file, line);
      if (line.empty() || line[0] == '#')
      {
        continue;
      }

      // Read image info line.
      Mat3 K, R;
      Vec3 t;
      Eigen::Quaterniond quaternionf_rotation;
      std::istringstream image_stream(line);
      int image_id,camera_id;
      std::string image_name;
      image_stream >> image_id
                   >> quaternionf_rotation.w()
                   >> quaternionf_rotation.x()
                   >> quaternionf_rotation.y()
                   >> quaternionf_rotation.z()
                   >> t(0,0)
                   >> t(1,0)
                   >> t(2,0)
                   >> camera_id
                   >> image_name;

      R = quaternionf_rotation.toRotationMatrix();

      Cameras_Data temp_camera = camera_datas[camera_id];
      double focus = (temp_camera.parameter_[0]+temp_camera.parameter_[1])/2;
      K << focus, 0, temp_camera.parameter_[2],
           0, focus, temp_camera.parameter_[3],
           0, 0, 1;

      // Read feature observations line.
      // No use for us
      std::getline(gt_file, line);

      Mat34 P;
      P_From_KRt(K, R, t, &P);
      cameras_data_.emplace_back(P);

      // Parse image name
      images_.emplace_back( stlplus::filename_part(image_name) );
    }
    gt_file.close();

    return true;
  }

  bool loadImages() override
  {
    return LoadImages(this->image_dir_, this->images_, this->cameras_data_, this->sfm_data_);
  }
};


#endif // IO_READ_GT_ETH3D_HPP
