// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong,Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_ETH3D_HPP
#define IO_READ_GT_ETH3D_HPP

#include "io_readGTInterface.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/geometry/pose3.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

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
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the pose information
  std::map<int,Cameras_Data> camera_datas;
public:
  bool loadGT()
  {   
    // Check all the files under the path 
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );
    std::sort(gt_files.begin(), gt_files.end());

    // Because the ETH3D's Data store all the data in one file
    // So make sure there is only two file under the gt_dir
    if (gt_files.size()!=2)
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
      while (line_stream) 
      {
        temp_camera.parameter_.emplace_back();
        line_stream >> temp_camera.parameter_.back();
      }
      camera_datas.insert({temp_camera.id_,temp_camera});
    }
    camera_data_file.close();

    // Load the gt_data from the file
    std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_,"images.txt"), std::ifstream::in);
    if (!gt_file) 
    {
      return false;
    }
    int image_number_count = 0;
    while (gt_file) 
    {
      std::string line;
      std::getline(gt_file, line);
      if (line.size() == 0 || line[0] == '#') 
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
      if (line.size() == 0 || line[0] == '#') 
      {
        std::cout<<line<<std::endl;
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
           0,0,1;                        
      
      // Read feature observations line.
      // No use for us
      std::getline(gt_file, line);
      
      cameras::PinholeCamera camera_temp(K,R,t);        
      cameras_data_.push_back(camera_temp);

      // Parse image name
      images_.push_back( stlplus::filename_part(image_name) );
    }
    gt_file.close();

    return true;
  };

  bool loadImages()
  {
    this->sfm_data_.s_root_path = this->image_dir_; // Setup main image root_path

    Views & views = this->sfm_data_.views;
    Poses & poses = this->sfm_data_.poses;
    Intrinsics & intrinsics = this->sfm_data_.intrinsics;

    C_Progress_display my_progress_bar( images_.size(), std::cout, "\n- Image listing From Known Poses- \n" );
    std::ostringstream error_report_stream;
    std::vector<cameras::PinholeCamera>::const_iterator iter_camera = cameras_data_.begin();
    for ( std::vector<std::string>::const_iterator iter_image = images_.begin();
      iter_image != images_.end();
      ++iter_image,++iter_camera, ++my_progress_bar )
    {
      const std::string sImageFilename = stlplus::create_filespec( this->image_dir_, *iter_image );
      const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

      // Test if the image format is supported
      if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
      {
        error_report_stream << sImFilenamePart << ": Unkown image file format." << "\n";
        continue; // Image cannot be opened
      }

      if (sImFilenamePart.find("mask.png") != std::string::npos
      || sImFilenamePart.find("_mask.png") != std::string::npos)
      {
        error_report_stream << sImFilenamePart << " is a mask image" << "\n";
        continue;
      }

      // Test if this image can be read
      ImageHeader imgHeader;
      if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
        continue; // Image cannot be read

      Mat3 K = iter_camera->_K;
      double focal = (K(0,0)+K(1,1))/2.0; // Assume K(0,0)==K(1,1)
      double pxx   = K(0,2);
      double pyy   = K(1,2);

      std::shared_ptr<View> view                   = std::make_shared<View>(*iter_image, views.size(), views.size(), views.size(), imgHeader.width, imgHeader.height);
      Pose3 pose                                   = Pose3(iter_camera->_R,iter_camera->_C);
      std::shared_ptr<Pinhole_Intrinsic> intrinsic = std::make_shared<Pinhole_Intrinsic>(imgHeader.width,imgHeader.height,focal,pxx,pyy);

      // Add the view to the sfm_container
      views[view->id_view] = view;
      // Add the pose to the sfm_container
      poses[view->id_pose] = pose;
      // Add the intrinsic to the sfm_container
      intrinsics[view->id_intrinsic] = intrinsic;

    }

    // Display saved warning & error messages if any.
    if (!error_report_stream.str().empty())
    {
      std::cerr
        << "\nWarning & Error messages:" << std::endl
        << error_report_stream.str() << std::endl;
    }

    // Check group
    GroupSharedIntrinsics(this->sfm_data_);

    return true;
  }; 
};


#endif // IO_READ_GT_ETH3D_HPP
