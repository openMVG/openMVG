// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong,Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_KITTI_HPP
#define IO_READ_GT_KITTI_HPP


#include "io_readGTInterface.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/geometry/pose3.hpp"

#include <sstream>

std::string int2str(const int int_temp)
{
    std::stringstream stream;
    stream<<int_temp;
    return stream.str();
}

// The feature of the Kitti's Data:
// 1. all the gt information are store in one file
// 2. the gt information's line number is the same with the corresponded image file name
class SfM_Data_GT_Loader_Kitti : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the pose information
public:
  bool loadGT()
  {   
    // Check all the files under the path
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );
        
    // Because the Kitti's Data store all the data in one file
    // So make sure there is only one file under the gt_dir
    if (gt_files.size()!=1)
    {
      std::cerr << "Error: Maybe give wrong gt_dir!"<<std::endl
        <<"Make sure there in only one file under the gt_dir!"<<std::endl;
      return false;
    }

    // Load the gt_data from the file
    std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_,gt_files[0]), std::ifstream::in);
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
      image_number_count++;
    }
    cameras_data_.reserve(image_number_count);

    gt_file.clear(std::ios::goodbit);
    gt_file.seekg(std::ios::beg);
    
    int frame_index = 0;
    while (gt_file) 
    {
      std::vector<double> val;
      for(int i=0;i<12;i++)
      {
        double valT;
        gt_file >> valT;
        val.push_back(valT);
      }

      Mat34 P;
      P << val[0], val[1], val[2], val[3],
           val[4], val[5], val[6], val[7],
           val[8], val[9], val[10], val[11];

      cameras::PinholeCamera camera_temp = cameras::PinholeCamera(P);        
      cameras_data_.push_back(camera_temp);

      // Parse image name
      std::string image_name = "";
      std::string image_name_code = int2str(frame_index)+".png";
	    image_name.append(10-image_name_code.size(), '0');
	    image_name += image_name_code;
      images_.push_back( stlplus::filename_part(image_name) );

      frame_index++;
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
    
    C_Progress_display my_progress_bar( images_.size(),std::cout, "\n- Image listing From Known Poses- \n" );
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
        error_report_stream
          << sImFilenamePart << ": Unkown image file format." << "\n";
        continue; // Image cannot be opened
      }

      if (sImFilenamePart.find("mask.png") != std::string::npos
        || sImFilenamePart.find("_mask.png") != std::string::npos)
      {
        error_report_stream
          << sImFilenamePart << " is a mask image" << "\n";
        continue;
      }

      // Test if this Image can be read
      ImageHeader imgHeader;
      if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
        continue; // Image cannot be read

      Mat3 K = iter_camera->_K;
      double focal = (K(0,0)+K(1,1))/2.0; //Assume K(0,0)==K(1,1)
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


#endif // IO_READ_GT_KITTI_HPP
