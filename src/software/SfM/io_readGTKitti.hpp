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
    // Look into the image dir to detect gray or rgb
    // Load the first image and detect the channel number
    std::cout<<stlplus::create_filespec(this->image_dir_,"image_0")<<std::endl;
    bool left_camera = stlplus::folder_exists(stlplus::create_filespec(this->image_dir_,"image_0"));
    bool left_camera_gray = true;

    ImageHeader imgHeader;
    if (openMVG::image::ReadImageHeader(stlplus::create_filespec( stlplus::create_filespec(this->image_dir_,"image_0"), "000000.png" ).c_str(), &imgHeader))
    {
      if(imgHeader.channels==3)
        left_camera_gray = false;
    }

    bool right_camera = stlplus::folder_exists(stlplus::create_filespec(this->image_dir_,"image_1"));
    bool right_camera_gray = true;
    if (openMVG::image::ReadImageHeader(stlplus::create_filespec( stlplus::create_filespec(this->image_dir_ , "image_1"), "000000.png" ).c_str(), &imgHeader))
    {
      if(imgHeader.channels==3)
        right_camera_gray = false;
    }

    // Check all the files under the path
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );
        
    // Because the Kitti's Data store all the data in one file
    // So make sure there is only two files under the gt_dir
    // One for calib.txt and one for pos information
    if (gt_files.size()!=2)
    {
      std::cerr << "Error: Maybe give wrong gt_dir!"<<std::endl
        <<"Make sure there in only two files under the gt_dir!"<<std::endl;
      return false;
    }

    // Load calib information
    std::vector<Mat34> camera_calib;
    std::ifstream calib_file(stlplus::create_filespec(this->gt_dir_,"calib.txt"),std::ifstream::in);
    if(!calib_file)
    {
      std::cerr<< "Error: Failed to open file '" << "calib.txt" << "' for reading" << std::endl;
      return false;
    }
    while(calib_file.good())
    {
      std::string name_temp;
      Mat34 calib_temp;
      calib_file>>name_temp
                >>calib_temp(0,0)>>calib_temp(0,1)>>calib_temp(0,2)>>calib_temp(0,3)
                >>calib_temp(1,0)>>calib_temp(1,1)>>calib_temp(1,2)>>calib_temp(1,3)
                >>calib_temp(2,0)>>calib_temp(2,1)>>calib_temp(2,2)>>calib_temp(2,3);
      camera_calib.emplace_back(calib_temp);
    }
    calib_file.close();


    std::vector<std::string>::iterator pos_file_iter = find( gt_files.begin( ), gt_files.end( ), "calib.txt" );
    pos_file_iter = pos_file_iter==gt_files.begin() ? gt_files.end() : gt_files.begin();
    // Load the gt_data from the file
    std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_,*pos_file_iter), std::ifstream::in);
    if (!gt_file) 
    {
      std::cerr << "Error: Failed to open file '" << *pos_file_iter << "' for reading" << std::endl;
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
    while (gt_file.good()) 
    {
      std::vector<double> val;
      for(int i=0;i<12;i++)
      {
        double valT;
        gt_file >> valT;
        val.push_back(valT);
      }
      
      Mat3 R;
      R << val[0], val[1], val[2],
           val[4], val[5], val[6],
           val[8], val[9], val[10];
      Vec3 c (val[3], val[7], val[11]);

      const Vec3 t (-R.transpose() * c);
      R.transposeInPlace();

      Eigen::Matrix<double, 4, 4> P;
      P(0,0) = R(0,0),P(0,1) = R(0,1),P(0,2) = R(0,2),P(0,3) = t(0,0);
      P(1,0) = R(1,0),P(1,1) = R(1,1),P(1,2) = R(1,2),P(1,3) = t(1,0);
      P(2,0) = R(2,0),P(2,1) = R(2,1),P(2,2) = R(2,2),P(2,3) = t(2,0);
      P(3,0) = 0,     P(3,1) = 0     ,P(3,2) = 0,     P(3,3) = 1     ;

      if(left_camera)
      {
        if(left_camera_gray)
        {
          Mat34 P_gray = camera_calib[0]*P;
          cameras::PinholeCamera camera_temp = cameras::PinholeCamera(P_gray);   
          cameras_data_.push_back(camera_temp);
        }
        else
        {
          Mat34 P_color = camera_calib[2]*P;
          cameras::PinholeCamera camera_temp = cameras::PinholeCamera(P_color);  
          cameras_data_.push_back(camera_temp);
        }
        
      }
      
      if(right_camera)
      {
        if(right_camera_gray)
        {
          Mat34 P_gray = camera_calib[1]*P;
          cameras::PinholeCamera camera_temp = cameras::PinholeCamera(P_gray);   
          cameras_data_.push_back(camera_temp);
        }
        else
        {
          Mat34 P_color = camera_calib[3]*P;
          cameras::PinholeCamera camera_temp = cameras::PinholeCamera(P_color);  
          cameras_data_.push_back(camera_temp);
        }
      }
      
      // Parse image name
      std::string image_name = "";
      std::string image_name_code = int2str(frame_index)+".png";
	    image_name.append(10-image_name_code.size(), '0');
	    image_name += image_name_code;
      if(left_camera)
        images_.push_back( "image_0/"+image_name );
      if(right_camera)
        images_.push_back( "image_1/"+image_name );

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
      std::cout<<sImageFilename<<std::endl;
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
