// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_STRECHA_HPP
#define IO_READ_GT_STRECHA_HPP


#include "io_readGTInterface.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/geometry/pose3.hpp"


// The feature of the Strecha's Data:
// 1. each gt file only stores one view's pose information
// 2. the gt file's name is image file name plus ".camera"
class SfM_Data_GT_Loader_Strecha : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the pose information
public:
  bool loadGT()
  {   
    // Check all the files under the path 
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );
    std::sort(gt_files.begin(), gt_files.end());

    cameras_data_.reserve(gt_files.size());

    // Read all the gt files and store them in the system
    for ( std::vector<std::string>::const_iterator iter_gt_file = gt_files.begin();
      iter_gt_file != gt_files.end();
      ++iter_gt_file)
    {
      std::ifstream ifs;
      ifs.open( stlplus::create_filespec(this->gt_dir_,(*iter_gt_file)).c_str(), std::ifstream::in);
      if (!ifs.is_open()) 
      {
        std::cerr << "Failed to open file '" << *iter_gt_file << "' for reading" << std::endl;
        continue;
      }
      std::vector<double> val;
      while (ifs.good() && !ifs.eof())
      {
        double valT;
        ifs >> valT;
        if (!ifs.fail())
        val.push_back(valT);
      }
      ifs.close();

      if (val.size() == 3*3 +3 +3*3 +3 + 3 || val.size() == 26) // Strecha cam
      {
        Mat3 K, R;
        K << val[0], val[1], val[2],
        val[3], val[4], val[5],
        val[6], val[7], val[8];
        R << val[12], val[13], val[14],
        val[15], val[16], val[17],
        val[18], val[19], val[20];

        Vec3 C (val[21], val[22], val[23]);
        // Strecha model is P = K[R^T|-R^T t];
        // My model is P = K[R|t], t = - RC
        const Vec3 t (-R.transpose() * C);
        R.transposeInPlace();

        PinholeCamera cam = cameras::PinholeCamera(K, R, t);
        cameras_data_.push_back(cam);

        images_.push_back( stlplus::basename_part(stlplus::create_filespec(this->gt_dir_,(*iter_gt_file))) ); // Parse image name
      }
      else
      {
        continue;
      }
    }
  
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


#endif // IO_READ_GT_STRECHA_HPP
