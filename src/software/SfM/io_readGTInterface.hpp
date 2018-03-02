// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_INTERFACE_HPP
#define IO_READ_GT_INTERFACE_HPP

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/image/image_io.hpp"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <vector>
#include <string>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::image;
using namespace openMVG::sfm;


class SfM_Data_GT_Loader_Interface
{
protected:
  SfM_Data sfm_data_;
  std::string gt_dir_; // Path store the ground truth 
  std::string image_dir_; // Path store the images
  std::vector<std::string> images_; // Store image lists from the ground truth 

  virtual bool loadGT() = 0; // Load ground truth from gt_dir and get all the images name need to load
  virtual bool loadImages() = 0; // Load images from image_dir 

public:
  bool run(const std::string &gt_dir,const std::string &image_dir)
  {
  
    this->gt_dir_ = gt_dir;
    this->image_dir_ = image_dir;

    if (!loadGT())
    {
      std::cerr<<"Failed to Load Ground Truth!"<<std::endl;
      return false;
    };
    if (!loadImages())
    {
      std::cerr<<"Failed to Load Images!"<<std::endl;
      return false;
    };

    return true;
  }

  SfM_Data GetSfMData()
  {
    return this->sfm_data_;
  };

  int GetImageNumber()
  {
    return this->images_.size();
  };
};

#endif // IO_READ_GT_INTERFACE_HPP
