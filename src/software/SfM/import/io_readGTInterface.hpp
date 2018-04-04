// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_INTERFACE_HPP
#define IO_READ_GT_INTERFACE_HPP

#include "openMVG/sfm/sfm_data.hpp"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::sfm;

class SfM_Data_GT_Loader_Interface
{
protected:
  SfM_Data sfm_data_;
  std::string gt_dir_; // Path to the GT (ground truth) files
  std::string image_dir_; // Path store to the images
  std::vector<std::string> images_; // Store GT image names list

  virtual bool loadGT() = 0; // Load GT from gt_dir and get all the image names
  virtual bool loadImages() = 0; // Load images from image_dir

public:
  bool run(const std::string &gt_dir,const std::string &image_dir)
  {
    this->gt_dir_ = gt_dir;
    this->image_dir_ = image_dir;

    if (!loadGT())
    {
      std::cerr << "Error: Failed to Load Ground Truth!" << std::endl;
      return false;
    };
    if (!loadImages())
    {
      std::cerr << "Error: Failed to Load Images!" << std::endl;
      return false;
    };

    return true;
  }

  SfM_Data GetSfMData() const
  {
    return this->sfm_data_;
  }

  int GetImageNumber() const
  {
    return this->images_.size();
  }
};

#endif // IO_READ_GT_INTERFACE_HPP
