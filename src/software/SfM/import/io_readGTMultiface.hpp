// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"



class SfM_Data_GT_Loader_Multiface : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
public:
  bool loadGT() override
  {
    this->images_ = stlplus::folder_wildcard(this->image_dir_, "*.png");
    std::sort(this->images_.begin(), this->images_.end());

    const size_t nImages = images_.size();

    cameras_data_.resize(nImages);

    std::ifstream krt_file(this->gt_dir_ + "/KRT");
    if (!krt_file)
    {
      OPENMVG_LOG_ERROR << "File << " << (this->gt_dir_ + "/KRT") << " could not be found.";
      return false;
    }

    std::istringstream iss;

    for (size_t i = 0; i < nImages; ++i)
    {
      std::string line[9];
      for (int j = 0; j < 9; ++j)
        std::getline(krt_file, line[j]);

      bool found = false;
      for (size_t j = 0; j < nImages; ++j)
      {
        if (stlplus::basename_part(this->images_[j]) == line[0])
        {
          Mat3 R, K;
          Vec3 t;
          iss = std::istringstream{line[1]}; iss >> K(0, 0) >> K(0, 1) >> K(0, 2);
          iss = std::istringstream{line[2]}; iss >> K(1, 0) >> K(1, 1) >> K(1, 2);
          iss = std::istringstream{line[3]}; iss >> K(2, 0) >> K(2, 1) >> K(2, 2);
          iss = std::istringstream{line[5]}; iss >> R(0, 0) >> R(0, 1) >> R(0, 2) >> t[0];
          iss = std::istringstream{line[6]}; iss >> R(1, 0) >> R(1, 1) >> R(1, 2) >> t[1];
          iss = std::istringstream{line[7]}; iss >> R(2, 0) >> R(2, 1) >> R(2, 2) >> t[2];

          cameras_data_[j] = cameras::PinholeCamera(K, R, t);
          found = true;
          break;
        }
      }
      if (!found)
      {
        OPENMVG_LOG_ERROR << "Camera " << line[0] << " could not find its corresponding image.";
        return false;
      }
    }

    return true;
  }

  bool loadImages() override
  {
    return LoadImages(this->image_dir_, this->images_, this->cameras_data_, this->sfm_data_);
  }
};
