// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_STRECHA_HPP
#define IO_READ_GT_STRECHA_HPP

#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"

#include <vector>

// The feature of the Strecha's Data:
// 1. each gt file only stores one view's pose information
// 2. the gt file's name is the image file name plus ".camera"
class SfM_Data_GT_Loader_Strecha : public SfM_Data_GT_Loader_Interface
{
private:
  std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
public:
  bool loadGT() override
  {
    // Check all the files under the path
    std::vector<std::string> gt_files = stlplus::folder_wildcard(this->gt_dir_, "*.camera");
    std::sort(gt_files.begin(), gt_files.end());

    if (gt_files.empty())
    {
      std::cerr << "The provided GT directory does not have any *.camera file." << std::endl;
      return false;
    }

    cameras_data_.reserve(gt_files.size());

    // Load the gt_data from the file
    for ( const auto & iter_gt_file : gt_files)
    {
      std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_, iter_gt_file).c_str(), std::ifstream::in);
      if (!gt_file)
      {
        std::cerr << "Error: Failed to open file '" << iter_gt_file << "' for reading" << std::endl;
        continue;
      }
      std::vector<double> val;
      while (gt_file)
      {
        double valT;
        gt_file >> valT;
        val.push_back(valT);
      }
      gt_file.close();

      if (val.size() == 3*3 + 3 +3*3 + 3 + 3 || val.size() == 26) // Strecha cam
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

        cameras_data_.emplace_back(K, R, t);

        // Parse image name
        images_.emplace_back( stlplus::basename_part(stlplus::create_filespec(this->gt_dir_, iter_gt_file)) );
      }
      else
      {
        continue;
      }
    }

    return true;
  }

  bool loadImages() override
  {
    return LoadImages(this->image_dir_, this->images_, this->cameras_data_, this->sfm_data_);
  }
};

#endif // IO_READ_GT_STRECHA_HPP
