// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_READ_GT_DTUMVS_HPP
#define IO_READ_GT_DTUMVS_HPP


#include "io_readGTInterface.hpp"
#include "io_loadImages.hpp"

#include "openMVG/cameras/PinholeCamera.hpp"

// The feature of the DTU_MVS's Data:
// 1. each gt file only stores one view's pose information
// 2. the gt file's name contains number that matches with the image
// 3. the gt files' number(64) doesn't match the number of the images(49,64)
class SfM_Data_GT_Loader_DTU_MVS : public SfM_Data_GT_Loader_Interface
{
private:
    std::vector<cameras::PinholeCamera> cameras_data_; // Store all the camera information
public:
  bool loadGT() override
  {
    // Check all the files under the path
    std::vector<std::string> gt_files = stlplus::folder_files( this->gt_dir_ );
    std::sort(gt_files.begin(), gt_files.end());

    // Check the image_dir to find the image names
    std::vector<std::string> image_files = stlplus::folder_files( this->image_dir_ );
    std::sort(image_files.begin(), image_files.end());

    cameras_data_.reserve(gt_files.size());

    // Load the gt_data from the file
    for ( const auto & gt_file_it : gt_files)
    {
      std::ifstream gt_file( stlplus::create_filespec(this->gt_dir_, gt_file_it), std::ifstream::in);
      if (!gt_file)
      {
        std::cerr << "Error: Failed to open file '" << gt_file_it << "' for reading" << std::endl;
        continue;
      }

      std::vector<double> val;
      while (gt_file)
      {
        double valT;
        gt_file >> valT;
        if (!gt_file.fail())
          val.push_back(valT);
      }

      gt_file.close();
      if (val.size() == 12)
      {
        Mat34 P;
        P << val[0], val[1], val[2], val[3],
             val[4], val[5], val[6], val[7],
             val[8], val[9], val[10], val[11];

        cameras_data_.emplace_back(P);

        // Parse image name
        const std::string index = stlplus::basename_part(stlplus::create_filespec(this->gt_dir_, gt_file_it)).substr(4, 3);

        // Check image exist
        for ( const auto image_file_it : image_files)
        {
          const std::string image_index = image_file_it.substr(5,3); // Match code
          if (image_index == index)
          {
            images_.push_back( image_file_it );
            break;
          }
        }
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


#endif // IO_READ_GT_DTUMVS_HPP
