// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IO_LOAD_IMAGES_HPP
#define IO_LOAD_IMAGES_HPP

#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>
#include <vector>

bool LoadImages
(
  const std::string & image_dir,
  const std::vector<std::string> & images,
  const std::vector<openMVG::cameras::PinholeCamera> & cameras,
  openMVG::sfm::SfM_Data & sfm_data
)
{
  if (image_dir.empty() || !stlplus::is_folder(image_dir))
  {
    std::cerr << "Invalid input image directory" << std::endl;
    return false;
  }
  if (images.empty())
  {
    std::cerr << "Invalid input image sequence" << std::endl;
    return false;
  }
  if (cameras.empty())
  {
    std::cerr << "Invalid input camera data" << std::endl;
    return false;
  }

  sfm_data.s_root_path = image_dir; // Setup main image root_path

  Views & views = sfm_data.views;
  Poses & poses = sfm_data.poses;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  C_Progress_display my_progress_bar( images.size(), std::cout, "\n- Loading dataset images -\n" );
  std::ostringstream error_report_stream;
  auto iter_camera = cameras.cbegin();
  for ( auto iter_image = images.cbegin();
    iter_image != images.cend();
    ++iter_image, ++iter_camera, ++my_progress_bar )
  {
    const std::string sImageFilename = stlplus::create_filespec( image_dir, *iter_image );
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

    // Test if this image can be read
    openMVG::image::ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // Image cannot be read

    const Mat3 K = iter_camera->_K;
    const double focal = (K(0,0) + K(1,1)) / 2.0; //Assume K(0,0)==K(1,1)
    const double pxx   = K(0,2);
    const double pyy   = K(1,2);

    const Pose3 pose(iter_camera->_R, iter_camera->_C);
    const auto view      = std::make_shared<sfm::View>(
      *iter_image,
      views.size(), views.size(), views.size(),
      imgHeader.width, imgHeader.height);
    const auto intrinsic = std::make_shared<openMVG::cameras::Pinhole_Intrinsic>(
      imgHeader.width, imgHeader.height,
      focal, pxx, pyy);

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

  // Group the camera that share the same set of camera parameters
  GroupSharedIntrinsics(sfm_data);

  return true;
}

#endif // IO_LOAD_IMAGES_HPP
