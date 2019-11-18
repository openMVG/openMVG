// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <tuple>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::image;
using namespace openMVG::sfm;

std::tuple<double, double, double, double, openMVG::geometry::Pose3> LoadCameraParameters(
  const std::string &filename) {
  // fx, fy, ppx, ppy, pose
  std::tuple<double, double, double, double, openMVG::geometry::Pose3> result{
    0, 0, 0, 0, openMVG::geometry::Pose3()};

  std::ifstream is(filename);
  if (!is.good())
    return result;

  auto &pose = std::get<4>(result);
  is >> pose.center().x() >> pose.center().y() >> pose.center().z();
  Quaternion q;
  is >> q.x() >> q.y() >> q.z() >> q.w();
  pose.rotation() = q.toRotationMatrix().transpose();

  is >> std::get<0>(result);
  is >> std::get<1>(result);
  is >> std::get<2>(result);
  is >> std::get<3>(result);

  return result;
}

//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv) {
  CmdLine cmd;

  std::string sImageDir,
    sOutputDir = "";

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  cmd.add(make_option('i', sImageDir, "imageDirectory"));
  cmd.add(make_option('o', sOutputDir, "outputDirectory"));
  cmd.add(make_option('c', i_User_camera_model, "camera_model"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-o|--outputDirectory]\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t 1: Pinhole\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole brown 2\n"
      << "\t 5: Pinhole with a simple Fish-eye distortion\n"
      << "\t 8: Spherical camera\n"
      << "\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << " You called : " << std::endl
    << argv[0] << std::endl
    << "--imageDirectory " << sImageDir << std::endl
    << "--outputDirectory " << sOutputDir << std::endl
    << "--camera_model " << i_User_camera_model << std::endl;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if (!stlplus::folder_exists(sImageDir)) {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty()) {
    std::cerr << "\nInvalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutputDir)) {
    if (!stlplus::folder_create(sOutputDir)) {
      std::cerr << "\nCannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::vector<std::string> vec_image = stlplus::folder_files(sImageDir);
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Poses & poses = sfm_data.poses;

  C_Progress_display my_progress_bar(vec_image.size(),
      std::cout, "\n- Image listing -\n");
  std::ostringstream error_report_stream;
  // Images should be of the same width and height
  int width = 0, height = 0;
  // We are using the same camera model, so accumulate intrinsic values and use their averages
  double accumulated_focal = 0, accumulated_ppx = 0, accumulated_ppy = 0;
  for (auto iter_image = vec_image.begin(); iter_image != vec_image.end(); ++iter_image) {
    ++my_progress_bar;

    const std::string sImageFilename = stlplus::create_filespec(sImageDir, *iter_image);
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown) {
      error_report_stream
        << sImFilenamePart << ": Unkown image file format." << "\n";
      continue; // image cannot be opened
    }

    if (sImFilenamePart.find("mask.png") != std::string::npos
       || sImFilenamePart.find("_mask.png") != std::string::npos) {
      error_report_stream
        << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

    if (width == 0 && height == 0) {
      width = imgHeader.width;
      height = imgHeader.height;
    } else {
      if (width != imgHeader.width || height != imgHeader.height) {
        error_report_stream << sImFilenamePart << " dimensions dismatch\n";
        continue;
      }
    }

    // Test if corresponding camera file exists
    const std::string sCameraFilename = stlplus::create_filespec(sImageDir,
      stlplus::basename_part(sImFilenamePart), ".cam");
    if (!stlplus::file_exists(sCameraFilename)) {
      error_report_stream << sImFilenamePart << " misses camera file\n";
    }

    // Load intrinsic and extrinsic
    double fx = 0, fy = 0, ppx = 0, ppy = 0;
    openMVG::geometry::Pose3 pose;
    std::tie(fx, fy, ppx, ppy, pose) = LoadCameraParameters(sCameraFilename);

    // Build the view corresponding to the image
    View v(*iter_image, views.size(), 0, views.size(), width, height);

    if (fx > 0 && fy > 0 && ppx > 0 && ppy > 0) {
      // Add extrinsic related to the image (if any)
      poses[v.id_pose] = pose;
      // Accumulate intrinsic
      accumulated_focal += fx;
      accumulated_focal += fy;
      accumulated_ppx += ppx;
      accumulated_ppy += ppy;
    }

    // Add the view to the sfm_container
    views[v.id_view] = std::make_shared<View>(v);
  }

  // Build intrinsic parameter
  std::shared_ptr<IntrinsicBase> intrinsic;

  if (accumulated_focal > 0 && accumulated_ppx > 0 && accumulated_ppy > 0
    && width > 0 && height > 0 && views.size() > 0) {
    // Calculate average intrinsic
    double focal = accumulated_focal / views.size() / 2;
    double ppx = accumulated_ppx / views.size();
    double ppy = accumulated_ppy / views.size();

    // Create the desired camera type
    switch (e_User_camera_model) {
    case PINHOLE_CAMERA:
      intrinsic = std::make_shared<Pinhole_Intrinsic>
        (width, height, focal, ppx, ppy);
      break;
    case PINHOLE_CAMERA_RADIAL1:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
        (width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
      break;
    case PINHOLE_CAMERA_RADIAL3:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
        (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
      break;
    case PINHOLE_CAMERA_BROWN:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
        (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
    case PINHOLE_CAMERA_FISHEYE:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
        (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
    case CAMERA_SPHERICAL:
      intrinsic = std::make_shared<Intrinsic_Spherical>
        (width, height);
      break;
    case PINHOLE_CAMERA_RADIAL1_PBA:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1_PBA>
        (width, height, focal, ppx, ppy, 0.0);
      break;
    default:
      std::cerr << "Error: unknown camera model: " << (int)e_User_camera_model << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    std::cerr << "Error: No usable view\n";
    return EXIT_FAILURE;
  }

  // Add intrinsic to the sfm_container
  sfm_data.intrinsics[0] = intrinsic;

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty()) {
    std::cerr
      << "\nWarning & Error messages:" << std::endl
      << error_report_stream.str() << std::endl;
  }

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec(sOutputDir, "sfm_data.json").c_str(),
    ESfM_Data(VIEWS | INTRINSICS | EXTRINSICS))) {
    return EXIT_FAILURE;
  }

  std::cout << std::endl
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;

  return EXIT_SUCCESS;
}
