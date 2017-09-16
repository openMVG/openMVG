// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/stl/split.hpp"

#include "openMVG/sfm/sfm.hpp"

#include "software/SfM/SfMIOHelper.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

// Compatibility binary
// Convert an openMVG v0.x=>7 lists.txt file to the new SfM_Data file format (v0.8)
// - Export a SfM_Data file with View & Intrinsic data
//

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir,
    sOutputDir = "";

  int i_User_camera_model = PINHOLE_CAMERA;
  bool b_Group_camera_model = true;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-o|--outputDirectory]\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t 1: Pinhole (default)\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3\n"
      << "\t 4: Pinhole brown with radial 3 and tangential 2\n"
      << "\t 5: Pinhole with a simple Fish-eye distortion\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> views can share some camera intrinsic parameters (default)\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    std::cerr << "\nInvalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      std::cerr << "\nCannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Check if the sOutputDir/lists.txt file is present
  const std::string sListsFile = stlplus::create_filespec( sOutputDir, "lists.txt" );
  if ( !stlplus::is_file( sListsFile ) )
  {
    std::cerr << "\nThe input lists.txt file doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  // Parse the lists.txt file and create the corresponding View & Intrinsic data

  std::vector<openMVG::SfMIO::CameraInfo> vec_camImageNames;
  std::vector<openMVG::SfMIO::IntrinsicCameraInfo> vec_intrinsicGroups;
  if (!openMVG::SfMIO::loadImageList( vec_camImageNames,
                                      vec_intrinsicGroups,
                                      sListsFile) )
  {
    std::cerr << "\nEmpty image list." << std::endl;
    return false;
  }

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;


  // Find to which intrinsic groups each image belong
  for (std::vector<openMVG::SfMIO::CameraInfo>::const_iterator iter = vec_camImageNames.begin();
    iter != vec_camImageNames.end(); ++iter)
  {
    const openMVG::SfMIO::CameraInfo & camInfo = *iter;
    // Find the index of the corresponding cameraInfo
    const size_t idx = std::distance((std::vector<openMVG::SfMIO::CameraInfo>::const_iterator)vec_camImageNames.begin(), iter);

    // Expected properties for each image
    double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

    std::shared_ptr<IntrinsicBase> intrinsic;

    const openMVG::SfMIO::IntrinsicCameraInfo & camIntrinsic = vec_intrinsicGroups[camInfo.m_intrinsicId];
    width = camIntrinsic.m_w;
    height = camIntrinsic.m_h;
    if (camIntrinsic.m_bKnownIntrinsic)
    {
      ppx = camIntrinsic.m_K(0,2);
      ppy = camIntrinsic.m_K(1,2);
      focal = camIntrinsic.m_focal;
      // Create the user defined camera model corresponding to the intrinsic loaded data
      switch (e_User_camera_model)
      {
        case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>(width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_RADIAL1:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>(width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_RADIAL3:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>(width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_BROWN:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>(width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_FISHEYE:
        intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>(width, height, focal, ppx, ppy);
        break;
        default:
          std::cerr << "Invalid camera model." << std::endl;
          return EXIT_FAILURE;
      }
    }
    // else -> no intrinsic data

    // Create the corresponding view
    const IndexT id_view = idx;
    const IndexT id_intrinsic =
      (camIntrinsic.m_bKnownIntrinsic) ?
        (b_Group_camera_model ? camInfo.m_intrinsicId : id_view)
        : UndefinedIndexT;
    const IndexT id_pose = id_view;

    // Build the view corresponding to the image
    View v( iter->m_sImageName, id_view, id_intrinsic, id_pose, width, height);
    // Add the view to the sfm_container
    views[v.id_view] = std::make_shared<View>(v);

    // Add intrinsic related to the image (if any)
    if (intrinsic == nullptr)
    {
      //Since the view have invalid intrinsic data
      // (export the view, with an invalid intrinsic field value)
      v.id_intrinsic = UndefinedIndexT;
    }
    else
    {
      // Add the intrinsic to the sfm_container
      intrinsics[v.id_intrinsic] = intrinsic;
    }
  }

  // Store SfM_Data views & intrinsic data
  if (Save(
      sfm_data,
      stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
      ESfM_Data(VIEWS|INTRINSICS)))
  {
  return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
