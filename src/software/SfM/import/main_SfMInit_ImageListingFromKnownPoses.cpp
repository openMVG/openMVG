// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <utility>

#include "io_readGTBlendedMVS.hpp"
#include "io_readGTDTUMVS.hpp"
#include "io_readGTETH3D.hpp"
#include "io_readGTKitti.hpp"
#include "io_readGTMiddleBury.hpp"
#include "io_readGTStrecha.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir,
    sGroundTruthDir,
    sOutputDir = "";
  int ground_truth_Dataset = 0;


  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('g', sGroundTruthDir, "groundTruthDirectory") );
  cmd.add( make_option('t', ground_truth_Dataset, "groundTruthDataset"));
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );


  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      OPENMVG_LOG_INFO
        << "Usage: " << argv[0] << '\n'
        << "[-i|--imageDirectory]\n"
        << "[-g|--groundTruthDirectory]\n"
        << "[-t]--groundTruthDataset]\n"
        << "\t 1: Strecha's Dataset\n"
        << "\t 2: MiddleBury's Dataset\n"
        << "\t 3: DTU MVS Dataset\n"
        << "\t 4: ETH 3D Dataset\n"
        << "\t 5: Kitti Odometry Dataset\n"
        << "\t 6: BlendedMVS Dataset\n"
        << "[-o|--outputDirectory]\n";

      OPENMVG_LOG_ERROR << s;
      return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO
    << " You called : \n"
    << argv[0] << "\n"
    << "--imageDirectory " << sImageDir << "\n"
    << "--groundTruthDirectory " << sGroundTruthDir << "\n"
    << "--groundTruthDataset " << ground_truth_Dataset << "\n"
    << "--outputDirectory " << sOutputDir;


  if ( !stlplus::folder_exists( sImageDir ) )
  {
    OPENMVG_LOG_ERROR << "\nThe input directory doesn't exist";
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sGroundTruthDir ) )
  {
    OPENMVG_LOG_ERROR << "\nThe input ground truth directory doesn't exist";
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    OPENMVG_LOG_ERROR << "\nInvalid output directory";
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      OPENMVG_LOG_ERROR << "\nCannot create output directory";
      return EXIT_FAILURE;
    }
  }


  // Configure an empty scene with Views and their corresponding cameras
  std::shared_ptr<SfM_Data_GT_Loader_Interface> sfm_data_gt;
  switch(ground_truth_Dataset)
  {
    case 1:
      sfm_data_gt = std::make_shared<SfM_Data_GT_Loader_Strecha>();
      break;
    case 2:
      sfm_data_gt = std::make_shared<SfM_Data_GT_Loader_MiddleBury>();
      break;
    case 3:
      sfm_data_gt = std::make_shared<SfM_Data_GT_Loader_DTU_MVS>();
      break;
    case 4:
      sfm_data_gt = std::make_shared<SfM_Data_GT_Loader_ETH_3D>();
      break;
    case 5:
      sfm_data_gt = std::make_shared<SfM_Data_GT_Loader_Kitti>();
      break;
    case 6:
      sfm_data_gt = std::make_shared<SfM_Data_GT_Loader_BlendedMVS>();
      break;
    default:
      OPENMVG_LOG_ERROR << "Error: Not Support Dataset";
      return EXIT_FAILURE;
  }

  // Load the gt data
  if (!sfm_data_gt->run(sGroundTruthDir, sImageDir))
  {
    return EXIT_FAILURE;
  }

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data_gt->GetSfMData(),
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO
    << "\nSfMInit_ImageListing report:\n"
    << "listed #File(s): " << sfm_data_gt->GetImageNumber() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data_gt->GetSfMData().GetViews().size() << "\n"
    << "usable #Pose(s) listed in sfm_data: " << sfm_data_gt->GetSfMData().GetPoses().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data_gt->GetSfMData().GetIntrinsics().size();
  return EXIT_SUCCESS;
}
