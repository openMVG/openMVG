// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/pipelines/sequential/sequential_SfM2.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerMaxPair.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerStellar.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

enum class ESfMSceneInitializer
{
  INITIALIZE_EXISTING_POSES,
  INITIALIZE_MAX_PAIR,
  INITIALIZE_AUTO_PAIR,
  INITIALIZE_STELLAR
};

bool StringToEnum_ESfMSceneInitializer
(
  const std::string & str,
  ESfMSceneInitializer & scene_initializer
)
{
  const std::map<std::string, ESfMSceneInitializer> string_to_enum_mapping =
  {
    {"EXISTING_POSE", ESfMSceneInitializer::INITIALIZE_EXISTING_POSES},
    {"MAX_PAIR", ESfMSceneInitializer::INITIALIZE_MAX_PAIR},
    {"AUTO_PAIR", ESfMSceneInitializer::INITIALIZE_AUTO_PAIR},
    {"STELLAR", ESfMSceneInitializer::INITIALIZE_STELLAR},
  };
  auto it  = string_to_enum_mapping.find(str);
  if (it == string_to_enum_mapping.end())
    return false;
  scene_initializer = it->second;
  return true;
}

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Sequential/Incremental reconstruction (Engine 2)" << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir, sMatchFilename;
  std::string sOutDir = "";
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  std::string sSfMInitializer_method = "STELLAR";
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
  bool b_use_motion_priors = false;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('M', sMatchFilename, "match_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_option('S', sSfMInitializer_method, "sfm_initializer") );
  cmd.add( make_switch('P', "prior_usage") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n[Optional]\n"
    << "[-c|--camera_model] Camera model type for view with unknown intrinsic:\n"
      << "\t 1: Pinhole \n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole radial 3 + tangential 2\n"
      << "\t 5: Pinhole fisheye\n"
    << "[-S|--sfm_initializer] Choose the SfM initializer method:\n"
      << "\t 'EXISTING_POSE'-> Initialize the reconstruction from the existing sfm_data camera poses\n"
      << "\t 'MAX_PAIR'-> Initialize the reconstruction from the pair that has the most of matches\n"
      << "\t 'AUTO_PAIR'-> Initialize the reconstruction with a pair selected automatically\n"
      << "\t 'STELLAR'-> Initialize the reconstruction with a 'stellar' reconstruction.\n"
    << "[-f|--refineIntrinsics] Intrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> intrinsic parameters are held as constant\n"
      << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
      << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
      << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
      << "\t -> NOTE: options can be combined thanks to '|'\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
      <<      "\t\t-> refine the focal length & the principal point position\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
      << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
    << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions) (default: false)\n"
    << "[-M|--match_file] path to the match file to use.\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if ( !isValid(openMVG::cameras::EINTRINSIC(i_User_camera_model)) )  {
    std::cerr << "\n Invalid camera type" << std::endl;
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
    cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
  {
    std::cerr << "Invalid input for the Bundle Adjusment Intrinsic parameter refinement option" << std::endl;
    return EXIT_FAILURE;
  }

  ESfMSceneInitializer scene_initializer_enum;
  if (!StringToEnum_ESfMSceneInitializer(sSfMInitializer_method, scene_initializer_enum))
  {
    std::cerr << "Invalid input for the SfM initializer option" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  // Matches reading
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if // Try to read the provided match filename or the default one (matches.f.txt/bin)
  (
    !(matches_provider->load(sfm_data, sMatchFilename) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.f.txt")) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.f.bin")))
  )
  {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "\nCannot create the output directory" << std::endl;
    }
  }

  //---------------------------------------
  // Sequential reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;

  std::unique_ptr<SfMSceneInitializer> scene_initializer;
  switch(scene_initializer_enum)
  {
    case ESfMSceneInitializer::INITIALIZE_AUTO_PAIR:
      std::cerr << "Not yet implemented." << std::endl;
      return EXIT_FAILURE;
    break;
    case ESfMSceneInitializer::INITIALIZE_MAX_PAIR:
      scene_initializer.reset(new SfMSceneInitializerMaxPair(sfm_data,
        feats_provider.get(),
        matches_provider.get()));
    break;
    case ESfMSceneInitializer::INITIALIZE_EXISTING_POSES:
      scene_initializer.reset(new SfMSceneInitializer(sfm_data,
        feats_provider.get(),
        matches_provider.get()));
    break;
    case ESfMSceneInitializer::INITIALIZE_STELLAR:
      scene_initializer.reset(new SfMSceneInitializerStellar(sfm_data,
        feats_provider.get(),
        matches_provider.get()));
    break;
    default:
      return EXIT_FAILURE;
  }
  if (!scene_initializer)
  {
    std::cerr << "Invalid scene initializer." << std::endl;
    return EXIT_FAILURE;
  }

  SequentialSfMReconstructionEngine2 sfmEngine(
    scene_initializer.get(),
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  sfmEngine.SetUnknownCameraType(EINTRINSIC(i_User_camera_model));
  b_use_motion_priors = cmd.used('P');
  sfmEngine.Set_Use_Motion_Prior(b_use_motion_priors);

  if (sfmEngine.Process())
  {
    std::cout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;

    std::cout << "...Generating SfM_Report.html" << std::endl;
    Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

    //-- Export to disk computed scene (data & visualizable results)
    std::cout << "...Export SfM_Data to disk." << std::endl;
    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "sfm_data", ".bin"),
      ESfM_Data(ALL));

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
      ESfM_Data(ALL));

    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
