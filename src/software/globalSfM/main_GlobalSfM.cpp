
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/system/timer.hpp"
using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Global Structure from Motion:\n"
    << "-----------------------------------------------------------\n"
    << "Open Source implementation of the paper:\n"
    << "\"Global Fusion of Relative Motions for "
    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
    << " ICCV 2013." << std::endl
    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  bool bColoredPointCloud = false;
  int iRotationAveragingMethod = 2;
  int iTranslationAveragingMethod = 1;
  bool bRefineIntrinsics = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('c', bColoredPointCloud, "coloredPointCloud") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', bRefineIntrinsics, "refineIntrinsics") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "[-c|--coloredPointCloud 0(default) or 1]\n"
    << "[-r|--rotationAveraging 2(default L2) or 1 (L1)]\n"
    << "[-t|--translationAveraging 1(default L1) or 2 (L2)]\n"
    << "[-f|--refineIntrinsics \n"
    << "\t 0-> intrinsic parameters are kept as constant\n"
    << "\t 1-> refine intrinsic parameters (default).] \n"
    << "\n"
    << " ICCV 2013: => -r 2 -t 1"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (iRotationAveragingMethod < globalSfM::ROTATION_AVERAGING_L1 ||
      iRotationAveragingMethod > globalSfM::ROTATION_AVERAGING_L2 )  {
    std::cerr << "\n Rotation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

    if (iTranslationAveragingMethod < globalSfM::TRANSLATION_AVERAGING_L1 ||
      iTranslationAveragingMethod > globalSfM::TRANSLATION_AVERAGING_L2 )  {
    std::cerr << "\n Translation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Prepare the features and matches provider
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }

  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.txt"))) {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------

  openMVG::Timer timer;
  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_bFixedIntrinsics(!bRefineIntrinsics);

  if (sfmEngine.Process())
  {
    std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;

    //-- Export to disk computed scene (data & visualizable results)
    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "sfm_data", ".json"),
      ESfM_Data(ALL));

    Save(sfmEngine.Get_SfM_Data(),
      stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
      ESfM_Data(ALL));

    if (bColoredPointCloud)
    {
      std::cout << "Point cloud colorization is not yet implemented." << std::endl;
    }

    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
