
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

bool computeIndexFromImageName(
  const SfM_Data & sfm_data,
  const std::string& initialName,
  IndexT& initialIndex)
{  
  initialIndex = UndefinedIndexT;

  bool isName = (initialName == stlplus::filename_part(initialName));
  
  /// List views filenames and find the one that correspond to the user ones:
  for(Views::const_iterator it = sfm_data.GetViews().begin();
    it != sfm_data.GetViews().end(); ++it)
  {
    const View * v = it->second.get();
    std::string filename;
    
    if(isName)
    {
      filename = stlplus::filename_part(v->s_Img_path);
    }
    else
    {
      if(stlplus::is_full_path(v->s_Img_path))
      {
        filename = v->s_Img_path;
      }
      else
      {
        filename = sfm_data.s_root_path + v->s_Img_path;
      }
    }
    
    if (filename == initialName)
    {
      if(initialIndex == UndefinedIndexT)
      {
          initialIndex = v->id_view;
      }
      else
      {
        std::cout<<"Error: Two pictures named :" << initialName << " !" << std::endl;
      }
    }
  }
  return initialIndex != UndefinedIndexT;
}


int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Sequential/Incremental reconstruction" << std::endl
            << " Perform incremental SfM (Initial Pair Essential + Resection)." << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  std::string sOutSfMDataFilepath = "";
  std::string sOutInterFileExtension = ".ply";
  std::pair<std::string,std::string> initialPairString("","");
  bool bRefineIntrinsics = true;
  int minInputTrackLength = 2;
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
  bool allowUserInteraction = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('s', sOutSfMDataFilepath, "out_sfmdata_file") );
  cmd.add( make_option('e', sOutInterFileExtension, "inter_file_extension") );
  cmd.add( make_option('a', initialPairString.first, "initialPairA") );
  cmd.add( make_option('b', initialPairString.second, "initialPairB") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('f', bRefineIntrinsics, "refineIntrinsics") );
  cmd.add( make_option('t', minInputTrackLength, "minInputTrackLength") );
  cmd.add( make_option('u', allowUserInteraction, "allowUserInteraction") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "[-s|--out_sfmdata_file] path of the output sfmdata file (default: $outdir/sfm_data.json)\n"
    << "[-e|--inter_file_extension] extension of the intermediate file export (default: .ply)\n"
    << "[-a|--initialPairA] filename of the first image (without path)\n"
    << "[-b|--initialPairB] filename of the second image (without path)\n"
    << "[-c|--camera_model] Camera model type for view with unknown intrinsic:\n"
      << "\t 1: Pinhole \n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
    << "[-f|--refineIntrinsics] \n"
    << "\t 0-> intrinsic parameters are kept as constant\n"
    << "\t 1-> refine intrinsic parameters (default). \n"
    << "[-t|--minInputTrackLength N] minimum track length in input of SfM (default: 2)\n"
    << "[-p|--matchFilePerImage] \n"
    << "\t To use one match file per image instead of a global file.\n"
    << "[-u|--allowUserInteraction] Enable/Disable user interactions. (default: true)\n"
    << "\t If the process is done on renderfarm, it doesn't make sense to wait for user inputs.\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if(sOutSfMDataFilepath.empty())
    sOutSfMDataFilepath = stlplus::create_filespec(sOutDir, "sfm_data", "json");

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
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

  if(!matches_provider->load(sfm_data, sMatchesDir, "f"))
  {
    std::cerr << std::endl << "Unable to load matches file from: " << sMatchesDir << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())
  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Sequential reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;
  SequentialSfMReconstructionEngine sfmEngine(
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_bFixedIntrinsics(!bRefineIntrinsics);
  sfmEngine.SetUnknownCameraType(EINTRINSIC(i_User_camera_model));
  sfmEngine.setMinInputTrackLength(minInputTrackLength);
  sfmEngine.setSfmdataInterFileExtension(sOutInterFileExtension);
  sfmEngine.setAllowUserInteraction(allowUserInteraction);

  if (initialPairString.first == initialPairString.second)
  {
    std::cerr << "\nInvalid image names. You cannot use the same image to initialize a pair." << std::endl;
    return false;
  }
  
  // Handle Initial pair parameter
  if (!initialPairString.first.empty() && !initialPairString.second.empty())
  {
    Pair initialPairIndex;
    if(!computeIndexFromImageName(sfm_data, initialPairString.first, initialPairIndex.first)
            || !computeIndexFromImageName(sfm_data, initialPairString.second, initialPairIndex.second))
    {
        std::cerr << "Could not find the initial pairs <" << initialPairString.first
          <<  ", " << initialPairString.second << ">!\n";
        return EXIT_FAILURE;
    }
 
    sfmEngine.setInitialPair(initialPairIndex);
  }

  if (!sfmEngine.Process())
  {
    return EXIT_FAILURE;
  }

  // get the color for the 3D points
  sfmEngine.Colorize();

  std::cout << std::endl << " Total Ac-Sfm took (s): " << timer.elapsed() << std::endl;

  std::cout << "...Generating SfM_Report.html" << std::endl;
  Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
    stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

  //-- Export to disk computed scene (data & visualizable results)
  std::cout << "...Export SfM_Data to disk:" << std::endl;
  std::cout << "   " << sOutSfMDataFilepath << std::endl;

  Save(sfmEngine.Get_SfM_Data(), sOutSfMDataFilepath, ESfM_Data(ALL));

  Save(sfmEngine.Get_SfM_Data(), stlplus::create_filespec(sOutDir, "cloud_and_poses", sOutInterFileExtension), ESfM_Data(VIEWS | EXTRINSICS | INTRINSICS | STRUCTURE));

  return EXIT_SUCCESS;
}
