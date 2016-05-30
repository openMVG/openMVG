#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "software/SfM/io_readGT.hpp"
#include "software/SfM/tools_precisionEvaluationToGt.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include <vector>
#include <iostream>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  // Command line configuration
  //
  CmdLine cmd;

  std::string
    sSfmFile,
    sGTDirectory,
    sOutFile = "";
  bool mayaTransform = false;


  cmd.add( make_option('i', sSfmFile, "sfm") );
  cmd.add( make_option('g', sGTDirectory, "gt") );
  cmd.add( make_option('o', sOutFile, "output") );
  cmd.add( make_option('m', mayaTransform, "maya") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfm] path (OpenMVG Json file, ouput from SfMInit_ImageListing)\n"
      << "[-g|--gt] path (where ground truth camera trajectory are saved)\n"
      << "[-o|--output] path (where export file, .json or .abc, will be save)\n"
      << "[-m|--maya] enable [1,-1,-1] transformation if ground truth comes from Maya\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Setup the camera type and the appropriate camera reader
  //
  bool (*fcnReadCamPtr)(const std::string&, Pinhole_Intrinsic&, geometry::Pose3&);
  int camType;
  std::string suffix;

  if (!stlplus::folder_wildcard(sGTDirectory, "*.bin", false, true).empty())
    camType = 1;
  else if (!stlplus::folder_wildcard(sGTDirectory, "*.png.camera", false, true).empty())
    camType = 2;
  else if (!stlplus::folder_wildcard(sGTDirectory, "*.jpg.camera", false, true).empty())
    camType = 3;
  else if (!stlplus::folder_wildcard(sGTDirectory, "*.PNG.camera", false, true).empty())
    camType = 4;
  else if (!stlplus::folder_wildcard(sGTDirectory, "*.JPG.camera", false, true).empty())
    camType = 5;
  else
    camType = std::numeric_limits<int>::infinity();
  switch (camType)
  {
    case 1:
      std::cout << "\nusing openMVG Camera";
      fcnReadCamPtr = &read_openMVG_Camera;
      suffix = "bin";
      break;
    case 2:
      std::cout << "\nusing Strechas Camera (png)";
      fcnReadCamPtr = &read_Strecha_Camera;
      suffix = "png.camera";
      break;
    case 3:
      std::cout << "\nusing Strechas Camera (jpg)";
      fcnReadCamPtr = &read_Strecha_Camera;
      suffix = "jpg.camera";
      break;
    case 4:
      std::cout << "\nusing Strechas Camera (PNG)";
      fcnReadCamPtr = &read_Strecha_Camera;
      suffix = "PNG.camera";
      break;
    case 5:
      std::cout << "\nusing Strechas Camera (JPG)";
      fcnReadCamPtr = &read_Strecha_Camera;
      suffix = "JPG.camera";
      break;
    default:
      std::cerr << "Unsupported camera type. Please write your camera reader." << std::endl;
      return EXIT_FAILURE;
  }

  // Load JSON
  //
  SfM_Data sfm_data_in;
  if (!Load(sfm_data_in, sSfmFile, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS)))
  {
    std::cerr << "ERROR" << std::endl;
    std::cerr << "The input SfM_Data file \"" << sSfmFile << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Load GT
  //
  SfM_Data sfm_data_gt;
  std::vector<std::string> vec_fileNames;
  readGt(fcnReadCamPtr, sGTDirectory, suffix, vec_fileNames, sfm_data_gt.poses, sfm_data_gt.intrinsics);
  std:cout << sfm_data_gt.poses.size() << " gt cameras have been found" << std::endl;
  
  // Fill sfm_data_in with poses from sfm_data_gt
  //
  for(const auto &iter : sfm_data_in.GetViews())
  {
    const auto &view = iter.second;
    const std::string sImageName = stlplus::filename_part(view->s_Img_path);
    const int idGT = findIdGT(sImageName, vec_fileNames);
    if(idGT == -1)
    {
      std::cerr << "No ground truth for file: " << sImageName << std::endl;
      continue;
    }
    else
    {
      const geometry::Pose3 poseGT = sfm_data_gt.GetPoses().at(idGT);
      Vec3 vecMaya;
      if(mayaTransform)
        vecMaya = {1,-1,-1};
      else
        vecMaya = {1,1,1};
      geometry::Pose3 poseIN(vecMaya.asDiagonal() * poseGT.rotation(), poseGT.center());

      sfm_data_in.poses.emplace(view->id_pose, poseIN);
      sfm_data_in.intrinsics[view->id_intrinsic] = sfm_data_gt.intrinsics.at(idGT);
    }
  }

  std::cout << "Saved: " << Save(sfm_data_in, sOutFile, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS)) << std::endl;

  return EXIT_SUCCESS;
}

