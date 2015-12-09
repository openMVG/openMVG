
// Copyright (c) 2012, 2013, 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"

#include "software/SfM/io_readGT.hpp"
#include "software/SfM/tools_precisionEvaluationToGt.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string
    sGTDirectory,
    sComputedDirectory,
    sOutDir = "";
  int camType = -1; //1: openMVG cam, 2,3: Strechas cam


  cmd.add( make_option('i', sGTDirectory, "gt") );
  cmd.add( make_option('c', sComputedDirectory, "computed") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('t', camType, "camtype") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--gt] path (where ground truth camera trajectory are saved)\n"
      << "[-c|--computed] path (openMVG SfM_Output directory)\n"
      << "[-o|--output] path (where statistics will be saved)\n"
      << "[-t|--camtype] Type of the camera:\n"
      << " -1: autoguess (try 1,2,3),\n"
      << "  1: openMVG (bin),\n"
      << "  2: Strechas 'png.camera' \n"
      << "  3: Strechas 'jpg.camera' ]\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //Setup the camera type and the appropriate camera reader
  bool (*fcnReadCamPtr)(const std::string &, PinholeCamera &);
  std::string suffix;

  switch (camType)
  {
    case -1:  // handle auto guess
    {
      if (!stlplus::folder_wildcard(sGTDirectory, "*.bin", false, true).empty())
        camType = 1;
      else if (!stlplus::folder_wildcard(sGTDirectory, "*.png.camera", false, true).empty())
        camType = 2;
      else if (!stlplus::folder_wildcard(sGTDirectory, "*.jpg.camera", false, true).empty())
        camType = 3;
      else
        camType = std::numeric_limits<int>::infinity();
    }
    break;
  }
  switch (camType)
  {
    case 1:
      std::cout << "\nusing openMVG Camera";
      fcnReadCamPtr = &read_openMVG_Camera;
      suffix = "bin";
      break;
    case 2:
    case 3:
      std::cout << "\nusing Strechas Camera";
      fcnReadCamPtr = &read_Strecha_Camera;
      suffix = (camType == 2) ? "png.camera" : "jpg.camera";
      break;
    default:
      std::cerr << "Unsupported camera type. Please write your camera reader." << std::endl;
      return EXIT_FAILURE;
  }

  //---------------------------------------
  // Quality evaluation
  //---------------------------------------

  // Load GT camera rotations & positions [R|C]:

  std::map< std::string, std::pair<Mat3, Vec3> > map_Rt_gt;
  std::map< size_t, PinholeCamera> map_Cam_gt;
  // READ DATA FROM GT
  {
    std::cout << "\nTry to read data from GT";
    std::vector<std::string> vec_fileNames;
    readGt(fcnReadCamPtr, sGTDirectory, suffix, vec_fileNames, map_Rt_gt, map_Cam_gt);
    std::cout << map_Cam_gt.size() << " gt cameras have been found" << std::endl;
  }

  //-- Load the camera that we have to evaluate
  SfM_Data sfm_data;
  if (!Load(sfm_data, sComputedDirectory, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sComputedDirectory << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }
  // Assert that GT and loaded scene have the same camera count
  if (map_Cam_gt.size() != sfm_data.GetPoses().size())
  {
    std::cerr << std::endl
      << "There is missing camera in the loaded scene." << "\n"
      << "#GT poses: " << map_Cam_gt.size() << "\n"
      << "#Scene poses: " << sfm_data.GetPoses().size() << std::endl;
    return EXIT_FAILURE;
  }

  // Prepare data for comparison (corresponding camera centers and rotations)
  Poses::const_iterator iter_loaded_poses = sfm_data.GetPoses().begin();
  std::vector<Vec3> vec_camPosGT, vec_C;
  std::vector<Mat3> vec_camRotGT, vec_camRot;
  for(std::map< size_t, PinholeCamera>::const_iterator iterGT = map_Cam_gt.begin();
    iterGT != map_Cam_gt.end(); ++iterGT, ++iter_loaded_poses)
  {
    // GT
    vec_camPosGT.push_back(iterGT->second._C);
    vec_camRotGT.push_back(iterGT->second._R);

    //-- Computed
    vec_C.push_back(iter_loaded_poses->second.center());
    vec_camRot.push_back(iter_loaded_poses->second.rotation());
  }

  // Visual output of the camera location
  plyHelper::exportToPly(vec_camPosGT, string(stlplus::folder_append_separator(sOutDir) + "camGT.ply").c_str());
  plyHelper::exportToPly(vec_C, string(stlplus::folder_append_separator(sOutDir) + "camComputed.ply").c_str());

  // Evaluation
  htmlDocument::htmlDocumentStream _htmlDocStream("openMVG Quality evaluation.");
  EvaluteToGT(vec_camPosGT, vec_C, vec_camRotGT, vec_camRot, sOutDir, &_htmlDocStream);

  ofstream htmlFileStream( string(stlplus::folder_append_separator(sOutDir) +
    "ExternalCalib_Report.html").c_str());
  htmlFileStream << _htmlDocStream.getDoc();

  return EXIT_SUCCESS;
}

