
// Copyright (c) 2012, 2013, 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"

#include "software/SfM/tools_precisionEvaluationToGt.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include <cstdlib>
#include <iostream>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string
    sGTFile,
    sComputedFile,
    sOutDirectory = "";


  cmd.add( make_option('i', sGTFile, "gt") );
  cmd.add( make_option('c', sComputedFile, "computed") );
  cmd.add( make_option('o', sOutDirectory, "outdir") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--gt] path (where ground truth camera trajectory are saved)\n"
      << "[-c|--computed] path (openMVG sfm_data.json file)\n"
      << "[-o|--output] path (where statistics will be saved)\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDirectory.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDirectory))
    stlplus::folder_create(sOutDirectory);

  //---------------------------------------
  // Quality evaluation
  //---------------------------------------

  //-- Load GT camera rotations & positions [R|C]:
  SfM_Data sfm_data_gt;
  // READ DATA FROM GT
  std::cout << "Try to read data from GT" << std::endl;
  if (!Load(sfm_data_gt, sGTFile, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS)))
  {
    std::cerr << "The input SfM_Data file \""<< sComputedFile << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << sfm_data_gt.poses.size() << " gt cameras have been found" << std::endl;

  //-- Load the camera that we have to evaluate
  SfM_Data sfm_data;
  if (!Load(sfm_data, sComputedFile, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sComputedFile << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Fill vectors of valid views for evaluation
  std::vector<Vec3> vec_camPosGT, vec_C;
  std::vector<Mat3> vec_camRotGT, vec_camRot;
  for(const auto &iter : sfm_data.GetViews())
  {
    const auto &view = iter.second;
    // Jump to next view if there is no correponding pose in reconstruction
    if(sfm_data.GetPoses().find(view->id_pose) == sfm_data.GetPoses().end())
    {
      std::cout << "no pose in input (" << view->id_pose << ")" << std::endl;
      continue;
    }

    // Jump to next view if there is no corresponding view in GT
    if(sfm_data_gt.GetViews().find(view->id_view) == sfm_data_gt.GetViews().end())
    {
      std::cout << "no view in GT (" << view->id_view << ")" << std::endl;
      continue;
    }
    const int idPoseGT = sfm_data_gt.GetViews().at(view->id_view)->id_pose;

    //-- GT
    const geometry::Pose3 pose_gt = sfm_data_gt.GetPoses().at(idPoseGT);
    vec_camPosGT.push_back(pose_gt.center());
    vec_camRotGT.push_back(pose_gt.rotation());

    //-- Data to evaluate
    const geometry::Pose3 pose_eval = sfm_data.GetPoses().at(view->id_pose);
    vec_C.push_back(pose_eval.center());
    vec_camRot.push_back(pose_eval.rotation());
  }

  // Visual output of the camera location
  plyHelper::exportToPly(vec_camPosGT, string(stlplus::folder_append_separator(sOutDirectory) + "camGT.ply").c_str());
  plyHelper::exportToPly(vec_C, string(stlplus::folder_append_separator(sOutDirectory) + "camComputed.ply").c_str());

  // Evaluation
  htmlDocument::htmlDocumentStream _htmlDocStream("openMVG Quality evaluation.");
  EvaluteToGT(vec_camPosGT, vec_C, vec_camRotGT, vec_camRot, sOutDirectory, &_htmlDocStream);

  ofstream htmlFileStream( string(stlplus::folder_append_separator(sOutDirectory) +
    "ExternalCalib_Report.html"));
  htmlFileStream << _htmlDocStream.getDoc();

  return EXIT_SUCCESS;
}

