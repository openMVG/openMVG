// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stlMap.hpp"

#include "software/SfM/tools_precisionEvaluationToGt.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

#include "software/SfM/import/io_readGTInterface.hpp"
#include "software/SfM/import/io_readGTStrecha.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string
    sGTDirectory,
    sComputedDirectory,
    sOutDir = "";

  cmd.add( make_option('i', sGTDirectory, "gt") );
  cmd.add( make_option('c', sComputedDirectory, "computed") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--gt] path (where the ground truth dataset is saved).\n"
      << "[-c|--computed] path to the sfm_data file to compare.\n"
      << "[-o|--output] path (where the registration statistics will be saved).\n"
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

  //---------------------------------------
  // Quality evaluation
  //---------------------------------------

  // 1. Initialize the GT:
  //
  std::shared_ptr<SfM_Data_GT_Loader_Interface> sfm_data_gt_loader =
      std::make_shared<SfM_Data_GT_Loader_Strecha>();

  // Load the gt data
  if (!sfm_data_gt_loader->run(
        stlplus::folder_append_separator(sGTDirectory) + "gt_dense_cameras",
        stlplus::folder_append_separator(sGTDirectory) + "images"))
  {
    std::cerr << "Cannot initialize the GT dataset." << std::endl;
    return EXIT_FAILURE;
  }
  const SfM_Data sfm_data_gt = sfm_data_gt_loader->GetSfMData();

  // 2. Load the scene to compare
  //
  SfM_Data sfm_data_to_compare;
  if (!Load(sfm_data_to_compare, sComputedDirectory, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sComputedDirectory << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }
  // Assert that GT and loaded scene have the same camera count
  if (sfm_data_gt.GetPoses().size() != sfm_data_to_compare.GetPoses().size())
  {
    std::cerr << std::endl
      << "Dataset poses count does not match." << "\n"
      << "#GT poses: " << sfm_data_gt.GetPoses().size() << "\n"
      << "#Scene poses: " << sfm_data_to_compare.GetPoses().size() << std::endl;
    return EXIT_FAILURE;
  }

  // Assert that the sfm_data files have the same camera's id_pose.
  {
    // Collect GT pose_ids.
    std::set<IndexT> pose_id_gt;
    std::transform(
      sfm_data_gt.GetPoses().cbegin(),
      sfm_data_gt.GetPoses().cend(),
      std::inserter(pose_id_gt, pose_id_gt.begin()),
      stl::RetrieveKey());

    // Collect the sfm_data to compare pose_ids.
    std::set<IndexT> pose_id_to_compare;
    std::transform(
      sfm_data_to_compare.GetPoses().cbegin(),
      sfm_data_to_compare.GetPoses().cend(),
      std::inserter(pose_id_to_compare, pose_id_to_compare.begin()),
      stl::RetrieveKey());

    // Check if the pose_id intersect or not
    std::vector<IndexT> pose_id_intersection;
    std::set_intersection(pose_id_gt.cbegin(), pose_id_gt.cend(),
                          pose_id_to_compare.cbegin(), pose_id_to_compare.cend(),
                          std::back_inserter(pose_id_intersection));

    if (pose_id_gt.empty() || pose_id_to_compare.empty()
        || pose_id_intersection.size() != sfm_data_gt.GetPoses().size())
    {
      std::cerr << "Invalid data input. "
       << "The dataset does not have corresponding camera Ids." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // 3. Find corresponding camera pose data:
  //
  std::vector<Vec3> camera_pos_gt, camera_pos_to_compare;
  std::vector<Mat3> camera_rot_gt, camera_rot_to_compare;

  for (const auto & pose_gt_it : sfm_data_gt.GetPoses())
  {
    const IndexT pose_id = pose_gt_it.first;
    //sfm_data_gt.GetPoses().count(pose_id)
    const auto & pose_to_compare_it = sfm_data_gt.GetPoses().at(pose_id);

    camera_pos_gt.push_back(pose_gt_it.second.center());
    camera_rot_gt.push_back(pose_gt_it.second.rotation());

    camera_pos_to_compare.push_back(pose_to_compare_it.center());
    camera_rot_to_compare.push_back(pose_to_compare_it.rotation());
  }

  // Visual output of the camera location
  plyHelper::exportToPly(camera_pos_gt, string(stlplus::folder_append_separator(sOutDir) + "camGT.ply").c_str());
  plyHelper::exportToPly(camera_pos_to_compare, string(stlplus::folder_append_separator(sOutDir) + "camComputed.ply").c_str());

  // Evaluation
  htmlDocument::htmlDocumentStream _htmlDocStream("openMVG Quality evaluation.");
  EvaluteToGT(camera_pos_gt, camera_pos_to_compare,
    camera_rot_gt, camera_rot_to_compare,
    sOutDir, &_htmlDocStream);

  ofstream htmlFileStream( string(stlplus::folder_append_separator(sOutDir) +
    "ExternalCalib_Report.html").c_str());
  htmlFileStream << _htmlDocStream.getDoc();

  return EXIT_SUCCESS;
}
