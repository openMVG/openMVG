// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/stl/stlMap.hpp"

#include "software/SfM/tools_precisionEvaluationToGt.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

#include "software/SfM/import/io_readGTStrecha.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cstdlib>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string
    sGTPath,
    sComputedSfmDataFilePath,
    sOutDir = "";

  cmd.add( make_option('i', sGTPath, "gt") );
  cmd.add( make_option('c', sComputedSfmDataFilePath, "computed") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
      << "[-i|--gt] path (where the ground truth dataset is saved (i.e sfm_data.bin/json/xml).\n"
      << "[-c|--computed] path to the sfm_data file to compare (i.e sfm_data.bin/json/xml).\n"
      << "[-o|--output] path (where the registration statistics will be saved).\n";

    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    OPENMVG_LOG_ERROR << "\nIt is an invalid output directory";
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Quality evaluation
  //---------------------------------------

  // 1. Initialize the GT:
  //
  SfM_Data sfm_data_gt;
  {
    if (stlplus::is_file(sGTPath))
    {
      if (!Load(sfm_data_gt, sGTPath, ESfM_Data(VIEWS|EXTRINSICS))) {
        OPENMVG_LOG_ERROR << "The input GT path file \""<< sGTPath << "\" cannot be read.";
        return EXIT_FAILURE;
      }
    }
    else // Compatibility with previous configuration (Strecha dataset)
    {
      if (stlplus::is_folder(sGTPath))
      {
        SfM_Data_GT_Loader_Strecha sfm_data_gt_loader;
        if (!sfm_data_gt_loader.run(
             stlplus::folder_append_separator(sGTPath) + "gt_dense_cameras",
             stlplus::folder_append_separator(sGTPath) + "images"))
        {
          OPENMVG_LOG_ERROR << "The input GT path directory \""<< sGTPath << "\" cannot be read.";
          return EXIT_FAILURE;
        }
        sfm_data_gt = sfm_data_gt_loader.GetSfMData();
      }
      else
      {
        OPENMVG_LOG_ERROR << "Cannot read GT: "<< sGTPath;
        return EXIT_FAILURE;
      }
    }
  }


  // 2. Load the scene to compare
  //
  SfM_Data sfm_data_to_compare;
  {
    if (!Load(sfm_data_to_compare, sComputedSfmDataFilePath, ESfM_Data(VIEWS|EXTRINSICS))) {
      OPENMVG_LOG_ERROR << "The input SfM_Data file \""<< sComputedSfmDataFilePath << "\" cannot be read.";
      return EXIT_FAILURE;
    }
  }

  // Check if GT and loaded scene have the same camera count
  if (sfm_data_gt.GetPoses().size() != sfm_data_to_compare.GetPoses().size())
  {
    OPENMVG_LOG_WARNING
      << "Dataset poses count does not match." << "\n"
      << "#GT poses: " << sfm_data_gt.GetPoses().size() << "\n"
      << "#Scene poses: " << sfm_data_to_compare.GetPoses().size();
  }

  // Collect corresponding poses ids
  std::vector<IndexT> pose_id_intersection;
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
    std::set_intersection(pose_id_gt.cbegin(), pose_id_gt.cend(),
                          pose_id_to_compare.cbegin(), pose_id_to_compare.cend(),
                          std::back_inserter(pose_id_intersection));

    if (pose_id_intersection.size() != sfm_data_gt.GetPoses().size())
    {
      OPENMVG_LOG_WARNING << "The datasets' poses id count does not match.";
    }
  }

  if (pose_id_intersection.empty())
  {
    OPENMVG_LOG_WARNING << "No common poses Id found.";
  }

  // 3. Find corresponding camera pose data:
  //
  std::vector<Vec3> camera_pos_gt, camera_pos_to_compare;
  std::vector<Mat3> camera_rot_gt, camera_rot_to_compare;

  for (const auto & pose_id : pose_id_intersection)
  {
    const auto & pose_to_compare_it = sfm_data_to_compare.GetPoses().at(pose_id);
    const auto & pose_gt_it = sfm_data_gt.GetPoses().at(pose_id);

    camera_pos_gt.push_back(pose_gt_it.center());
    camera_rot_gt.push_back(pose_gt_it.rotation());

    camera_pos_to_compare.push_back(pose_to_compare_it.center());
    camera_rot_to_compare.push_back(pose_to_compare_it.rotation());
  }

  // Visual output of the camera location
  plyHelper::exportToPly(camera_pos_gt,
    string(stlplus::folder_append_separator(sOutDir) + "camGT.ply").c_str());
  plyHelper::exportToPly(camera_pos_to_compare,
    string(stlplus::folder_append_separator(sOutDir) + "camComputed.ply").c_str());

  // Evaluation
  htmlDocument::htmlDocumentStream _htmlDocStream("openMVG Quality evaluation.");
  std::vector<double> vec_distance_residuals, vec_rotation_angular_residuals;
  EvaluateToGT(camera_pos_gt, camera_pos_to_compare,
    camera_rot_gt, camera_rot_to_compare,
    sOutDir, &_htmlDocStream, vec_distance_residuals, vec_rotation_angular_residuals);

  ofstream htmlFileStream( string(stlplus::folder_append_separator(sOutDir) +
    "ExternalCalib_Report.html").c_str());
  htmlFileStream << _htmlDocStream.getDoc();


  // Export a JSON blob with all the statistics:
  const std::string out_json_filename =
    stlplus::create_filespec(sOutDir, "gt_eval_stats_blob", "json");
  std::ofstream os(out_json_filename);
  if (os)
  {
    cereal::JSONOutputArchive ar( os );
    ar << cereal::make_nvp("gt_dataset", sGTPath)
      << cereal::make_nvp("eval_dataset", sComputedSfmDataFilePath)
      << cereal::make_nvp("gt_num_poses", sfm_data_gt.GetPoses().size())
      << cereal::make_nvp("eval_num_poses", sfm_data_to_compare.GetPoses().size())
      << cereal::make_nvp("distance_residuals", vec_distance_residuals)
      << cereal::make_nvp("rotation_angular_residuals", vec_rotation_angular_residuals);
  }
  else
  {
    OPENMVG_LOG_ERROR << "Cannot export JSON stats blob output: " << out_json_filename;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
