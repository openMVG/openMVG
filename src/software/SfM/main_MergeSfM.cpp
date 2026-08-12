// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2021 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"

#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"

#include "openMVG/sfm/sfm_data_merge.hpp"

#include "openMVG/sfm/sfm_data_utils.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/system/timer.hpp"

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

enum class ESfMSceneInitializer
{
  INITIALIZE_EXISTING_POSES,
  INITIALIZE_MAX_PAIR,
  INITIALIZE_AUTO_PAIR,
  INITIALIZE_STELLAR
};

enum class ESfMEngine
{
  INCREMENTAL,
  INCREMENTALV2,
  GLOBAL
};

bool StringToEnum
(
  const std::string & str,
  ESfMEngine & sfm_engine
)
{
  const std::map<std::string, ESfMEngine> string_to_enum_mapping =
  {
    {"INCREMENTAL", ESfMEngine::INCREMENTAL},
    {"INCREMENTALV2", ESfMEngine::INCREMENTALV2},
    {"GLOBAL", ESfMEngine::GLOBAL},
  };
  const auto it  = string_to_enum_mapping.find(str);
  if (it == string_to_enum_mapping.end())
    return false;
  sfm_engine = it->second;
  return true;
}

int main(int argc, char **argv)
{
  using namespace std;
  OPENMVG_LOG_INFO
      << "\n-----------------------------------------------------------"
      << "\n Structure from Motion:"
      << "\n-----------------------------------------------------------";
  CmdLine cmd;

  // Common options:
  std::string
      filename_first_sfm_scene,
      filename_first_sfm_scene_child,
      directory_match,
      filename_match,
      directory_output,
      engine_name = "GLOBAL";

  // Bundle adjustment options:
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  std::string sExtrinsic_refinement_options = "ADJUST_ALL";
  bool b_use_motion_priors = false;
  bool preform_final_ba = false;

  // Incremental SfM options
  /*
  int triangulation_method = static_cast<int>(ETriangulationMethod::DEFAULT);
  int resection_method  = static_cast<int>(resection::SolverType::DEFAULT);
  */
  int user_camera_model = PINHOLE_CAMERA_RADIAL3;
  

  // Global SfM
  int rotation_averaging_method = int (ROTATION_AVERAGING_L2);
  int translation_averaging_method = int (TRANSLATION_AVERAGING_SOFTL1);


  // Common options
  cmd.add( make_option('i', filename_first_sfm_scene, "main_sfm_file") );
  cmd.add( make_option('c', filename_first_sfm_scene_child, "second_sfm_file") );
  cmd.add( make_option('o', directory_output, "output") );
  
  cmd.add( make_option('s', engine_name, "sfm_engine") );

  // Bundle adjustment options
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refine_intrinsic_config") );
  cmd.add( make_option('e', sExtrinsic_refinement_options, "refine_extrinsic_config") );
  cmd.add( make_switch('P', "prior_usage") );
  cmd.add( make_switch('B', "final_ba") );

  // Global SfM
  cmd.add( make_option('R', rotation_averaging_method, "rotationAveraging") );
  cmd.add( make_option('T', translation_averaging_method, "translationAveraging") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
      << "[Required]\n"
      << "[-i|--main_sfm_file] path to the parent first_sfm_scene scene\n"
      << "[-c|--second_sfm_file] path to the child first_sfm_scene scene to merge\n"
      << "[-o|--output_dir] path where the output data will be stored\n"
      << "[-s|--sfm_engine] Type of SfM Engine to use for the reconstruction\n"
      << "\t GLOBAL    : initialize globally rotation and translations\n"
      << "\n\n"
      << "[Optional parameters]\n"
      << "\n\n"
      << "[Common]\n"
      << "[-M|--match_file] path to the match file to use (i.e matches.f.txt or matches.f.bin)\n"
      << "[-f|--refine_extrinsic_config] Intrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> intrinsic parameters are held as constant\n"
      << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
      << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
      << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
      << "\t -> NOTE: options can be combined thanks to '|'\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
      <<    "\t\t-> refine the focal length & the principal point position\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
      <<    "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
      << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
      <<    "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
      << "[-e|--refine_extrinsic_config] Extrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> extrinsic parameters are held as constant\n"
      << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions) (default: false)\n"
      << "[-B|--preform_final_ba] Run a final bundle adjustment on the merge operation (default: false)\n"
      << "\n\n"
      << "[Engine specifics]\n"
      << "[GLOBAL]\n"
      << "\t[-R|--rotationAveraging]\n"
      << "\t\t 1 -> L1 minimization\n"
      << "\t\t 2 -> L2 minimization (default)\n"
      << "\t[-T|--translationAveraging]:\n"
      << "\t\t 1 -> L1 minimization\n"
      << "\t\t 2 -> L2 minimization of sum of squared Chordal distances\n"
      << "\t\t 3 -> SoftL1 minimization (default)\n"
      << "\t\t 4 -> LiGT: Linear Global Translation constraints from rotation and matches\n";

    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  b_use_motion_priors = cmd.used('P');
  preform_final_ba = cmd.used('B');

/*
  // Check validity of command line parameters:
  if ( !isValid(static_cast<ETriangulationMethod>(triangulation_method))) {
    OPENMVG_LOG_ERROR << "Invalid triangulation method";
    return EXIT_FAILURE;
  }
  */

  if ( !isValid(openMVG::cameras::EINTRINSIC(user_camera_model)) )  {
    OPENMVG_LOG_ERROR << "Invalid camera type";
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
      cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
  {
    OPENMVG_LOG_ERROR << "Invalid input for Bundle Adjustment Intrinsic parameter refinement option";
    return EXIT_FAILURE;
  }

  const sfm::Extrinsic_Parameter_Type extrinsic_refinement_options =
      sfm::StringTo_Extrinsic_Parameter_Type(sExtrinsic_refinement_options);
  if (extrinsic_refinement_options == static_cast<sfm::Extrinsic_Parameter_Type>(0) )
  {
    OPENMVG_LOG_ERROR << "Invalid input for the Bundle Adjustment Extrinsic parameter refinement option";
    return EXIT_FAILURE;
  }

  ESfMEngine sfm_engine_type;
  if (!StringToEnum(engine_name, sfm_engine_type))
  {
    OPENMVG_LOG_ERROR << "Invalid input for the SfM Engine type";
    return EXIT_FAILURE;
  }

  if (rotation_averaging_method < ROTATION_AVERAGING_L1 ||
      rotation_averaging_method > ROTATION_AVERAGING_L2 )  {
    OPENMVG_LOG_ERROR << "Rotation averaging method is invalid";
    return EXIT_FAILURE;
  }

#ifndef USE_PATENTED_LIGT
  if (translation_averaging_method == TRANSLATION_LIGT) {
    OPENMVG_LOG_ERROR << "OpenMVG was not compiled with USE_PATENTED_LIGT cmake option";
    return EXIT_FAILURE;
  }
#endif
  if (translation_averaging_method < TRANSLATION_AVERAGING_L1 ||
      translation_averaging_method > TRANSLATION_LIGT )  {
    OPENMVG_LOG_ERROR << "Translation averaging method is invalid";
    return EXIT_FAILURE;
  }

  if (directory_output.empty())  {
    OPENMVG_LOG_ERROR << "It is an invalid output directory";
    return EXIT_FAILURE;
  }

  // SfM related
  OPENMVG_LOG_INFO << "main SfM:\n"<< filename_first_sfm_scene;
  // Load input first_sfm_scene scene
  SfM_Data first_sfm_scene,second_sfm_scene;
  if (!Load(first_sfm_scene, filename_first_sfm_scene, openMVG::sfm::ESfM_Data(openMVG::sfm::ALL))) {
    OPENMVG_LOG_ERROR << "The input first_sfm_scene file \""<< filename_first_sfm_scene << "\" cannot be read.";
    return EXIT_FAILURE;
  }
  OPENMVG_LOG_INFO << "second SfM:\n"<< filename_first_sfm_scene_child;

  if (!Load(second_sfm_scene, filename_first_sfm_scene_child, openMVG::sfm::ESfM_Data(openMVG::sfm::ALL))) {
    OPENMVG_LOG_ERROR << "The input first_sfm_scene file \""<< filename_first_sfm_scene_child << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(directory_output))
  {
    if (!stlplus::folder_create(directory_output))
    {
      OPENMVG_LOG_ERROR << "Cannot create the output directory";
      return EXIT_FAILURE;
    }
  }

  //---------------------------------------
  // GLOBAL reconstruction merging operation
  //---------------------------------------

  IndexT start_views = first_sfm_scene.views.size();
  IndexT start_poses = first_sfm_scene.poses.size();
  IndexT start_tracks = first_sfm_scene.structure.size();

  std::cout << "#Main SfM Views:" << start_views  << " #Secondary SfM Views:" << second_sfm_scene.views.size() << std::endl;

  openMVG::system::Timer timer;

  // Merge every camera intrinsic present in either scene (not just camera 0): ids shared
  // by both scenes are assumed to be the same physical camera and the better of the two
  // is kept, ids only present in the second scene are brought in under a fresh id.
  const std::unordered_map<IndexT, IndexT> intrinsic_id_remap = mergeIntrinsics(first_sfm_scene, second_sfm_scene);

  std::cout << "Cameras merged: " << intrinsic_id_remap.size() << std::endl;

  // this contains the name of the file, whether it's new for the first scene and the indexes it occurs in each SfM scene
  OverlapMap sfm_filenames_indexes;

  // we'll read in the filenames and id's within the SfM scenes
  bool has_overlap = getOverlappingImages(first_sfm_scene, second_sfm_scene, sfm_filenames_indexes);
  // first is do a simple overlap test
  if(has_overlap){
    std::cout << "Do datasets have overlap: true" << std::endl;
  }else{
    std::cout << "Do datasets have overlap: false" << std::endl;
  }
  // next we'll see if it's enough to register the two scenes
  IndexT overlap_amount = getNumberOverlapping(sfm_filenames_indexes);

  OPENMVG_LOG_INFO  << "Overlap amount: " << overlap_amount;

  if(overlap_amount < 6){
    OPENMVG_LOG_ERROR << "Too few images to merge correctly, only " << overlap_amount << " overlapping images found";
    return EXIT_FAILURE;
  }

  // this could be used to print a report of the overlap information
  //std::string overlap_printout = printOverlapInformation(sfm_filenames_indexes);

  //OPENMVG_LOG_INFO << overlap_printout;

  std::vector<openMVG::Vec3> first_vecs,second_vecs,result_vecs;

  bool has_alignableVecs = getVecs2Align(first_sfm_scene, second_sfm_scene, &first_vecs, &second_vecs, sfm_filenames_indexes);
  
  OPENMVG_LOG_INFO << first_vecs.size() << " " << second_vecs.size();

  if(!has_alignableVecs){
    OPENMVG_LOG_ERROR << "Not enough poses in both SfM scenes for images that exist in both scenes";
    OPENMVG_LOG_ERROR << "maybe try a higher overlap to increase likelyhood of success";
    return EXIT_FAILURE;
  }

  Mat3 R;
  Vec3 T;
  double S;

  OPENMVG_LOG_INFO << "Calculating similarity transform:";

  bool p_res = computeSimilarity(first_vecs, second_vecs, result_vecs, &S, &R, &T);

  if(!p_res){
    OPENMVG_LOG_ERROR << " SRT transform failed to get a stable lock ";
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO << "Scale: " << S;
  OPENMVG_LOG_INFO << "Rotation:" << R;
  OPENMVG_LOG_INFO << "Translation:" << T;

  OPENMVG_LOG_INFO << "SRT calculated";

  OPENMVG_LOG_INFO << "Beginning the SRT transforms of views in scene two";

  bool merge_success = mergeSfMScenes(first_sfm_scene, second_sfm_scene, S, R, T, sfm_filenames_indexes, intrinsic_id_remap);

  if(!merge_success){
    OPENMVG_LOG_ERROR << " Merging the scenes failed";
    return EXIT_FAILURE;
  }

  // fix a new directory for the accumulation of all the images
  //std::string root_directory = directory_output.substr(0, directory_output.find_last_of("/\\"));
  //root_directory = root_directory.substr(0, root_directory.find_last_of("/\\"));

  //first_sfm_scene.s_root_path = stlplus::create_filespec(root_directory, "Originals");

  // group any intrinsics that ended up identical after the merge
  GroupSharedIntrinsics(first_sfm_scene);

  if(preform_final_ba){
    OPENMVG_LOG_INFO << "Starting final BA";

    Bundle_Adjustment_Ceres::BA_Ceres_options options;
    if ( first_sfm_scene.GetPoses().size() > 100 &&
      (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
      ){
      // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
      options.preconditioner_type_ = ceres::JACOBI;
      options.linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
    else{
      options.linear_solver_type_ = ceres::DENSE_SCHUR;
    }

    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);

    Optimize_Options ba_refine1_options(
      Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
      Extrinsic_Parameter_Type::ADJUST_ALL,// rotations are constant ADJUST_TRANSLATION
      Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure ADJUST_ALL
      Control_Point_Parameter(),
      b_use_motion_priors // Use motion priors
    );

    Optimize_Options ba_refine2_options(
      Intrinsic_Parameter_Type::NONE,
      Extrinsic_Parameter_Type::NONE,
      Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
      Control_Point_Parameter(),
      b_use_motion_priors // Use motion priors
    );

    // final polish pass: now that the scene has been cleaned of outliers, also let
    // intrinsics adjust so systematic focal/distortion error picked up by the merge
    // (e.g. two scenes shot with slightly different settings) gets corrected
    Optimize_Options ba_refine3_options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL,
      Control_Point_Parameter(),
      b_use_motion_priors // Use motion priors
    );

    // note : parameters copied from sequential sfm
    const double requiredPixelResidualError = 4.0;
    const double angle_error = 2.0;
    const size_t outlierNumberThreshold = 100;


    if(engine_name=="GLOBAL"){
      // do the initial adjustment with no changes to intrinsic to remove excess noise
      bundle_adjustment_obj.Adjust(first_sfm_scene,ba_refine1_options);

      const size_t pointcount_initial = first_sfm_scene.structure.size();
      RemoveOutliers_PixelResidualError(first_sfm_scene, requiredPixelResidualError);
      const size_t pointcount_pixelresidual_filter = first_sfm_scene.structure.size();
      RemoveOutliers_AngleError(first_sfm_scene, angle_error);
      const size_t pointcount_angular_filter = first_sfm_scene.structure.size();

      // Check that poses & intrinsic cover some measures (after outlier removal)
      const IndexT minPointPerPose = 6; // 6 min , 12
      const IndexT minTrackLength = 2; // 2 min , 3
      if (eraseUnstablePosesAndObservations(first_sfm_scene, minPointPerPose, minTrackLength))
      {
        // TODO: must ensure that track graph is producing a single connected component
        const size_t pointcount_cleaning = first_sfm_scene.structure.size();
        OPENMVG_LOG_INFO << "Point_cloud cleaning:\n"
          << "\t #3DPoints: " << pointcount_initial << " -> " << pointcount_pixelresidual_filter
          << " -> " << pointcount_angular_filter << " -> " << pointcount_cleaning << "\n";
      }

      // refine structure again now that outliers/unstable poses have been removed
      bundle_adjustment_obj.Adjust(first_sfm_scene,ba_refine2_options);

      // final pass with intrinsics unlocked
      bundle_adjustment_obj.Adjust(first_sfm_scene,ba_refine3_options);
    }
    else if(engine_name=="STELLAR"){
      OPENMVG_LOG_WARNING << "INCREMENTALV2 not implemented yet";
    }
    else if(engine_name=="INCREMENTALV2"){
      OPENMVG_LOG_WARNING << "INCREMENTALV2 not implemented yet";
    }
    else{
      // sequential method for pose recitifcation
      do{
        bundle_adjustment_obj.Adjust(first_sfm_scene,ba_refine1_options);
      }
      while (badTrackRejector(requiredPixelResidualError, outlierNumberThreshold, first_sfm_scene));
      eraseUnstablePosesAndObservations(first_sfm_scene);
    }
  }
  
  OPENMVG_LOG_INFO << "...Export first_sfm_scene to disk.";

  Generate_SfM_Report(first_sfm_scene,
    stlplus::create_filespec(directory_output, "SfMReconstruction_Report.html"));

  Save(first_sfm_scene,
    stlplus::create_filespec(directory_output, "first_sfm_scene", ".bin"),
    ESfM_Data(ALL));

  Save(first_sfm_scene,
     stlplus::create_filespec(directory_output, "cloud_and_poses", "ply"),
     ESfM_Data(ALL));

  //return EXIT_SUCCESS;
  //return EXIT_FAILURE
  return EXIT_SUCCESS;
}
