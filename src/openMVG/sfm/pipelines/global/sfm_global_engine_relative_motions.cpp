
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/multiview/essential.hpp"

#include "openMVG/sfm/pipelines/non_central_cameras/sfm_robust_relative_pose_rig.hpp"

#include "third_party/progress/progress.hpp"

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

GlobalSfMReconstructionEngine_RelativeMotions::GlobalSfMReconstructionEngine_RelativeMotions(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory), sLogging_file_(sloggingFile)
{

  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("GlobalReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("GlobalSfMReconstructionEngine_RelativeMotions")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }

  // Set default motion Averaging methods
  eRotation_averaging_method_ = ROTATION_AVERAGING_L2;
  eTranslation_averaging_method_ = TRANSLATION_AVERAGING_L1;
}

GlobalSfMReconstructionEngine_RelativeMotions::~GlobalSfMReconstructionEngine_RelativeMotions()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetRotationAveragingMethod
(
  ERotationAveragingMethod eRotationAveragingMethod
)
{
  eRotation_averaging_method_ = eRotationAveragingMethod;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetTranslationAveragingMethod
(
  ETranslationAveragingMethod eTranslationAveragingMethod
)
{
  eTranslation_averaging_method_ = eTranslationAveragingMethod;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process() {

  //-------------------
  // Keep only the largest biedge connected pose subgraph
  //-------------------
  {
    // Build the relative poses pair list
    Pair_Set relative_pose_pairs;
    for (const auto & matches_it : matches_provider_->pairWise_matches_)
    {
      const Pair pair = matches_it.first;
      const View * v1 = sfm_data_.GetViews().at(pair.first).get();
      const View * v2 = sfm_data_.GetViews().at(pair.second).get();
      if (v1->id_pose != v2->id_pose)
      {
        relative_pose_pairs.insert(
          Pair(
            std::min(v1->id_pose, v2->id_pose),
            std::max(v1->id_pose, v2->id_pose))
          );
      }
    }

    const std::set<IndexT> set_remaining_poseIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(relative_pose_pairs);
    if(set_remaining_poseIds.empty())
    {
      std::cerr << "Invalid input pose graph for global SfM" << std::endl;
      return false;
    }
    openMVG::matching::PairWiseMatches infered_matches;
    for (const auto & matches_it : matches_provider_->pairWise_matches_)
    {
      const Pair pair = matches_it.first;
      const View * v1 = sfm_data_.GetViews().at(pair.first).get();
      const View * v2 = sfm_data_.GetViews().at(pair.second).get();
      if (set_remaining_poseIds.count(v1->id_pose) &&
          set_remaining_poseIds.count(v2->id_pose))
      {
        infered_matches.insert(matches_it);
      }
    }
    matches_provider_->pairWise_matches_.swap(infered_matches);
  }

  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Compute_Relative_Rotations(relatives_R);

  Hash_Map<IndexT, Mat3> global_rotations;
  if (!Compute_Global_Rotations(relatives_R, global_rotations))
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }

  matching::PairWiseMatches  tripletWise_matches;
  if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure(tripletWise_matches))
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << sfm_data_.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
      << "-- Track count: "  << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  return true;
}

/// Compute from relative rotations the global rotations of the camera poses
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Rotations
(
  const rotation_averaging::RelativeRotations & relatives_R,
  Hash_Map<IndexT, Mat3> & global_rotations
)
{
  if(relatives_R.empty())
    return false;
  // Log statistics about the relative rotation graph
  {
    std::set<IndexT> set_pose_ids;
    for (const auto & relative_R : relatives_R)
    {
      set_pose_ids.insert(relative_R.i);
      set_pose_ids.insert(relative_R.j);
    }

    std::cout << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << relatives_R.size() << "\n"
      << "  #global rotations: " << set_pose_ids.size() << std::endl;
  }

  // Global Rotation solver:
  const ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod =
    TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    //TRIPLET_ROTATION_INFERENCE_NONE;

  system::Timer t;
  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  const bool b_rotation_averaging = rotation_averaging_solver.Run(
    eRotation_averaging_method_, eRelativeRotationInferenceMethod,
    relatives_R, global_rotations);

  std::cout
    << "Found #global_rotations: " << global_rotations.size() << "\n"
    << "Timing: " << t.elapsed() << " seconds" << std::endl;


  if (b_rotation_averaging)
  {
    // Compute & display rotation fitting residual errors
    std::vector<float> vec_rotation_fitting_error;
    vec_rotation_fitting_error.reserve(relatives_R.size());
    for (const auto & relative_R : relatives_R)
    {
      const Mat3 & Rij = relative_R.Rij;
      const IndexT i = relative_R.i;
      const IndexT j = relative_R.j;
      if (global_rotations.count(i)==0 || global_rotations.count(j)==0)
        continue;
      const Mat3 & Ri = global_rotations[i];
      const Mat3 & Rj = global_rotations[j];
      const Mat3 eRij(Rj.transpose()*Rij*Ri);
      const double angularErrorDegree = R2D(getRotationMagnitude(eRij));
      vec_rotation_fitting_error.push_back(angularErrorDegree);
    }

    if (!vec_rotation_fitting_error.empty())
    {
      const float error_max = *max_element(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
      Histogram<float> histo(0.0f,error_max, 20);
      histo.Add(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
      std::cout
        << "\nRelative/Global degree rotations residual errors {0," << error_max<< "}:"
        << histo.ToString() << std::endl;
      {
        Histogram<float> histo(0.0f, 5.0f, 20);
        histo.Add(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
        std::cout
          << "\nRelative/Global degree rotations residual errors {0,5}:"
          << histo.ToString() << std::endl;
      }
      std::cout << "\nStatistics about global rotation evaluation:" << std::endl;
      minMaxMeanMedian<float>(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
    }

    // Log input graph to the HTML report
    if (!sLogging_file_.empty() && !sOut_directory_.empty())
    {
      // Log a relative pose graph
      {
        std::set<IndexT> set_pose_ids;
        Pair_Set relative_pose_pairs;
        for (const auto & view : sfm_data_.GetViews())
        {
          const IndexT pose_id = view.second->id_pose;
          set_pose_ids.insert(pose_id);
        }
        const std::string sGraph_name = "global_relative_rotation_pose_graph_final";
        graph::indexedGraph putativeGraph(set_pose_ids, rotation_averaging_solver.GetUsedPairs());
        graph::exportToGraphvizData(
          stlplus::create_filespec(sOut_directory_, sGraph_name),
          putativeGraph);

        using namespace htmlDocument;
        std::ostringstream os;

        os << "<br>" << sGraph_name << "<br>"
           << "<img src=\""
           << stlplus::create_filespec(sOut_directory_, sGraph_name, "svg")
           << "\" height=\"600\">\n";

        html_doc_stream_->pushInfo(os.str());
      }
    }
  }
  return b_rotation_averaging;
}

/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run(
    eTranslation_averaging_method_,
    sfm_data_,
    features_provider_,
    matches_provider_,
    global_rotations,
    tripletWise_matches);

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

/// Compute the initial structure of the scene
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Initial_Structure
(
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Build tracks from selected triplets (Union of all the validated triplet tracks (_tripletWise_matches))
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
#if defined USE_ALL_VALID_MATCHES // not used by default
    matching::PairWiseMatches pose_supported_matches;
    for (const std::pair< Pair, IndMatches > & match_info :  matches_provider_->pairWise_matches_)
    {
      const View * vI = sfm_data_.GetViews().at(match_info.first.first).get();
      const View * vJ = sfm_data_.GetViews().at(match_info.first.second).get();
      if (sfm_data_.IsPoseAndIntrinsicDefined(vI) && sfm_data_.IsPoseAndIntrinsicDefined(vJ))
      {
        pose_supported_matches.insert(match_info);
      }
    }
    tracksBuilder.Build(pose_supported_matches);
#else
    // Use triplet validated matches
    tracksBuilder.Build(tripletWise_matches);
#endif
    tracksBuilder.Filter(3);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = sfm_data_.structure;
    IndexT idx(0);
    for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
      itTracks != map_selectedTracks.end();
      ++itTracks, ++idx)
    {
      const submapTrack & track = itTracks->second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
      {
        const size_t imaIndex = it->first;
        const size_t featIndex = it->second;
        const PointFeature & pt = features_provider_->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

    std::cout << std::endl << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats:
      //    - number of images
      //    - number of tracks
      std::set<size_t> set_imagesId;
      TracksUtilsMap::ImageIdInTracks(map_selectedTracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(map_selectedTracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)  {
        osTrack << "\t" << iter->first << "\t" << iter->second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Compute 3D position of the landmark of the structure by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = sfm_data_.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(sfm_data_);

    std::cout << "\n#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(sfm_data_.GetLandmarks().size()) << std::endl;
    std::cout << std::endl << "  Triangulation took (s): " << timer.elapsed() << std::endl;

    // Export initial structure
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !sfm_data_.structure.empty();
}

// Adjust the scene (& remove outliers)
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL)
    );
  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL); // adjust scene structure

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

void GlobalSfMReconstructionEngine_RelativeMotions::Compute_Relative_Rotations
(
  rotation_averaging::RelativeRotations & vec_relatives_R
)
{
  //
  // Build the Relative pose graph from matches:
  //
  /// pairwise view relation between poseIds
  typedef std::map< Pair, Pair_Set > PoseWiseMatches;

  // List shared correspondences (pairs) between poses
  PoseWiseMatches poseWiseMatches;
  for (const auto & iterMatches : matches_provider_->pairWise_matches_)
  {
    const Pair pair = iterMatches.first;
    const View * v1 = sfm_data_.GetViews().at(pair.first).get();
    const View * v2 = sfm_data_.GetViews().at(pair.second).get();
    poseWiseMatches[Pair(v1->id_pose, v2->id_pose)].insert(pair);
  }

  C_Progress_display my_progress_bar( poseWiseMatches.size(),
      std::cout, "\n- Relative pose computation -\n" );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  // Compute the relative pose from pairwise point matches:
  for (int i = 0; i < static_cast<int>(poseWiseMatches.size()); ++i)
  {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
    {
      ++my_progress_bar;
    }

    PoseWiseMatches::const_iterator iter (poseWiseMatches.begin());
    std::advance(iter, i);
    const auto & relative_pose_iterator(*iter);
    const Pair relative_pose_pair = relative_pose_iterator.first;
    const Pair_Set & match_pairs = relative_pose_iterator.second;

    // If a pair has the same ID, discard it
    if (relative_pose_pair.first == relative_pose_pair.second)
    {
      continue;
    }

    // Compute the relative pose:
    // - non central camera
    // - central camera
    if (match_pairs.size() > 1)
    {
      // We have a non central cameras
      //
      // - Select used intrinsics
      std::set<IndexT> used_intrinsic_ids;
      for (const auto & pair_it : match_pairs)
      {
        const View * view_1 = sfm_data_.GetViews().at(pair_it.first).get();
        const View * view_2 = sfm_data_.GetViews().at(pair_it.second).get();

        // Check that valid camera are existing for the pair view
        if (sfm_data_.GetIntrinsics().count(view_1->id_intrinsic) == 0 ||
            sfm_data_.GetIntrinsics().count(view_2->id_intrinsic) == 0)
          continue;

        used_intrinsic_ids.insert(view_1->id_intrinsic);
        used_intrinsic_ids.insert(view_2->id_intrinsic);
      }
      // Setup the local poses used by the non central camera
      opengv::translations_t  rigOffsets;
      opengv::rotations_t     rigRotations;
      double                  minFocal = std::numeric_limits<double>::max();
      for (const auto & intrinsic_id_it : used_intrinsic_ids)
      {
        const cameras::IntrinsicBase * intrinsic_ptr = sfm_data_.GetIntrinsics().at(intrinsic_id_it).get();
        if (
              intrinsic_ptr != nullptr &&
              intrinsic_ptr->getType() == cameras::PINHOLE_CAMERA_SUBPOSE
           )
        {
          // retrieve camera information from shared pointer
          const cameras::Pinhole_Intrinsic_Subpose * subpose_intrinsic_ptr =
            dynamic_cast< const cameras::Pinhole_Intrinsic_Subpose * > (intrinsic_ptr);
          const geometry::Pose3 sub_pose = subpose_intrinsic_ptr->get_subpose();

          rigOffsets.emplace_back(sub_pose.center());
          rigRotations.emplace_back(sub_pose.rotation().transpose());

          minFocal = std::min( minFocal , subpose_intrinsic_ptr->focal() );
        }
      }

      //--
      // Setup bearing vector correspondences by using feature tracking
      //--
      opengv::bearingVectors_t bI, bJ;
      std::vector<int> intrinsic_index_I, intrinsic_index_J;

      for (const auto & pair_it : match_pairs)
      {
        const IndexT I = pair_it.first;
        const IndexT J = pair_it.second;

        const View * view_I = sfm_data_.GetViews().at(I).get();
        const View * view_J = sfm_data_.GetViews().at(J).get();

        // Check that valid camera are existing for the pair view
        if (used_intrinsic_ids.count(view_I->id_intrinsic) == 0 ||
            used_intrinsic_ids.count(view_J->id_intrinsic) == 0)
          continue;

        // add bearing vectors if they do not belong to the same pose
        if ( view_I->id_pose == view_J->id_pose )
          continue;

        const cameras::IntrinsicBase * intrinsic_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
        const cameras::IntrinsicBase * intrinsic_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

        for (const auto & match_it : matches_provider_->pairWise_matches_.at(pair_it))
        {
          opengv::bearingVector_t bearing_vector;
          bearing_vector << (*intrinsic_I)(intrinsic_I->get_ud_pixel(features_provider_->feats_per_view[I][match_it.i_].coords().cast<double>()));
          bI.emplace_back( bearing_vector.normalized() );
          intrinsic_index_I.emplace_back( std::distance( used_intrinsic_ids.begin(), used_intrinsic_ids.find(view_I->id_intrinsic) ) );

          bearing_vector << (*intrinsic_J)(intrinsic_J->get_ud_pixel(features_provider_->feats_per_view[J][match_it.j_].coords().cast<double>()));
          bJ.emplace_back( bearing_vector.normalized() );
          intrinsic_index_J.emplace_back( std::distance( used_intrinsic_ids.begin(), used_intrinsic_ids.find(view_J->id_intrinsic) ) );
        }
      }

      //--
      // Robust pose estimation
      //--

      //--> Estimate the best possible Rotation/Translation from correspondences
      double errorMax = std::numeric_limits<double>::max();
      //const double upper_bound_pixel_threshold = 2.5 * sqrt(2.0);
      const double maxExpectedError =
        // AC Ransac converge to better solution (in angular residual error mode) when no apriori threshold is set
        // It avoids pre-optimization of a false solution
        std::numeric_limits<double>::infinity();
        //1.0 - cos ( atan ( upper_bound_pixel_threshold / minFocal ) );

      opengv::transformation_t relative_pose;
      std::vector<size_t> vec_inliers;

      const bool isPoseUsable =
        noncentral::robust_relative_pose::non_central_cameras_robust_relative_pose(
        bI, bJ,
        intrinsic_index_I, intrinsic_index_J,
        rigOffsets, rigRotations,
        &relative_pose,
        &vec_inliers,
        &errorMax,
        maxExpectedError);

      if ( isPoseUsable )
      {
        // set output model
        geometry::Pose3 relativePose(relative_pose.block<3,3>(0,0).transpose(), relative_pose.col(3));

        // Build a tiny SfM scene with only the geometry of the relative pose
        //  for parameters refinement: 3D points & camera poses.
        SfM_Data tiny_scene;

        const IndexT indexRig1 = relative_pose_pair.first;
        const IndexT indexRig2 = relative_pose_pair.second;
        // intialize poses (shared by a group of images)
        tiny_scene.poses[indexRig1] = geometry::Pose3(Mat3::Identity(), Vec3::Zero());
        tiny_scene.poses[indexRig2] = relativePose;

        // insert views used by the relative pose pairs
        for (const auto & pairIterator : match_pairs)
        {
          // initialize camera indexes
          const IndexT I = pairIterator.first;
          const IndexT J = pairIterator.second;

          // add views
          tiny_scene.views.insert(*sfm_data_.GetViews().find(pairIterator.first));
          tiny_scene.views.insert(*sfm_data_.GetViews().find(pairIterator.second));

          // add intrinsics
          const View * view_I = sfm_data_.GetViews().at(I).get();
          const View * view_J = sfm_data_.GetViews().at(J).get();
          if (tiny_scene.intrinsics.count(view_I->id_intrinsic) == 0)
            tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
          if (tiny_scene.intrinsics.count(view_J->id_intrinsic) == 0)
            tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));
        }

        const bool b_use_tracks = true;
        if (b_use_tracks)
        {
          matching::PairWiseMatches relative_pose_matches;
          for (const auto & pair_it : match_pairs)
          {
            const IndexT I = pair_it.first;
            const IndexT J = pair_it.second;

            const View * view_I = tiny_scene.GetViews().at(I).get();
            const View * view_J = tiny_scene.GetViews().at(J).get();

            // Check that valid camera are existing for the pair view
            if (tiny_scene.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
                tiny_scene.GetIntrinsics().count(view_J->id_intrinsic) == 0)
              continue;

            // Consider this pair
            // - If they define a minimal of parallax ?
            // - if they do not belong to the same pose
            if ( view_I->id_pose == view_J->id_pose )
              continue;

            relative_pose_matches.insert(
              std::make_pair(
                pair_it,
                matches_provider_->pairWise_matches_.at( pair_it )));
          }
          // Compute tracks
          using namespace openMVG::tracks;
          TracksBuilder tracksBuilder;
          tracksBuilder.Build( relative_pose_matches );
          tracksBuilder.Filter( 2 ); // matches must be seen by at least two view/pose.
          STLMAPTracks map_tracks;   // reconstructed track (visibility per 3D point)
          tracksBuilder.ExportToSTL(map_tracks);

          // List tracks as observations
          for (const auto & tracks_it : map_tracks)
          {
            Landmark landmark;
            Observations & obs = landmark.obs;
            for (const auto & track_it : tracks_it.second)
            {
              const size_t view_index = track_it.first;
              const size_t feat_index = track_it.second;
              obs[view_index] = Observation(features_provider_->feats_per_view[view_index][feat_index].coords().cast<double>(), feat_index);
            }
            tiny_scene.structure[tracks_it.first] = std::move(landmark);
          }// end loop on tracks
        }
        else
        {
          // Add 2-view observations
          IndexT cpt = 0;
          Landmarks & structure = tiny_scene.structure;
          const std::set<size_t> set_inliers(vec_inliers.begin(), vec_inliers.end());

          for (const auto & pair_it : match_pairs)
          {
            const IndexT I = pair_it.first;
            const IndexT J = pair_it.second;

            const View * view_I = tiny_scene.GetViews().at(I).get();
            const View * view_J = tiny_scene.GetViews().at(J).get();

            // Check that valid camera are existing for the pair view
            if (tiny_scene.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
                tiny_scene.GetIntrinsics().count(view_J->id_intrinsic) == 0)
              continue;

            // add bearing vectors if they do not belong to the same pose
            if ( view_I->id_pose == view_J->id_pose )
              continue;

            if (set_inliers.count(cpt))
            {
              for (const auto & match_it : matches_provider_->pairWise_matches_.at(pair_it))
              {
                Observations & obs = structure[cpt].obs;
                PointFeature pt = features_provider_->feats_per_view.at(I)[match_it.i_];
                obs[I] = Observation(pt.coords().cast<double>(), match_it.i_);

                pt = features_provider_->feats_per_view.at(J)[match_it.j_];
                obs[J] = Observation(pt.coords().cast<double>(), match_it.j_);

                ++cpt;
              }
            }
          }
        }

        // Compute 3D position of the landmarks (triangulation of the observations)
        {
          SfM_Data_Structure_Computation_Blind structure_estimator(false);
          structure_estimator.triangulate(tiny_scene);
        }

        /*std::ostringstream os;
        os << this->sOut_directory_ << "/" << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
        Save(tiny_scene, os.str(), ESfM_Data(ALL));*/

        Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
          (Intrinsic_Parameter_Type::NONE, // -> Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL);// adjust scene structure
        if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
        {
          // --> to debug: save relative pair geometry on disk
          /*std::ostringstream os;
          os << this->sOut_directory_ << "/" << relative_pose_pair.first << "_" << relative_pose_pair.second << "_BA" << ".ply";
          Save(tiny_scene, os.str(), ESfM_Data(ALL));*/

          // Remove outliers (max_angle, residual error)
          RemoveOutliers_PixelResidualError(tiny_scene, 4.0);
          RemoveOutliers_AngleError(tiny_scene, 2.0);

          if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
          {
            const Mat3 R1 = tiny_scene.poses[relative_pose_pair.first].rotation();
            const Mat3 R2 = tiny_scene.poses[relative_pose_pair.second].rotation();
            const Vec3 t1 = tiny_scene.poses[relative_pose_pair.first].translation();
            const Vec3 t2 = tiny_scene.poses[relative_pose_pair.second].translation();
            // Compute relative motion and save it
            Mat3 Rrel;
            Vec3 trel;
            RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
            // Update found relative pose
            relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
            /*std::ostringstream os;
            os << this->sOut_directory_ << "/" << relative_pose_pair.first << "_" << relative_pose_pair.second << "_BA2" << ".ply";
            Save(tiny_scene, os.str(), ESfM_Data(ALL));*/
          }
        }
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          // Add the relative rotation to the relative 'rotation' pose graph
          using namespace openMVG::rotation_averaging;
            vec_relatives_R.emplace_back(
              relative_pose_pair.first, relative_pose_pair.second,
              relativePose.rotation(), 1);
        }
      }
    }
    else
    // Central camera relative pose estimation
    {
      const Pair pairIterator = *(match_pairs.begin());

      const IndexT I = pairIterator.first;
      const IndexT J = pairIterator.second;

      const View * view_I = sfm_data_.views[I].get();
      const View * view_J = sfm_data_.views[J].get();

      // Check that valid cameras are existing for the pair of view
      if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
        sfm_data_.GetIntrinsics().count(view_J->id_intrinsic) == 0)
        continue;

      // Setup corresponding bearing vectors
      const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
      const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

      const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(pairIterator);
      size_t nbBearing = matches.size();
      Mat x1(2, nbBearing), x2(2, nbBearing);
      nbBearing = 0;
      for (const auto & match : matches)
      {
        x1.col(nbBearing) = (*cam_I)(cam_I->get_ud_pixel(features_provider_->feats_per_view[I][match.i_].coords().cast<double>())).hnormalized();
        x2.col(nbBearing++) = (*cam_J)(cam_J->get_ud_pixel(features_provider_->feats_per_view[J][match.j_].coords().cast<double>())).hnormalized();
      }

      RelativePose_Info relativePose_info;
      // Compute max authorized error as geometric mean of camera plane tolerated residual error
      relativePose_info.initial_residual_tolerance = std::pow(
        cam_I->imagePlane_toCameraPlaneError(2.5) *
        cam_J->imagePlane_toCameraPlaneError(2.5),
        1./2.);

      // Since we use normalized features, we will use unit image size and intrinsic matrix:
      const std::pair<size_t, size_t> imageSize(1., 1.);
      const Mat3 K  = Mat3::Identity();

      if (!robustRelativePose(K, K, x1, x2, relativePose_info, imageSize, imageSize, 256))
      {
        continue;
      }
      const bool bRefine_using_BA = true;
      if (bRefine_using_BA)
      {
        // Refine the defined scene
        SfM_Data tiny_scene;
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

        // Init poses
        const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
        const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

        // Init structure
        const Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
        const Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
        Landmarks & landmarks = tiny_scene.structure;
        for (Mat::Index k = 0; k < x1.cols(); ++k) {
          const Vec2 x1_ = features_provider_->feats_per_view[I][matches[k].i_].coords().cast<double>();
          const Vec2 x2_ = features_provider_->feats_per_view[J][matches[k].j_].coords().cast<double>();
          Vec3 X;
          TriangulateDLT(P1, x1_, P2, x2_, &X);
          Observations obs;
          obs[view_I->id_view] = Observation(x1_, matches[k].i_);
          obs[view_J->id_view] = Observation(x2_, matches[k].j_);
          landmarks[k].obs = obs;
          landmarks[k].X = X;
        }
        // - refine only Structure and Rotations & translations (keep intrinsic constant)
        Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
          (Intrinsic_Parameter_Type::NONE, // -> Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL);// adjust scene structure
        if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
        {
          // --> to debug: save relative pair geometry on disk
          // std::ostringstream os;
          // os << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
          // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
          //
          const Mat3 R1 = tiny_scene.poses[view_I->id_pose].rotation();
          const Mat3 R2 = tiny_scene.poses[view_J->id_pose].rotation();
          const Vec3 t1 = tiny_scene.poses[view_I->id_pose].translation();
          const Vec3 t2 = tiny_scene.poses[view_J->id_pose].translation();
          // Compute relative motion and save it
          Mat3 Rrel;
          Vec3 trel;
          RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
          // Update found relative pose
          relativePose_info.relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
        }
      }
#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
      {
        // Add the relative rotation to the relative 'rotation' pose graph
        using namespace openMVG::rotation_averaging;
          vec_relatives_R.emplace_back(
            relative_pose_pair.first, relative_pose_pair.second,
            relativePose_info.relativePose.rotation(),
            1.f);
      }
    }
  } // for all relative pose

  // Log input graph to the HTML report
  if (!sLogging_file_.empty() && !sOut_directory_.empty())
  {
    // Log a relative view graph
    {
      std::set<IndexT> set_ViewIds;
      std::transform(sfm_data_.GetViews().begin(), sfm_data_.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(matches_provider_->pairWise_matches_));
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, "global_relative_rotation_view_graph"),
        putativeGraph);
    }

    // Log a relative pose graph
    {
      std::set<IndexT> set_pose_ids;
      Pair_Set relative_pose_pairs;
      for (const auto & relative_R : vec_relatives_R)
      {
        const Pair relative_pose_indices(relative_R.i, relative_R.j);
        relative_pose_pairs.insert(relative_pose_indices);
        set_pose_ids.insert(relative_R.i);
        set_pose_ids.insert(relative_R.j);
      }
      const std::string sGraph_name = "global_relative_rotation_pose_graph";
      graph::indexedGraph putativeGraph(set_pose_ids, relative_pose_pairs);
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, sGraph_name),
        putativeGraph);
      using namespace htmlDocument;
      std::ostringstream os;

      os << "<br>" << "global_relative_rotation_pose_graph" << "<br>"
         << "<img src=\""
         << stlplus::create_filespec(sOut_directory_, "global_relative_rotation_pose_graph", "svg")
         << "\" height=\"600\">\n";

      html_doc_stream_->pushInfo(os.str());
    }
  }
}

} // namespace sfm
} // namespace openMVG

