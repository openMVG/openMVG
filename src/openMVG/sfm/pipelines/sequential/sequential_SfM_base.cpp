// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON, Ricardo Fabbri and Gabriel Andrade.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <ceres/types.h>
#include <functional>
#include <iostream>
#include <utility>


#include "third_party/histogram/histogram_raw.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/sfm/base/sfm_data.hpp"
#include "openMVG/sfm/base/sfm_data_BA.hpp"
#include "openMVG/sfm/base/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/base/sfm_data_filters.hpp"
#include "openMVG/sfm/base/sfm_data_io.hpp"
#include "openMVG/sfm/base/SfM_Localizer.hpp"
#include "openMVG/sfm/base/sfm_features_provider.hpp"
#include "openMVG/sfm/base/sfm_matches_provider.hpp"
#include "openMVG/sfm/base/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/pipelines/sequential/sequential_SfM_base.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/system/loggerprogress.hpp"

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;
using namespace histogramming;

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

SequentialSfMReconstructionEngineBase::SequentialSfMReconstructionEngineBase(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory),
    sLogging_file_(sloggingFile),
    initial_pair_(0,0),
    initial_triplet_(0,0,0),
    cam_type_(EINTRINSIC(PINHOLE_CAMERA_RADIAL3))
{
  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("SequentialReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("SequentialSfMReconstructionEngineBase")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }
  // Init remaining image list
  for (Views::const_iterator itV = sfm_data.GetViews().begin();
    itV != sfm_data.GetViews().end(); ++itV)
  {
    set_remaining_view_id_.insert(itV->second->id_view);
  }
  Histogram<double> h;
}

SequentialSfMReconstructionEngineBase::~SequentialSfMReconstructionEngineBase()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_);
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void SequentialSfMReconstructionEngineBase::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void SequentialSfMReconstructionEngineBase::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}


bool SequentialSfMReconstructionEngineBase::InitLandmarkTracks() 
{
  // Compute tracks from matches
  tracks::TracksBuilder tracksBuilder;

  {
    // List of features matches for each pair of images
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
    OPENMVG_LOG_INFO << "Track building";

    tracksBuilder.Build(map_Matches);
    OPENMVG_LOG_INFO << "Track filtering";
    tracksBuilder.Filter();
    OPENMVG_LOG_INFO << "Track export to internal struct";
    //-- Build tracks with STL compliant type :
    tracksBuilder.ExportToSTL(map_tracks_);

    {
      std::ostringstream osTrack;
      //-- Display stats :
      //    - number of images
      //    - number of tracks
      std::set<uint32_t> set_imagesId;
      tracks::TracksUtilsMap::ImageIdInTracks(map_tracks_, set_imagesId);
      osTrack << "\n------------------\n"
        << "-- Tracks Stats --" << "\n"
        << " Number of tracks: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<uint32_t>(osTrack, ", "));
      osTrack << "\n------------------\n";

      std::map<uint32_t, uint32_t> map_Occurrence_TrackLength;
      tracks::TracksUtilsMap::TracksLength(map_tracks_, map_Occurrence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (const auto & it : map_Occurrence_TrackLength)  {
        osTrack << "\t" << it.first << "\t" << it.second << "\n";
      }
      OPENMVG_LOG_INFO << osTrack.str();
    }
  }
  // Initialize the shared track visibility helper
  shared_track_visibility_helper_.reset(new openMVG::tracks::SharedTrackVisibilityHelper(map_tracks_));
  return map_tracks_.size() > 0;
}

bool SequentialSfMReconstructionEngineBase::AutomaticInitialPairChoice(Pair & initial_pair) const
{
  // select a pair that have the largest baseline (mean angle between its bearing vectors).

  const unsigned iMin_inliers_count = 100;
  const float fRequired_min_angle = 3.0f;
  const float fLimit_max_angle = 60.0f; // More than 60 degree, we cannot rely on matches for initial pair seeding

  // List Views that support valid intrinsic (view that could be used for Essential matrix computation)
  std::set<IndexT> valid_views;
  for (Views::const_iterator it = sfm_data_.GetViews().begin();
    it != sfm_data_.GetViews().end(); ++it)
  {
    const View * v = it->second.get();
    if (sfm_data_.GetIntrinsics().count(v->id_intrinsic))
      valid_views.insert(v->id_view);
  }

  if (valid_views.size() < 2)
  {
    return false; // There is not view that support valid intrinsic data
  }

  std::vector<std::pair<double, Pair>> scoring_per_pair;

  // Compute the relative pose & the 'baseline score'
  system::LoggerProgress my_progress_bar( matches_provider_->pairWise_matches_.size(),
    "Selection of an initial pair");
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (const std::pair<Pair, IndMatches> & match_pair : matches_provider_->pairWise_matches_)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      ++my_progress_bar;

      const Pair current_pair = match_pair.first;

      const uint32_t I = std::min(current_pair.first, current_pair.second);
      const uint32_t J = std::max(current_pair.first, current_pair.second);
      if (valid_views.count(I) && valid_views.count(J))
      {
        const View
          * view_I = sfm_data_.GetViews().at(I).get(),
          * view_J = sfm_data_.GetViews().at(J).get();
        const Intrinsics::const_iterator
          iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
          iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

        const auto
          cam_I = iterIntrinsic_I->second.get(),
          cam_J = iterIntrinsic_J->second.get();
        if (cam_I && cam_J)
        {
          openMVG::tracks::STLMAPTracks map_tracksCommon;
          shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

          // Copy points correspondences to arrays for relative pose estimation
          const size_t n = map_tracksCommon.size();
          Mat xI(2,n), xJ(2,n);
          size_t cptIndex = 0;
          for (const auto & track_iter : map_tracksCommon)
          {
            auto iter = track_iter.second.cbegin();
            const uint32_t i = iter->second;
            const uint32_t j = (++iter)->second;

            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
            xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
            ++cptIndex;
          }

          // Robust estimation of the relative pose
          RelativePose_Info relativePose_info;
          relativePose_info.initial_residual_tolerance = Square(4.0);

          if (robustRelativePose(
                cam_I, cam_J,
                xI, xJ, relativePose_info,
                {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
                256)
              && relativePose_info.vec_inliers.size() > iMin_inliers_count)
          {
            // Triangulate inliers & compute angle between bearing vectors
            std::vector<float> vec_angles;
            vec_angles.reserve(relativePose_info.vec_inliers.size());
            const Pose3 pose_I = Pose3(Mat3::Identity(), Vec3::Zero());
            const Pose3 pose_J = relativePose_info.relativePose;
            for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
            {
              openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
              std::advance(iterT, inlier_idx);
              tracks::submapTrack::const_iterator iter = iterT->second.begin();
              const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
              const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
              vec_angles.push_back(AngleBetweenRay(pose_I, cam_I, pose_J, cam_J,
                cam_I->get_ud_pixel(featI), cam_J->get_ud_pixel(featJ)));
            }
            // Compute the median triangulation angle
            const unsigned median_index = vec_angles.size() / 2;
            std::nth_element(
              vec_angles.begin(),
              vec_angles.begin() + median_index,
              vec_angles.end());
            const float scoring_angle = vec_angles[median_index];
            // Store the pair iff the pair is in the asked angle range [fRequired_min_angle;fLimit_max_angle]
            if (scoring_angle > fRequired_min_angle &&
                scoring_angle < fLimit_max_angle)
            {
  #ifdef OPENMVG_USE_OPENMP
              #pragma omp critical
  #endif
              scoring_per_pair.emplace_back(scoring_angle, current_pair);
            }
          }
        }
      }
    } // omp section
  }
  std::sort(scoring_per_pair.begin(), scoring_per_pair.end());
  // Since scoring is ordered in increasing order, reverse the order
  std::reverse(scoring_per_pair.begin(), scoring_per_pair.end());
  if (!scoring_per_pair.empty())
  {
    initial_pair = scoring_per_pair.begin()->second;
    return true;
  }
  return false;
}

// Sketching automatic triplet generation!
// TODO XXX: MAKE IT WORK AND INTEGRATE INTO MAIN!!!
//
//bool SequentialSfMReconstructionEngineBase::AutomaticInitialTripletChoice(Triplet & initial_triplet) const
//{
//  // select a triplet that have the largest baseline (mean angle between its bearing vectors).
//  // if two of the views are close enough, discard it
//
//  const unsigned iMin_inliers_count = 100;
//  const float fRequired_min_angle = 3.0f;
//  const float fLimit_max_angle = 60.0f; // More than 60 degree, we cannot rely on matches for initial pair seeding
//
//  // List Views that support valid intrinsic (view that could be used for Essential matrix computation)
//  std::set<IndexT> valid_views;
//  for (Views::const_iterator it = sfm_data_.GetViews().begin();
//    it != sfm_data_.GetViews().end(); ++it)
//  {
//    const View * v = it->second.get();
//    if (sfm_data_.GetIntrinsics().count(v->id_intrinsic))
//      valid_views.insert(v->id_view);
//  }
//
//  if (valid_views.size() < 3)
//  {
//    return false; // There is not view that support valid intrinsic data
//  }
//
//  std::vector<std::pair<double, Pair>> scoring_per_pair;
//
//  // Compute the relative pose & the 'baseline score'
//  system::LoggerProgress my_progress_bar( matches_provider_->pairWise_matches_.size(),
//    "Selection of an initial pair");
//#ifdef OPENMVG_USE_OPENMP
//  #pragma omp parallel
//#endif
//  for (const std::pair<Pair, IndMatches> & match_pair : matches_provider_->pairWise_matches_)
//  {
//#ifdef OPENMVG_USE_OPENMP
//  #pragma omp single nowait
//#endif
//    {
//      ++my_progress_bar;
//
//      const Pair current_pair = match_pair.first;
//
//      const uint32_t I = std::min(current_pair.first, current_pair.second);
//      const uint32_t J = std::max(current_pair.first, current_pair.second);
//      if (valid_views.count(I) && valid_views.count(J))
//      {
//        const View
//          * view_I = sfm_data_.GetViews().at(I).get(),
//          * view_J = sfm_data_.GetViews().at(J).get();
//        const Intrinsics::const_iterator
//          iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
//          iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);
//
//        const auto
//          cam_I = iterIntrinsic_I->second.get(),
//          cam_J = iterIntrinsic_J->second.get();
//        if (cam_I && cam_J)
//        {
//          openMVG::tracks::STLMAPTracks map_tracksCommon;
//          shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);
//
//          // Copy points correspondences to arrays for relative pose estimation
//          const size_t n = map_tracksCommon.size();
//          Mat xI(2,n), xJ(2,n);
//          size_t cptIndex = 0;
//          for (const auto & track_iter : map_tracksCommon)
//          {
//            auto iter = track_iter.second.cbegin();
//            const uint32_t i = iter->second;
//            const uint32_t j = (++iter)->second;
//
//            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
//            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
//            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
//            xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
//            ++cptIndex;
//          }
//
//          // Robust estimation of the relative pose
//          RelativePose_Info relativePose_info;
//          relativePose_info.initial_residual_tolerance = Square(4.0);
//
//          if (robustRelativePose(
//                cam_I, cam_J,
//                xI, xJ, relativePose_info,
//                {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
//                256)
//              && relativePose_info.vec_inliers.size() > iMin_inliers_count)
//          {
//            // Triangulate inliers & compute angle between bearing vectors
//            std::vector<float> vec_angles;
//            vec_angles.reserve(relativePose_info.vec_inliers.size());
//            const Pose3 pose_I = Pose3(Mat3::Identity(), Vec3::Zero());
//            const Pose3 pose_J = relativePose_info.relativePose;
//            for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
//            {
//              openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
//              std::advance(iterT, inlier_idx);
//              tracks::submapTrack::const_iterator iter = iterT->second.begin();
//              const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
//              const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
//              vec_angles.push_back(AngleBetweenRay(pose_I, cam_I, pose_J, cam_J,
//                cam_I->get_ud_pixel(featI), cam_J->get_ud_pixel(featJ)));
//            }
//            // Compute the median triangulation angle
//            const unsigned median_index = vec_angles.size() / 2;
//            std::nth_element(
//              vec_angles.begin(),
//              vec_angles.begin() + median_index,
//              vec_angles.end());
//            const float scoring_angle = vec_angles[median_index];
//            // Store the pair iff the pair is in the asked angle range [fRequired_min_angle;fLimit_max_angle]
//            if (scoring_angle > fRequired_min_angle &&
//                scoring_angle < fLimit_max_angle)
//            {
//  #ifdef OPENMVG_USE_OPENMP
//              #pragma omp critical
//  #endif
//              scoring_per_pair.emplace_back(scoring_angle, current_pair);
//            }
//          }
//        }
//      }
//    } // omp section
//  }
//  std::sort(scoring_per_pair.begin(), scoring_per_pair.end());
//  // Since scoring is ordered in increasing order, reverse the order
//  std::reverse(scoring_per_pair.begin(), scoring_per_pair.end());
//  if (!scoring_per_pair.empty())
//  {
//    initial_pair = scoring_per_pair.begin()->second;
//    return true;
//  }
//  return false;
//}

/// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
//
// Output
//  - sfm_data_
//  - map_ACThreshold_: if ACRansac used
//  - set_remaining_view_id_: remaining views to reconstruct
bool SequentialSfMReconstructionEngineBase::MakeInitialPair3D(const Pair & current_pair)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  const uint32_t
    I = std::min(current_pair.first, current_pair.second),
    J = std::max(current_pair.first, current_pair.second);

  if (sfm_data_.GetViews().count(I) == 0)
  {
    OPENMVG_LOG_ERROR << "Cannot find the view corresponding to the view id: " << I;
    return false;
  }
  if (sfm_data_.GetViews().count(J) == 0)
  {
    OPENMVG_LOG_ERROR << "Cannot find the view corresponding to the view id: " << J;
    return false;
  }
  // a. Assert we have valid cameras
  const View
    * view_I = sfm_data_.GetViews().at(I).get(),
    * view_J = sfm_data_.GetViews().at(J).get();
  const Intrinsics::const_iterator
    iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic),
    iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

  if (iterIntrinsic_I == sfm_data_.GetIntrinsics().end() ||
      iterIntrinsic_J == sfm_data_.GetIntrinsics().end() )
  {
    OPENMVG_LOG_ERROR << "Views with valid intrinsic data are required.";
    return false;
  }

  const auto
    * cam_I = iterIntrinsic_I->second.get(),
    * cam_J = iterIntrinsic_J->second.get();
  if (!cam_I || !cam_J)
  {
    OPENMVG_LOG_ERROR << "Cannot get back the camera intrinsic model for the pair.";
    return false;
  }

  OPENMVG_LOG_INFO << "Putative starting pair info:"
    << "\nindex:(" << I << "," << J << ")"
    << "\nview basename:("
    << stlplus::basename_part(view_I->s_Img_path) << ","
    << stlplus::basename_part(view_J->s_Img_path) << ")";

  // b. Get common features between the two view
  // use the track to have a more dense match correspondence set
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

  //-- Copy point to arrays
  const size_t n = map_tracksCommon.size();
  Mat xI(2,n), xJ(2,n);
  uint32_t cptIndex = 0;
  for (const auto & track_iter : map_tracksCommon)
  {
    auto iter = track_iter.second.cbegin();
    const uint32_t
      i = iter->second,
      j = (++iter)->second;

    Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
    xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
    feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
    xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
    ++cptIndex;
  }

  // c. Robust estimation of the relative pose
  RelativePose_Info relativePose_info;

  const std::pair<size_t, size_t>
    imageSize_I(cam_I->w(), cam_I->h()),
    imageSize_J(cam_J->w(), cam_J->h());

  if (!robustRelativePose(
    cam_I, cam_J, xI, xJ, relativePose_info, imageSize_I, imageSize_J, 4096))
  {
    OPENMVG_LOG_ERROR << " /!\\ Robust estimation failed to compute E for this pair: "
      << "{"<< current_pair.first << "," << current_pair.second << "}";
    return false;
  }
  OPENMVG_LOG_INFO << "Relative pose a-contrario upper_bound residual is: "
    << relativePose_info.found_residual_precision;
  // Bound min precision at 1 pix.
  relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

  const bool bRefine_using_BA = true;
  if (bRefine_using_BA) {
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
    Landmarks & landmarks = tiny_scene.structure;
    for (const auto & track_iterator : map_tracksCommon)
    {
      // Get corresponding points
      auto iter = track_iterator.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second;

      const Vec2
        x1 = features_provider_->feats_per_view[I][i].coords().cast<double>(),
        x2 = features_provider_->feats_per_view[J][j].coords().cast<double>();

      Vec3 X;
      if (Triangulate2View(
            Pose_I.rotation(),
            Pose_I.translation(),
            (*cam_I)(cam_I->get_ud_pixel(x1)),
            Pose_J.rotation(),
            Pose_J.translation(),
            (*cam_J)(cam_J->get_ud_pixel(x2)),
            X,
            triangulation_method_))
      {
        Observations obs;
        obs[view_I->id_view] = Observation(x1, i);
        obs[view_J->id_view] = Observation(x2, j);
        landmarks[track_iterator.first].obs = obs;
        landmarks[track_iterator.first].X = X;
      }
    }
    Save(tiny_scene, stlplus::create_filespec(sOut_directory_, "initialPair.ply"), ESfM_Data(ALL));

    // - refine only Structure and Rotations & translations (keep intrinsic constant)
    Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.Adjust(tiny_scene,
        Optimize_Options
        (
          Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL) // Adjust structure
        )
      )
    {
      return false;
    }

    // Save computed data
    const Pose3 pose_I = sfm_data_.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
    const Pose3 pose_J = sfm_data_.poses[view_J->id_pose] = tiny_scene.poses[view_J->id_pose];
    map_ACThreshold_.insert({I, relativePose_info.found_residual_precision});
    map_ACThreshold_.insert({J, relativePose_info.found_residual_precision});
    set_remaining_view_id_.erase(view_I->id_view);
    set_remaining_view_id_.erase(view_J->id_view);

    // List inliers and save them
    for (const auto & landmark_entry : tiny_scene.GetLandmarks())
    {
      const IndexT trackId = landmark_entry.first;
      const Landmark & landmark = landmark_entry.second;
      const Observations & obs = landmark.obs;
      Observations::const_iterator
        iterObs_xI = obs.find(view_I->id_view),
        iterObs_xJ = obs.find(view_J->id_view);

      const Observation & ob_xI = iterObs_xI->second;
      const Observation & ob_xJ = iterObs_xJ->second;
      const Vec2
        ob_xI_ud = cam_I->get_ud_pixel(ob_xI.x),
        ob_xJ_ud = cam_J->get_ud_pixel(ob_xJ.x);

      const double angle = AngleBetweenRay(
        pose_I, cam_I, pose_J, cam_J, ob_xI_ud, ob_xJ_ud);
      const Vec2 residual_I = cam_I->residual(pose_I(landmark.X), ob_xI.x);
      const Vec2 residual_J = cam_J->residual(pose_J(landmark.X), ob_xJ.x);
      if (angle > 2.0 &&
          CheiralityTest((*cam_I)(ob_xI_ud), pose_I,
                         (*cam_J)(ob_xJ_ud), pose_J,
                         landmark.X) &&
          residual_I.norm() < relativePose_info.found_residual_precision &&
          residual_J.norm() < relativePose_info.found_residual_precision)
      {
        sfm_data_.structure[trackId] = landmarks[trackId];
      }
    }
    // Save outlier residual information
    Histogram<double> histoResiduals;
    OPENMVG_LOG_INFO
      << "\n=========================\n"
      << " MSE Residual InitialPair Inlier:\n";
    ComputeResidualsHistogram(&histoResiduals);
    std::cout << "passed Histogram\n";
    if (!sLogging_file_.empty())
    {
      using namespace htmlDocument;
      html_doc_stream_->pushInfo(htmlMarkup("h1","Essential Matrix."));
      std::ostringstream os;
      os
        << "-------------------------------" << "<br>"
        << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
        << view_I->s_Img_path << ","
        << view_J->s_Img_path << "<br>"
        << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
        << "-- Resection status: " << "OK" << "<br>"
        << "-- Nb points used for robust Essential matrix estimation: "
        << xI.cols() << "<br>"
        << "-- Nb points validated by robust estimation: "
        << sfm_data_.structure.size() << "<br>"
        << "-- % points validated: "
        << sfm_data_.structure.size()/static_cast<float>(xI.cols())
        << "<br>"
        << "-------------------------------" << "<br>";
      html_doc_stream_->pushInfo(os.str());

      html_doc_stream_->pushInfo(htmlMarkup("h2",
        "Residual of the robust estimation (Initial triangulation). Thresholded at: "
        + toString(relativePose_info.found_residual_precision)));

      html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of residuals"));

      const std::vector<double> xBin = histoResiduals.GetXbinsValue();
      const auto range = autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

      htmlDocument::JSXGraphWrapper jsxGraph;
      jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
      jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
      jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
        relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
      jsxGraph.UnsuspendUpdate();
      jsxGraph.setViewport(range);
      jsxGraph.close();
      html_doc_stream_->pushInfo(jsxGraph.toStr());

      html_doc_stream_->pushInfo("<hr>");

      std::ofstream htmlFileStream( std::string(stlplus::folder_append_separator(sOut_directory_) +
        "Reconstruction_Report.html"));
      htmlFileStream << html_doc_stream_->getDoc();
    }
  }
  return !sfm_data_.structure.empty();
}

double SequentialSfMReconstructionEngineBase::ComputeResidualsHistogram(Histogram<double> * histo) const
{
  // Collect residuals for each observation
  std::vector<float> vec_residuals;
  vec_residuals.reserve(sfm_data_.structure.size());
  // OPENMVG_LOG_INFO << "3D point info --------"; 
  for (const auto & landmark_entry : sfm_data_.GetLandmarks())
  {
    // std::cerr << "\tX" << landmark_entry.second.X.transpose();
    const Observations & obs = landmark_entry.second.obs;
    for (const auto & observation : obs)
    {
      const View * view = sfm_data_.GetViews().find(observation.first)->second.get();
      const Pose3 pose = sfm_data_.GetPoseOrDie(view);
      // OPENMVG_LOG_INFO << "Pose rc " << pose.rotation() << "\n" << pose.center().transpose();
      const auto intrinsic = sfm_data_.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose(landmark_entry.second.X), observation.second.x);
      // std::cerr << " residual " << residual.transpose();
      vec_residuals.emplace_back( std::abs(residual(0)) );
      vec_residuals.emplace_back( std::abs(residual(1)) );
    }
    std::cerr << std::endl;
  }
  // Display statistics
  if (vec_residuals.size() > 1)
  {
    float dMin, dMax, dMean, dMedian;
    minMaxMeanMedian<float>(vec_residuals.cbegin(), vec_residuals.cend(),
                            dMin, dMax, dMean, dMedian);
    if (histo)  {
      *histo = Histogram<double>(dMin, dMax, 10);
      histo->Add(vec_residuals.cbegin(), vec_residuals.cend());
    }

    OPENMVG_LOG_INFO
      << "\nSequentialSfMReconstructionEngineBase::ComputeResidualsMSE."
      << "\n\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size()
      << "\n\t-- Residual min:\t" << dMin
      << "\n\t-- Residual median:\t" << dMedian
      << "\n\t-- Residual max:\t "  << dMax
      << "\n\t-- Residual mean:\t " << dMean;

    return dMean;
  }
  return -1.0;
}

/// Functor to sort a vector of pair given the pair's second value
template<class T1, class T2, class Pred = std::less<T2>>
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left,
                    const std::pair<T1,T2>&right)
  {
    Pred p;
    return p(left.second, right.second);
  }
};

/**
 * @brief Discard tracks with too large residual error
 *
 * Remove observation/tracks that have:
 *  - too large residual error
 *  - too small angular value
 *
 * @return True if more than 'count' outliers have been removed.
 */
bool SequentialSfMReconstructionEngineBase::badTrackRejector(double dPrecision, size_t count)
{
  const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(sfm_data_, dPrecision, 2);
  const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(sfm_data_, 2.0);

  return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
}

void SequentialSfMReconstructionEngineBase::FinalStatistics()
{
  //-- Reconstruction done.
  //-- Display some statistics
  std::ostringstream os_sfm_stats;
  os_sfm_stats << "\n-------------------------------\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
    << " from " << sfm_data_.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  os_sfm_stats << "\nHistogram of residuals:\n" << h.ToString();

  OPENMVG_LOG_INFO << os_sfm_stats.str();

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion process finished.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Structure from Motion (statistics):<br>"
      << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
      << " from " <<sfm_data_.GetViews().size() << " input images.<br>"
      << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());

    html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

    const std::vector<double> xBin = h.GetXbinsValue();
    const auto range = autoJSXGraphViewport<double>(xBin, h.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("3DtoImageResiduals",600,300);
    jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());
  }
}

/**
 * @brief Estimate images on which we can compute the resectioning safely.
 *
 * @param[out] vec_possible_indexes: list of indexes we can use for resectioning.
 * @return False if there is no possible resection.
 *
 * Sort the images by the number of features id shared with the reconstruction.
 * Select the image I that share the most of correspondences.
 * Then keep all the images that have at least:
 *  0.75 * #correspondences(I) common correspondences to the reconstruction.
 */
bool SequentialSfMReconstructionEngineBase::FindImagesWithPossibleResection(
  std::vector<uint32_t> & vec_possible_indexes)
{
  // Threshold used to select the best images
  static const float dThresholdGroup = 0.75f;

  vec_possible_indexes.clear();

  if (set_remaining_view_id_.empty() || sfm_data_.GetLandmarks().empty())
    return false;

  // Collect tracksIds
  std::set<uint32_t> reconstructed_trackId;
  std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  Pair_Vec vec_putative; // ImageId, NbPutativeCommonPoint
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (std::set<uint32_t>::const_iterator iter = set_remaining_view_id_.begin();
        iter != set_remaining_view_id_.end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      const uint32_t viewId = *iter;

      // Compute 2D - 3D possible content
      openMVG::tracks::STLMAPTracks map_tracksCommon;
      shared_track_visibility_helper_->GetTracksInImages({viewId}, map_tracksCommon);

      if (!map_tracksCommon.empty())
      {
        std::set<uint32_t> set_tracksIds;
        tracks::TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

        // Count the common possible putative point
        //  with the already 3D reconstructed trackId
        std::vector<uint32_t> vec_trackIdForResection;
        std::set_intersection(set_tracksIds.cbegin(), set_tracksIds.cend(),
          reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
          std::back_inserter(vec_trackIdForResection));

#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          vec_putative.emplace_back(viewId, vec_trackIdForResection.size());
        }
      }
    }
  }

  // Sort by the number of matches to the 3D scene.
  std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<uint32_t, uint32_t, std::greater<uint32_t>>());

  // If the list is empty or if the list contains images with no correspdences
  // -> (no resection will be possible)
  if (vec_putative.empty() || vec_putative[0].second == 0)
  {
    // All remaining images cannot be used for pose estimation
    set_remaining_view_id_.clear();
    return false;
  }

  // Add the image view index that share the most of 2D-3D correspondences
  vec_possible_indexes.push_back(vec_putative[0].first);

  // Then, add all the image view indexes that have at least N% of the number of the matches of the best image.
  const IndexT M = vec_putative[0].second; // Number of 2D-3D correspondences
  const size_t threshold = static_cast<uint32_t>(dThresholdGroup * M);
  for (size_t i = 1; i < vec_putative.size() &&
    vec_putative[i].second > threshold; ++i)
  {
    vec_possible_indexes.push_back(vec_putative[i].first);
  }
  return true;
}

/// Bundle adjustment to refine Structure; Motion and Intrinsics
bool SequentialSfMReconstructionEngineBase::BundleAdjustment()
{
  Bundle_Adjustment_Ceres::BA_Ceres_options options;
  if ( sfm_data_.GetPoses().size() > 100 &&
      (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
      )
  // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
  {
    options.preconditioner_type_ = ceres::JACOBI;
    options.linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
  }
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  const Optimize_Options ba_refine_options
    ( ReconstructionEngine::intrinsic_refinement_options_,
      Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
      Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
      Control_Point_Parameter(),
      this->b_use_motion_prior_
    );
  return bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
}

} // namespace sfm
} // namespace openMVG
