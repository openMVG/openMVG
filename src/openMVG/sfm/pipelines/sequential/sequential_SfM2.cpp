// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sequential/sequential_SfM2.hpp"
#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializer.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/stl/stl.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#include <array>
#include <ceres/types.h>
#include <functional>
#include <iostream>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;

SequentialSfMReconstructionEngine2::SequentialSfMReconstructionEngine2(
  SfMSceneInitializer * scene_initializer,
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory),
    scene_initializer_(scene_initializer),
    sLogging_file_(sloggingFile),
    cam_type_(EINTRINSIC(PINHOLE_CAMERA_RADIAL3))
{
  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("SequentialReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("SequentialSfMReconstructionEngine")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
    html_doc_stream_->pushInfo( "Poses count: " +
      htmlDocument::toString( sfm_data.GetPoses().size()) + "<br>");
  }
}

SequentialSfMReconstructionEngine2::~SequentialSfMReconstructionEngine2()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void SequentialSfMReconstructionEngine2::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void SequentialSfMReconstructionEngine2::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

bool SequentialSfMReconstructionEngine2::Process() {

  //-------------------
  //-- Incremental reconstruction
  //-------------------
  //- 1. Init the reconstruction with a seed
  //- 2. While we can localize some cameras in the reconstruction
  //     a. Triangulation
  //     b. Bundle Adjustment and cleaning
  //- 3. Final bundle Adjustment
  //-------------------

  //--
  //- 1. Init the reconstruction with a Seed
  //--
  {
    if (!scene_initializer_ || !scene_initializer_->Process())
    {
      std::cout << "Initialization status: Failed" << std::endl;
      return false;
    }
    else
    {
      std::cout << "Initialization status : Success" << std::endl;
      sfm_data_.poses = scene_initializer_->Get_sfm_data().GetPoses();
    }

    if (!InitTracksAndLandmarks())
      return false;

    if (!sfm_data_.GetPoses().empty())
    {
      const bool bTriangulation = Triangulation();
      Save(sfm_data_, stlplus::create_filespec(sOut_directory_, "Initialization", ".ply"), ESfM_Data(ALL));
      RemoveOutliers_AngleError(sfm_data_, Square(2.0));
      RemoveOutliers_PixelResidualError(sfm_data_, Square(4.0));

      //-- Display some statistics
      std::cout << "\n\n-------------------------------" << "\n"
        << "-- Starting Structure from Motion (statistics) with:\n"
        << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
        << " from " << sfm_data_.GetViews().size() << " input images.\n"
        << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
        << "-------------------------------" << "\n";

      if (!bTriangulation)
      {
        // We are not able to triangulate enough 3D points for the given poses.
        return false;
      }
    }
  }

  if (sfm_data_.GetPoses().empty())
  {
    return false;
  }

  //--
  //- 2. While we can localize some cameras in the reconstruction
  //     a. Triangulate the landmarks
  //     b. Perform Bundle Adjustment and cleaning
  //--
  IndexT resection_round = 0;

  // Incrementally estimate the pose of the cameras based on a confidence score.
  // The confidence score is based on the track_inlier_ratio.
  // First the camera with the most of 2D-3D overlap are added the we added
  // the one with lower confidence.
  const std::array<float, 2> track_inlier_ratios = {0.2, 0.0};
  for (const float track_inlier_ratio : track_inlier_ratios)
  {
    while (AddingMissingView(track_inlier_ratio))
    {
      // Create new 3D points
      Triangulation();
      // Adjust the scene
      BundleAdjustment();
      // Remove unstable triangulations and camera poses
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      eraseUnstablePosesAndObservations(sfm_data_);
      ++resection_round;
    }
  }

  //--
  //- 3. Final bundle Adjustment
  //--
  BundleAdjustment();

  //-- Reconstruction done.
  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
    << " from " << sfm_data_.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
    << "-- #Poses loop used: " << resection_round << "\n"
    << "-------------------------------" << "\n";

  return true;
}

bool SequentialSfMReconstructionEngine2::InitTracksAndLandmarks()
{
  // Compute tracks from matches
  tracks::TracksBuilder tracksBuilder;
  {
    tracksBuilder.Build(matches_provider_->pairWise_matches_);
    tracksBuilder.Filter();
    tracksBuilder.ExportToSTL(map_tracks_);

    std::cout << "\n" << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats :
      //    - number of images
      //    - number of tracks
      std::set<uint32_t> set_imagesId;
      tracks::TracksUtilsMap::ImageIdInTracks(map_tracks_, set_imagesId);
      osTrack
        << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.cbegin(),
        set_imagesId.cend(),
        std::ostream_iterator<uint32_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<uint32_t, uint32_t> map_Occurence_TrackLength;
      tracks::TracksUtilsMap::TracksLength(map_tracks_, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (const auto & it : map_Occurence_TrackLength)  {
        osTrack << "\t" << it.first << "\t" << it.second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Init the putative landmarks
  {
    // For every track add the obervations:
    // - views and feature positions that see this landmark
    for ( const auto & iterT : map_tracks_ )
    {
      Observations obs;
      const size_t track_Length = iterT.second.size();
      for (const auto & track_ids : iterT.second) // {ViewId, FeatureId}
      {
        const auto & view_id = track_ids.first;
        const auto & feat_id = track_ids.second;
        const Vec2 x = features_provider_->feats_per_view[view_id][feat_id].coords().cast<double>();
        obs.insert({view_id, Observation(x, feat_id)});
      }
      landmarks_[iterT.first].obs = std::move(obs);
    }
  }

  // Initialize the shared track visibility helper
  shared_track_visibility_helper_.reset(new openMVG::tracks::SharedTrackVisibilityHelper(map_tracks_));
  return map_tracks_.size() > 0;
}

bool SequentialSfMReconstructionEngine2::Triangulation()
{
  sfm_data_.structure = landmarks_;

  //--
  // Triangulation
  //--
  // Clean the structure:
  //  - keep observations that are linked to valid pose and intrinsic data.

  const double max_reprojection_error = 4.0;
  const IndexT min_required_inliers = 2;
  const IndexT min_sample_index = 2;
  eraseObservationsWithMissingPoses(sfm_data_, min_sample_index);
  SfM_Data_Structure_Computation_Robust triangulation_engine(
      max_reprojection_error,
      min_required_inliers,
      min_sample_index);

  triangulation_engine.triangulate(sfm_data_);

  return !sfm_data_.structure.empty();
}

bool SequentialSfMReconstructionEngine2::AddingMissingView
(
  const float & track_inlier_ratio
)
{
  if (sfm_data_.GetLandmarks().empty())
    return false;

  // Collect the views that does not have any 3D pose
  const std::set<IndexT> view_with_no_pose = [&]
  {
    std::set<IndexT> idx;
    for (const auto & view_it : sfm_data_.GetViews())
    {
      const View * v = view_it.second.get();
      const IndexT id_pose = v->id_pose;
      if (sfm_data_.GetPoses().count(id_pose) == 0)
        idx.insert(view_it.first);
    }
    return idx;
  }();

  const IndexT pose_before = sfm_data_.GetPoses().size();

  // Get the track ids of the reconstructed landmarks
  const std::vector<IndexT> reconstructed_trackId = [&]
  {
    std::vector<IndexT> tracks_ids;
    tracks_ids.reserve(sfm_data_.GetLandmarks().size());
    std::transform(sfm_data_.GetLandmarks().cbegin(), sfm_data_.GetLandmarks().cend(),
      std::back_inserter(tracks_ids),
      stl::RetrieveKey());
    std::sort(tracks_ids.begin(), tracks_ids.end());
    return tracks_ids;
  }();

  // List the view that have a sufficient 2D-3D coverage for robust pose estimation
#pragma omp parallel
  for (const auto & view_id : view_with_no_pose)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      // List the track related to the current view_id
      openMVG::tracks::STLMAPTracks view_tracks;
      shared_track_visibility_helper_->GetTracksInImages({view_id}, view_tracks);
      std::set<IndexT> view_tracks_ids;
      tracks::TracksUtilsMap::GetTracksIdVector(view_tracks, &view_tracks_ids);

      // Get the ids of the already reconstructed tracks
      const std::set<IndexT> track_id_for_resection = [&]
      {
        std::set<IndexT> track_id;
        std::set_intersection(view_tracks_ids.cbegin(), view_tracks_ids.cend(),
          reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
          std::inserter(track_id, track_id.begin()));
        return track_id;
      }();

      const double track_ratio = track_id_for_resection.size() / static_cast<float>(view_tracks_ids.size() + 1);
      std::cout
        << "ViewId: " << view_id
        << "; #number of 2D-3D matches: " << track_id_for_resection.size()
        << "; " << track_ratio * 100 << " % of the view track coverage."
        << std::endl;

      if (!track_id_for_resection.empty() && track_ratio > track_inlier_ratio)
      {
        // Get feat_id for the 2D/3D associations
        std::vector<IndexT> feature_id_for_resection;
        tracks::TracksUtilsMap::GetFeatIndexPerViewAndTrackId(
          view_tracks,
          track_id_for_resection,
          view_id,
          &feature_id_for_resection);

        // Localize the image inside the SfM reconstruction
        Image_Localizer_Match_Data resection_data;
        resection_data.pt2D.resize(2, track_id_for_resection.size());
        resection_data.pt3D.resize(3, track_id_for_resection.size());

        // Look if the intrinsic data is known or not
        const View * view = sfm_data_.GetViews().at(view_id).get();
        std::shared_ptr<cameras::IntrinsicBase> intrinsic;
        if (sfm_data_.GetIntrinsics().count(view->id_intrinsic))
        {
          intrinsic = sfm_data_.GetIntrinsics().at(view->id_intrinsic);
        }

        // Collect the feature observation
        Mat2X pt2D_original(2, track_id_for_resection.size());
        auto track_it = track_id_for_resection.cbegin();
        auto feat_it = feature_id_for_resection.cbegin();
        for (size_t cpt = 0; cpt < track_id_for_resection.size(); ++cpt, ++track_it, ++feat_it)
        {
          resection_data.pt3D.col(cpt) = sfm_data_.GetLandmarks().at(*track_it).X;
          resection_data.pt2D.col(cpt) = pt2D_original.col(cpt) =
            features_provider_->feats_per_view.at(view_id)[*feat_it].coords().cast<double>();
          // Handle image distortion if intrinsic is known (to ease the resection)
          if (intrinsic && intrinsic->have_disto())
          {
            resection_data.pt2D.col(cpt) = intrinsic->get_ud_pixel(resection_data.pt2D.col(cpt));
          }
        }

        geometry::Pose3 pose;
        const bool bResection = sfm::SfM_Localizer::Localize
        (
          intrinsic ? resection::SolverType::P3P_KE_CVPR17 : resection::SolverType::DLT_6POINTS,
          {view->ui_width, view->ui_height},
          intrinsic ? intrinsic.get() : nullptr,
          resection_data,
          pose
        );
        resection_data.pt2D = std::move(pt2D_original); // restore original image domain points

        const float inlier_ratio = resection_data.vec_inliers.size()/static_cast<float>(feature_id_for_resection.size());
        std::cout
          << std::endl
          << "-------------------------------" << "\n"
          << "-- Robust Resection of camera index: <" << view_id << "> image: "
          <<  view->s_Img_path <<"\n"
          << "-- Threshold: " << resection_data.error_max << "\n"
          << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "\n"
          << "-- Nb points used for Resection: " << feature_id_for_resection.size() << "\n"
          << "-- Nb points validated by robust estimation: " << resection_data.vec_inliers.size() << "\n"
          << "-- % points validated: "
          << inlier_ratio * 100 << "\n"
          << "-------------------------------" << std::endl;

        // Refine the pose of the found camera pose by using a BA and fix 3D points.
        if (bResection && inlier_ratio > 0.5)
        {
          // A valid pose has been found (try to refine it):
          // If no valid intrinsic as input:
          //  init a new one from the projection matrix decomposition
          // Else use the existing one and consider it as constant.
          if (!intrinsic)
          {
            // setup a default camera model from the found projection matrix
            Mat3 K, R;
            Vec3 t;
            KRt_From_P(resection_data.projection_matrix, &K, &R, &t);

            const double focal = (K(0,0) + K(1,1))/2.0;
            const Vec2 principal_point(K(0,2), K(1,2));

            // Create the new camera intrinsic
            switch (cam_type_)
            {
              case PINHOLE_CAMERA:
                intrinsic = std::make_shared<Pinhole_Intrinsic>
                  (view->ui_width, view->ui_height, focal, principal_point(0), principal_point(1));
              break;
              case PINHOLE_CAMERA_RADIAL1:
                intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
                  (view->ui_width, view->ui_height, focal, principal_point(0), principal_point(1));
              break;
              case PINHOLE_CAMERA_RADIAL3:
                intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
                  (view->ui_width, view->ui_height, focal, principal_point(0), principal_point(1));
              break;
              case PINHOLE_CAMERA_BROWN:
                intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
                  (view->ui_width, view->ui_height, focal, principal_point(0), principal_point(1));
              break;
              case PINHOLE_CAMERA_FISHEYE:
                intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
                  (view->ui_width, view->ui_height, focal, principal_point(0), principal_point(1));
              break;
              default:
                std::cerr << "Try to create an unknown camera type." << std::endl;
            }
          }
          const bool b_refine_pose = true;
          const bool b_refine_intrinsics = false;
          if (intrinsic && sfm::SfM_Localizer::RefinePose(
              intrinsic.get(), pose,
              resection_data, b_refine_pose, b_refine_intrinsics))
          {
            // - intrinsic parameters (if the view has no intrinsic group add a new one)
            if (sfm_data_.intrinsics.count(sfm_data_.views.at(view_id)->id_intrinsic) == 0)
            {
              // Since the view have not yet an intrinsic group before, create a new one
              IndexT new_intrinsic_id = 0;
              if (!sfm_data_.GetIntrinsics().empty())
              {
                // Since some intrinsic Id already exists,
                //  we have to create a new unique identifier following the existing one
                std::set<IndexT> existing_intrinsic_id;
                  std::transform(sfm_data_.GetIntrinsics().cbegin(), sfm_data_.GetIntrinsics().cend(),
                  std::inserter(existing_intrinsic_id, existing_intrinsic_id.begin()),
                  stl::RetrieveKey());
                new_intrinsic_id = (*existing_intrinsic_id.rbegin()) + 1;
              }
              #pragma omp critical
              {
                sfm_data_.views.at(view_id)->id_intrinsic = new_intrinsic_id;
                sfm_data_.intrinsics[new_intrinsic_id] = intrinsic;
              }
            }

            // Update the found camera pose
            #pragma omp critical
            sfm_data_.poses[view->id_pose] = pose;
          }
        }
      }
    }
  }

  const IndexT pose_after = sfm_data_.GetPoses().size();
  return (pose_after != pose_before);
}

bool SequentialSfMReconstructionEngine2::BundleAdjustment()
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
