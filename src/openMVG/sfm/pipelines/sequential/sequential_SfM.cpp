// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#include <ceres/types.h>
#include <functional>
#include <iostream>
#include <utility>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;

SequentialSfMReconstructionEngine::SequentialSfMReconstructionEngine(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory),
    sLogging_file_(sloggingFile),
    initial_pair_(Pair(0,0)),
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
  }
  // Init remaining image list
  for (Views::const_iterator itV = sfm_data.GetViews().begin();
    itV != sfm_data.GetViews().end(); ++itV)
  {
    set_remaining_view_id_.insert(itV->second->id_view);
  }
}

SequentialSfMReconstructionEngine::~SequentialSfMReconstructionEngine()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void SequentialSfMReconstructionEngine::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void SequentialSfMReconstructionEngine::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

bool SequentialSfMReconstructionEngine::Process() {

  //-------------------
  //-- Incremental reconstruction
  //-------------------

  if (!InitLandmarkTracks())
    return false;

  // Initial pair choice
  if (initial_pair_ == Pair(0,0))
  {
    if (!AutomaticInitialPairChoice(initial_pair_))
    {
      // Cannot find a valid initial pair, try to set it by hand?
      if (!ChooseInitialPair(initial_pair_))
      {
        return false;
      }
    }
  }
  // Else a starting pair was already initialized before

  // Initial pair Essential Matrix and [R|t] estimation.
  if (!MakeInitialPair3D(initial_pair_))
    return false;

  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  size_t resectionGroupIndex = 0;
  std::vector<uint32_t> vec_possible_resection_indexes;
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (const auto & iter : vec_possible_resection_indexes)
    {
      bImageAdded |= Resection(iter);
      set_remaining_view_id_.erase(iter);
    }

    if (bImageAdded)
    {
      // Scene logging as ply for visual debug
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
      Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

      // Perform BA until all point are under the given precision
      do
      {
        BundleAdjustment();
      }
      while (badTrackRejector(4.0, 50));
      eraseUnstablePosesAndObservations(sfm_data_);
    }
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  if (badTrackRejector(4.0, 0))
  {
    eraseUnstablePosesAndObservations(sfm_data_);
  }

  //-- Reconstruction done.
  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
    << " from " << sfm_data_.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  std::cout << "\nHistogram of residuals:" << h.ToString() << std::endl;

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
    std::pair< std::pair<double,double>, std::pair<double,double> > range =
      autoJSXGraphViewport<double>(xBin, h.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("3DtoImageResiduals",600,300);
    jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());
  }
  return true;
}

/// Select a candidate initial pair
bool SequentialSfMReconstructionEngine::ChooseInitialPair(Pair & initialPairIndex) const
{
  if (initial_pair_ != Pair(0,0))
  {
    // Internal initial pair is already initialized (so return it)
    initialPairIndex = initial_pair_;
  }
  else
  {
    // List Views that supports valid intrinsic
    std::set<IndexT> valid_views;
    for (Views::const_iterator it = sfm_data_.GetViews().begin();
      it != sfm_data_.GetViews().end(); ++it)
    {
      const View * v = it->second.get();
      if (sfm_data_.GetIntrinsics().find(v->id_intrinsic) != sfm_data_.GetIntrinsics().end())
        valid_views.insert(v->id_view);
    }

    if (sfm_data_.GetIntrinsics().empty() || valid_views.empty())
    {
      std::cerr
        << "There is no defined intrinsic data in order to compute an essential matrix for the initial pair."
        << std::endl;
      return false;
    }

    std::cout << std::endl
      << "----------------------------------------------------\n"
      << "SequentialSfMReconstructionEngine::ChooseInitialPair\n"
      << "----------------------------------------------------\n"
      << " Pairs that have valid intrinsic and high support of points are displayed:\n"
      << " Choose one pair manually by typing the two integer indexes\n"
      << "----------------------------------------------------\n"
      << std::endl;

    // Try to list the 10 top pairs that have:
    //  - valid intrinsics,
    //  - valid estimated Fundamental matrix.
    std::vector< uint32_t > vec_NbMatchesPerPair;
    std::vector<openMVG::matching::PairWiseMatches::const_iterator> vec_MatchesIterator;
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
    for (openMVG::matching::PairWiseMatches::const_iterator
      iter = map_Matches.begin();
      iter != map_Matches.end(); ++iter)
    {
      const Pair current_pair = iter->first;
      if (valid_views.count(current_pair.first) &&
        valid_views.count(current_pair.second) )
      {
        vec_NbMatchesPerPair.push_back(iter->second.size());
        vec_MatchesIterator.push_back(iter);
      }
    }
    // sort the Pairs in descending order according their correspondences count
    using namespace stl::indexed_sort;
    std::vector< sort_index_packet_descend< uint32_t, uint32_t> > packet_vec(vec_NbMatchesPerPair.size());
    sort_index_helper(packet_vec, &vec_NbMatchesPerPair[0], std::min((size_t)10, vec_NbMatchesPerPair.size()));

    for (size_t i = 0; i < std::min((size_t)10, vec_NbMatchesPerPair.size()); ++i) {
      const uint32_t index = packet_vec[i].index;
      openMVG::matching::PairWiseMatches::const_iterator iter = vec_MatchesIterator[index];
      std::cout << "(" << iter->first.first << "," << iter->first.second <<")\t\t"
        << iter->second.size() << " matches" << std::endl;
    }

    // Ask the user to choose an initial pair (by set some view ids)
    std::cout << std::endl << " type INITIAL pair ids: X enter Y enter\n";
    int val, val2;
    if ( std::cin >> val && std::cin >> val2) {
      initialPairIndex.first = val;
      initialPairIndex.second = val2;
    }
  }

  std::cout << "\nPutative starting pair is: (" << initialPairIndex.first
      << "," << initialPairIndex.second << ")" << std::endl;

  // Check validity of the initial pair indices:
  if (features_provider_->feats_per_view.find(initialPairIndex.first) == features_provider_->feats_per_view.end() ||
      features_provider_->feats_per_view.find(initialPairIndex.second) == features_provider_->feats_per_view.end())
  {
    std::cerr << "At least one of the initial pair indices is invalid."
      << std::endl;
    return false;
  }
  return true;
}

bool SequentialSfMReconstructionEngine::InitLandmarkTracks()
{
  // Compute tracks from matches
  tracks::TracksBuilder tracksBuilder;

  {
    // List of features matches for each couple of images
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
    std::cout << "\n" << "Track building" << std::endl;

    tracksBuilder.Build(map_Matches);
    std::cout << "\n" << "Track filtering" << std::endl;
    tracksBuilder.Filter();
    std::cout << "\n" << "Track export to internal struct" << std::endl;
    //-- Build tracks with STL compliant type :
    tracksBuilder.ExportToSTL(map_tracks_);

    std::cout << "\n" << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats :
      //    - number of images
      //    - number of tracks
      std::set<uint32_t> set_imagesId;
      tracks::TracksUtilsMap::ImageIdInTracks(map_tracks_, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
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
  // Initialize the shared track visibility helper
  shared_track_visibility_helper_.reset(new openMVG::tracks::SharedTrackVisibilityHelper(map_tracks_));
  return map_tracks_.size() > 0;
}

bool SequentialSfMReconstructionEngine::AutomaticInitialPairChoice(Pair & initial_pair) const
{
  // select a pair that have the largest baseline (mean angle between it's bearing vectors).

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

  std::vector<std::pair<double, Pair> > scoring_per_pair;

  // Compute the relative pose & the 'baseline score'
  C_Progress_display my_progress_bar( matches_provider_->pairWise_matches_.size(),
    std::cout,
    "Automatic selection of an initial pair:\n" );
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (const std::pair< Pair, IndMatches > & match_pair : matches_provider_->pairWise_matches_)
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
        const View * view_I = sfm_data_.GetViews().at(I).get();
        const Intrinsics::const_iterator iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic);
        const View * view_J = sfm_data_.GetViews().at(J).get();
        const Intrinsics::const_iterator iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

        const Pinhole_Intrinsic * cam_I = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_I->second.get());
        const Pinhole_Intrinsic * cam_J = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_J->second.get());
        if (cam_I != nullptr && cam_J != nullptr)
        {
          openMVG::tracks::STLMAPTracks map_tracksCommon;
          shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

          // Copy points correspondences to arrays for relative pose estimation
          const size_t n = map_tracksCommon.size();
          Mat xI(2,n), xJ(2,n);
          size_t cptIndex = 0;
          for (openMVG::tracks::STLMAPTracks::const_iterator
            iterT = map_tracksCommon.begin(); iterT != map_tracksCommon.end();
            ++iterT, ++cptIndex)
          {
            tracks::submapTrack::const_iterator iter = iterT->second.begin();
            const uint32_t i = iter->second;
            const uint32_t j = (++iter)->second;

            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
            xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
          }

          // Robust estimation of the relative pose
          RelativePose_Info relativePose_info;
          relativePose_info.initial_residual_tolerance = Square(4.0);

          if (robustRelativePose(
            cam_I->K(), cam_J->K(),
            xI, xJ, relativePose_info,
            {cam_I->w(), cam_I->h()}, {cam_J->w(), cam_J->h()},
            256) && relativePose_info.vec_inliers.size() > iMin_inliers_count)
          {
            // Triangulate inliers & compute angle between bearing vectors
            std::vector<float> vec_angles;
            vec_angles.reserve(relativePose_info.vec_inliers.size());
            const Pose3 pose_I = Pose3(Mat3::Identity(), Vec3::Zero());
            const Pose3 pose_J = relativePose_info.relativePose;
            const Mat34 PI = cam_I->get_projective_equivalent(pose_I);
            const Mat34 PJ = cam_J->get_projective_equivalent(pose_J);
            for (const uint32_t & inlier_idx : relativePose_info.vec_inliers)
            {
              Vec3 X;
              TriangulateDLT(
                PI, xI.col(inlier_idx).homogeneous(),
                PJ, xJ.col(inlier_idx).homogeneous(), &X);

              openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
              std::advance(iterT, inlier_idx);
              tracks::submapTrack::const_iterator iter = iterT->second.begin();
              const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
              const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
              vec_angles.push_back(AngleBetweenRay(pose_I, cam_I, pose_J, cam_J, featI, featJ));
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

/// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
bool SequentialSfMReconstructionEngine::MakeInitialPair3D(const Pair & current_pair)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  const uint32_t I = std::min(current_pair.first, current_pair.second);
  const uint32_t J = std::max(current_pair.first, current_pair.second);

  // a. Assert we have valid pinhole cameras
  const View * view_I = sfm_data_.GetViews().at(I).get();
  const Intrinsics::const_iterator iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic);
  const View * view_J = sfm_data_.GetViews().at(J).get();
  const Intrinsics::const_iterator iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

  if (iterIntrinsic_I == sfm_data_.GetIntrinsics().end() ||
      iterIntrinsic_J == sfm_data_.GetIntrinsics().end() )
  {
    return false;
  }

  const Pinhole_Intrinsic * cam_I = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_I->second.get());
  const Pinhole_Intrinsic * cam_J = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_J->second.get());
  if (cam_I == nullptr || cam_J == nullptr)
  {
    return false;
  }

  // b. Get common features between the two view
  // use the track to have a more dense match correspondence set
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  shared_track_visibility_helper_->GetTracksInImages({I, J}, map_tracksCommon);

  //-- Copy point to arrays
  const size_t n = map_tracksCommon.size();
  Mat xI(2,n), xJ(2,n);
  uint32_t cptIndex = 0;
  for (openMVG::tracks::STLMAPTracks::const_iterator
    iterT = map_tracksCommon.begin(); iterT != map_tracksCommon.end();
    ++iterT, ++cptIndex)
  {
    tracks::submapTrack::const_iterator iter = iterT->second.begin();
    const uint32_t i = iter->second;
    const uint32_t j = (++iter)->second;

    Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
    xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
    feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
    xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
  }

  // c. Robust estimation of the relative pose
  RelativePose_Info relativePose_info;

  const std::pair<size_t, size_t>
    imageSize_I(cam_I->w(), cam_I->h()),
    imageSize_J(cam_J->w(), cam_J->h());

  if (!robustRelativePose(
    cam_I->K(), cam_J->K(), xI, xJ, relativePose_info, imageSize_I, imageSize_J, 4096))
  {
    std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
      << std::endl;
    return false;
  }
  std::cout << "A-Contrario initial pair residual: "
    << relativePose_info.found_residual_precision << std::endl;
  // Bound min precision at 1 pix.
  relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

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

    for (openMVG::tracks::STLMAPTracks::const_iterator
      iterT = map_tracksCommon.begin();
      iterT != map_tracksCommon.end();
      ++iterT)
    {
      // Get corresponding points
      tracks::submapTrack::const_iterator iter = iterT->second.begin();
      const uint32_t i = iter->second;
      const uint32_t j = (++iter)->second;

      const Vec2 x1_ = features_provider_->feats_per_view[I][i].coords().cast<double>();
      const Vec2 x2_ = features_provider_->feats_per_view[J][j].coords().cast<double>();

      Vec3 X;
      TriangulateDLT(P1, x1_.homogeneous(), P2, x2_.homogeneous(), &X);
      Observations obs;
      obs[view_I->id_view] = Observation(x1_, i);
      obs[view_J->id_view] = Observation(x2_, j);
      landmarks[iterT->first].obs = std::move(obs);
      landmarks[iterT->first].X = X;
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
    for (Landmarks::const_iterator iter = tiny_scene.GetLandmarks().begin();
      iter != tiny_scene.GetLandmarks().end(); ++iter)
    {
      const IndexT trackId = iter->first;
      const Landmark & landmark = iter->second;
      const Observations & obs = landmark.obs;
      Observations::const_iterator iterObs_xI = obs.find(view_I->id_view);
      Observations::const_iterator iterObs_xJ = obs.find(view_J->id_view);

      const Observation & ob_xI = iterObs_xI->second;
      const Observation & ob_xJ = iterObs_xJ->second;

      const double angle = AngleBetweenRay(
        pose_I, cam_I, pose_J, cam_J, ob_xI.x, ob_xJ.x);
      const Vec2 residual_I = cam_I->residual(pose_I, landmark.X, ob_xI.x);
      const Vec2 residual_J = cam_J->residual(pose_J, landmark.X, ob_xJ.x);
      if ( angle > 2.0 &&
           pose_I.depth(landmark.X) > 0 &&
           pose_J.depth(landmark.X) > 0 &&
           residual_I.norm() < relativePose_info.found_residual_precision &&
           residual_J.norm() < relativePose_info.found_residual_precision)
      {
        sfm_data_.structure[trackId] = landmarks[trackId];
      }
    }
    // Save outlier residual information
    Histogram<double> histoResiduals;
    std::cout << std::endl
      << "=========================\n"
      << " MSE Residual InitialPair Inlier: " << ComputeResidualsHistogram(&histoResiduals) << "\n"
      << "=========================" << std::endl;

    if (!sLogging_file_.empty())
    {
      using namespace htmlDocument;
      html_doc_stream_->pushInfo(htmlMarkup("h1","Essential Matrix."));
      std::ostringstream os;
      os << std::endl
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

      std::vector<double> xBin = histoResiduals.GetXbinsValue();
      std::pair< std::pair<double,double>, std::pair<double,double> > range =
        autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

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
        "Reconstruction_Report.html").c_str());
      htmlFileStream << html_doc_stream_->getDoc();
    }
  }
  return !sfm_data_.structure.empty();
}

double SequentialSfMReconstructionEngine::ComputeResidualsHistogram(Histogram<double> * histo)
{
  // Collect residuals for each observation
  std::vector<float> vec_residuals;
  vec_residuals.reserve(sfm_data_.structure.size());
  for (Landmarks::const_iterator iterTracks = sfm_data_.GetLandmarks().begin();
      iterTracks != sfm_data_.GetLandmarks().end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data_.GetViews().find(itObs->first)->second.get();
      const Pose3 pose = sfm_data_.GetPoseOrDie(view);
      const std::shared_ptr<IntrinsicBase> intrinsic = sfm_data_.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      vec_residuals.push_back( std::abs(residual(0)) );
      vec_residuals.push_back( std::abs(residual(1)) );
    }
  }
  // Display statistics
  if (vec_residuals.size() > 1)
  {
    float dMin, dMax, dMean, dMedian;
    minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                            dMin, dMax, dMean, dMedian);
    if (histo)  {
      *histo = Histogram<double>(dMin, dMax, 10);
      histo->Add(vec_residuals.begin(), vec_residuals.end());
    }

    std::cout << std::endl << std::endl;
    std::cout << std::endl
      << "SequentialSfMReconstructionEngine::ComputeResidualsMSE." << "\n"
      << "\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size() << std::endl
      << "\t-- Residual min:\t" << dMin << std::endl
      << "\t-- Residual median:\t" << dMedian << std::endl
      << "\t-- Residual max:\t "  << dMax << std::endl
      << "\t-- Residual mean:\t " << dMean << std::endl;

    return dMean;
  }
  return -1.0;
}

/// Functor to sort a vector of pair given the pair's second value
template<class T1, class T2, class Pred = std::less<T2> >
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left,
                    const std::pair<T1,T2>&right)
  {
    Pred p;
    return p(left.second, right.second);
  }
};

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
bool SequentialSfMReconstructionEngine::FindImagesWithPossibleResection(
  std::vector<uint32_t> & vec_possible_indexes)
{
  // Threshold used to select the best images
  static const float dThresholdGroup = 0.75f;

  vec_possible_indexes.clear();

  if (set_remaining_view_id_.empty() || sfm_data_.GetLandmarks().empty())
    return false;

  // Collect tracksIds
  std::set<uint32_t> reconstructed_trackId;
  std::transform(sfm_data_.GetLandmarks().begin(), sfm_data_.GetLandmarks().end(),
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
  std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<uint32_t, uint32_t, std::greater<uint32_t> >());

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

/**
 * @brief Add one image to the 3D reconstruction. To the resectioning of
 * the camera and triangulate all the new possible tracks.
 * @param[in] viewIndex: image index to add to the reconstruction.
 *
 * A. Compute 2D/3D matches
 * B. Look if intrinsic data is known or not
 * C. Do the resectioning: compute the camera pose.
 * D. Refine the pose of the found camera
 * E. Update the global scene with the new camera
 * F. Update the observations into the global scene structure
 * G. Triangulate new possible 2D tracks
 */
bool SequentialSfMReconstructionEngine::Resection(const uint32_t viewIndex)
{
  using namespace tracks;

  // A. Compute 2D/3D matches
  // A1. list tracks ids used by the view
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  shared_track_visibility_helper_->GetTracksInImages({viewIndex}, map_tracksCommon);
  std::set<uint32_t> set_tracksIds;
  TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

  // A2. intersects the track list with the reconstructed
  std::set<uint32_t> reconstructed_trackId;
  std::transform(sfm_data_.GetLandmarks().begin(), sfm_data_.GetLandmarks().end(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  // Get the ids of the already reconstructed tracks
  std::set<uint32_t> set_trackIdForResection;
  std::set_intersection(set_tracksIds.cbegin(), set_tracksIds.cend(),
    reconstructed_trackId.cbegin(), reconstructed_trackId.cend(),
    std::inserter(set_trackIdForResection, set_trackIdForResection.begin()));

  if (set_trackIdForResection.empty())
  {
    // No match. The image has no connection with already reconstructed points.
    std::cout << std::endl
      << "-------------------------------" << "\n"
      << "-- Resection of camera index: " << viewIndex << "\n"
      << "-- Resection status: " << "FAILED" << "\n"
      << "-------------------------------" << std::endl;
    return false;
  }

  // Get back featId associated to a tracksID already reconstructed.
  // These 2D/3D associations will be used for the resection.
  std::vector<uint32_t> vec_featIdForResection;
  TracksUtilsMap::GetFeatIndexPerViewAndTrackId(map_tracksCommon,
    set_trackIdForResection,
    viewIndex,
    &vec_featIdForResection);

  // Localize the image inside the SfM reconstruction
  Image_Localizer_Match_Data resection_data;
  resection_data.pt2D.resize(2, set_trackIdForResection.size());
  resection_data.pt3D.resize(3, set_trackIdForResection.size());

  // B. Look if intrinsic data is known or not
  const View * view_I = sfm_data_.GetViews().at(viewIndex).get();
  std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
  if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic))
  {
    optional_intrinsic = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic);
  }

  // Setup the track 2d observation for this new view
  Mat2X pt2D_original(2, set_trackIdForResection.size());
  std::set<uint32_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
  std::vector<uint32_t>::const_iterator iterfeatId = vec_featIdForResection.begin();
  for (size_t cpt = 0; cpt < vec_featIdForResection.size(); ++cpt, ++iterTrackId, ++iterfeatId)
  {
    resection_data.pt3D.col(cpt) = sfm_data_.GetLandmarks().at(*iterTrackId).X;
    resection_data.pt2D.col(cpt) = pt2D_original.col(cpt) =
      features_provider_->feats_per_view.at(viewIndex)[*iterfeatId].coords().cast<double>();
    // Handle image distortion if intrinsic is known (to ease the resection)
    if (optional_intrinsic && optional_intrinsic->have_disto())
    {
      resection_data.pt2D.col(cpt) = optional_intrinsic->get_ud_pixel(resection_data.pt2D.col(cpt));
    }
  }

  // C. Do the resectioning: compute the camera pose
  std::cout << std::endl
    << "-------------------------------" << std::endl
    << "-- Robust Resection of view: " << viewIndex << std::endl;

  geometry::Pose3 pose;
  const bool bResection = sfm::SfM_Localizer::Localize
  (
    Pair(view_I->ui_width, view_I->ui_height),
    optional_intrinsic.get(),
    resection_data,
    pose
  );
  resection_data.pt2D = std::move(pt2D_original); // restore original image domain points

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Resection of Image index: <" << viewIndex << "> image: "
      << view_I->s_Img_path <<"<br> \n";
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << std::endl
      << "-------------------------------" << "<br>"
      << "-- Robust Resection of camera index: <" << viewIndex << "> image: "
      <<  view_I->s_Img_path <<"<br>"
      << "-- Threshold: " << resection_data.error_max << "<br>"
      << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "<br>"
      << "-- Nb points used for Resection: " << vec_featIdForResection.size() << "<br>"
      << "-- Nb points validated by robust estimation: " << resection_data.vec_inliers.size() << "<br>"
      << "-- % points validated: "
      << resection_data.vec_inliers.size()/static_cast<float>(vec_featIdForResection.size()) << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  if (!bResection)
    return false;

  // D. Refine the pose of the found camera.
  // We use a local scene with only the 3D points and the new camera.
  {
    const bool b_new_intrinsic = (optional_intrinsic == nullptr);
    // A valid pose has been found (try to refine it):
    // If no valid intrinsic as input:
    //  init a new one from the projection matrix decomposition
    // Else use the existing one and consider it as constant.
    if (b_new_intrinsic)
    {
      // setup a default camera model from the found projection matrix
      Mat3 K, R;
      Vec3 t;
      KRt_From_P(resection_data.projection_matrix, &K, &R, &t);

      const double focal = (K(0,0) + K(1,1))/2.0;
      const Vec2 principal_point(K(0,2), K(1,2));

      // Create the new camera intrinsic group
      switch (cam_type_)
      {
        case PINHOLE_CAMERA:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_RADIAL1:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_RADIAL3:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_BROWN:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_FISHEYE:
            optional_intrinsic =
                std::make_shared<Pinhole_Intrinsic_Fisheye>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        default:
          std::cerr << "Try to create an unknown camera type." << std::endl;
          return false;
      }
    }
    const bool b_refine_pose = true;
    const bool b_refine_intrinsics = false;
    if (!sfm::SfM_Localizer::RefinePose(
        optional_intrinsic.get(), pose,
        resection_data, b_refine_pose, b_refine_intrinsics))
    {
      return false;
    }

    // E. Update the global scene with:
    // - the new found camera pose
    sfm_data_.poses[view_I->id_pose] = pose;
    // - track the view's AContrario robust estimation found threshold
    map_ACThreshold_.insert({viewIndex, resection_data.error_max});
    // - intrinsic parameters (if the view has no intrinsic group add a new one)
    if (b_new_intrinsic)
    {
      // Since the view have not yet an intrinsic group before, create a new one
      IndexT new_intrinsic_id = 0;
      if (!sfm_data_.GetIntrinsics().empty())
      {
        // Since some intrinsic Id already exists,
        //  we have to create a new unique identifier following the existing one
        std::set<IndexT> existing_intrinsicId;
          std::transform(sfm_data_.GetIntrinsics().begin(), sfm_data_.GetIntrinsics().end(),
          std::inserter(existing_intrinsicId, existing_intrinsicId.begin()),
          stl::RetrieveKey());
        new_intrinsic_id = (*existing_intrinsicId.rbegin())+1;
      }
      sfm_data_.views.at(viewIndex)->id_intrinsic = new_intrinsic_id;
      sfm_data_.intrinsics[new_intrinsic_id] = optional_intrinsic;
    }
  }

  // F. List tracks that share content with this view and add observations and new 3D track if required.
  //    - If the track already exists (look if the new view tracks observation are valid)
  //    - If the track does not exists, try robust triangulation & add the new valid view track observation
  {
    // Get information of new view
    const IndexT I = viewIndex;
    const View * view_I = sfm_data_.GetViews().at(I).get();
    const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
    const Pose3 pose_I = sfm_data_.GetPoseOrDie(view_I);

    // Vector of all already reconstructed views
    const std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);

    // Go through each track and look if we must add new view observations or new 3D points
    for (const std::pair< uint32_t, tracks::submapTrack >& trackIt : map_tracksCommon)
    {
      const uint32_t trackId = trackIt.first;
      const tracks::submapTrack & track = trackIt.second;

      // List the potential view observations of the track
      const tracks::submapTrack & allViews_of_track = map_tracks_[trackId];

      // List to save the new view observations that must be added to the track
      std::set<IndexT> new_track_observations_valid_views;

      // If the track was already reconstructed
      if (sfm_data_.structure.count(trackId) != 0)
      {
        // Since the 3D point was triangulated before we add the new the Inth view observation
        new_track_observations_valid_views.insert(I);
      }
      else
      {
        // Go through the views that observe this track & look if a successful triangulation can be done
        for (const std::pair< IndexT, IndexT >& trackViewIt : allViews_of_track)
        {
          const IndexT & J = trackViewIt.first;
          // If view is valid try triangulation
          if (J!=I && valid_views.count(J) != 0 )
          {
            // If successfuly triangulated add the observation from J view
            if (sfm_data_.structure.count(trackId) != 0)
            {
              new_track_observations_valid_views.insert(J);
            }
            else
            {
              const View * view_J = sfm_data_.GetViews().at(J).get();
              const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
              const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
              const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();

              // Position of the point in view I
              const Vec2 xI = features_provider_->feats_per_view.at(I)[track.at(I)].coords().cast<double>();

              // Try to triangulate a 3D point from J view
              // A new 3D point must be added
              // Triangulate it
              const Vec2 xI_ud = cam_I->get_ud_pixel(xI);
              const Vec2 xJ_ud = cam_J->get_ud_pixel(xJ);
              const Mat34 P_I = cam_I->get_projective_equivalent(pose_I);
              const Mat34 P_J = cam_J->get_projective_equivalent(pose_J);
              Vec3 X = Vec3::Zero();
              TriangulateDLT(P_I, xI_ud.homogeneous(), P_J, xJ_ud.homogeneous(), &X);
              // Check triangulation result
              const double angle = AngleBetweenRay(pose_I, cam_I, pose_J, cam_J, xI, xJ);
              const Vec2 residual_I = cam_I->residual(pose_I, X, xI);
              const Vec2 residual_J = cam_J->residual(pose_J, X, xJ);
              if (
                  //  - Check angle (small angle leads to imprecise triangulation)
                  angle > 2.0 &&
                  //  - Check positive depth
                  pose_I.depth(X) > 0 &&
                  pose_J.depth(X) > 0 &&
                  //  - Check residual values (must be inferior to the found view's AContrario threshold)
                  residual_I.norm() < std::max(4.0, map_ACThreshold_.at(I)) &&
                  residual_J.norm() < std::max(4.0, map_ACThreshold_.at(J))
                 )
              {
                // Add a new track
                Landmark & landmark = sfm_data_.structure[trackId];
                landmark.X = X;
                new_track_observations_valid_views.insert(I);
                new_track_observations_valid_views.insert(J);
              } // 3D point is valid
              else
              {
                // We mark the view to add the observations once the point is triangulated
                new_track_observations_valid_views.insert(J);
              } // 3D point is invalid
            }
          }
        }// Go through all the views
      }// If new point

      // If successfuly triangulated, add the valid view observations
      if (sfm_data_.structure.count(trackId) != 0 &&
          !new_track_observations_valid_views.empty()
         )
      {
        Landmark & landmark = sfm_data_.structure[trackId];
        // Check if view feature point observations of the track are valid (residual, depth) or not
        for (const IndexT &J: new_track_observations_valid_views)
        {
          const View * view_J = sfm_data_.GetViews().at(J).get();
          const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
          const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
          const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();

          const Vec2 residual = cam_J->residual(pose_J, landmark.X, xJ);
          if (pose_J.depth(landmark.X) > 0 &&
              residual.norm() < std::max(4.0, map_ACThreshold_.at(J))
             )
          {
            landmark.obs[J] = Observation(xJ, allViews_of_track.at(J));
          }
        }
      }
    }// All the tracks in the view
  }
  return true;
}

/// Bundle adjustment to refine Structure; Motion and Intrinsics
bool SequentialSfMReconstructionEngine::BundleAdjustment()
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

/**
 * @brief Discard tracks with too large residual error
 *
 * Remove observation/tracks that have:
 *  - too large residual error
 *  - too small angular value
 *
 * @return True if more than 'count' outliers have been removed.
 */
bool SequentialSfMReconstructionEngine::badTrackRejector(double dPrecision, size_t count)
{
  const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(sfm_data_, dPrecision, 2);
  const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(sfm_data_, 2.0);

  return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
}

} // namespace sfm
} // namespace openMVG
