// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/structure_from_known_poses/structure_estimator.hpp"

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/robust_estimation/guided_matching.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/tracks/tracks.hpp"

#include "third_party/progress/progress_display.hpp"

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::features;
using namespace openMVG::geometry;
using namespace openMVG::matching;

/// Camera pair epipole (Projection of camera center 2 in the image plane 1)
inline Vec3 epipole_from_P(const Mat34& P1, const Pose3& P2)
{
  return P1 * P2.center().homogeneous();
}

/// Export point feature based vector to a matrix [(x,y)'T, (x,y)'T]
/// Use the camera intrinsics in order to get undistorted pixel coordinates
template<typename MatT >
void PointsToMat(
  const IntrinsicBase * cam,
  const PointFeatures & vec_feats,
  MatT & m)
{
  m.resize(2, vec_feats.size());

  Mat::Index i = 0;
  for (PointFeatures::const_iterator iter = vec_feats.begin();
    iter != vec_feats.end(); ++iter, ++i)
  {
    if (cam)
      m.col(i) = cam->get_ud_pixel({iter->x(), iter->y()});
    else
      m.col(i) << iter->x(), iter->y();
  }
}

SfM_Data_Structure_Estimation_From_Known_Poses::SfM_Data_Structure_Estimation_From_Known_Poses
(
  double max_reprojection_error // pixels
):
  max_reprojection_error_(max_reprojection_error)
{
}

/// Use geometry of the views to compute a putative structure from features and descriptors.
void SfM_Data_Structure_Estimation_From_Known_Poses::run(
  SfM_Data & sfm_data,
  const Pair_Set & pairs,
  const std::shared_ptr<Regions_Provider> & regions_provider)
{
  sfm_data.structure.clear();

  match(sfm_data, pairs, regions_provider);
  filter(sfm_data, pairs, regions_provider);
  triangulate(sfm_data, regions_provider);
}

/// Use guided matching to find corresponding 2-view correspondences
void SfM_Data_Structure_Estimation_From_Known_Poses::match(
  const SfM_Data & sfm_data,
  const Pair_Set & pairs,
  const std::shared_ptr<Regions_Provider> & regions_provider)
{
  C_Progress_display my_progress_bar( pairs.size(), std::cout,
    "Compute pairwise fundamental guided matching:\n" );
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif // OPENMVG_USE_OPENMP
  for (const Pair & pair : pairs)
  {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif // OPENMVG_USE_OPENMP
    {
    // --
    // Perform GUIDED MATCHING
    // --
    // Use the computed model to check valid correspondences
    // - by considering geometric error and descriptor distance ratio.
    std::vector<IndMatch> vec_corresponding_indexes;

    const View
      * viewL = sfm_data.GetViews().at(pair.first).get(),
      * viewR = sfm_data.GetViews().at(pair.second).get();

    const Pose3
      poseL = sfm_data.GetPoseOrDie(viewL),
      poseR = sfm_data.GetPoseOrDie(viewR);

    if (sfm_data.GetIntrinsics().count(viewL->id_intrinsic) != 0 ||
        sfm_data.GetIntrinsics().count(viewR->id_intrinsic) != 0)
    {
      const Intrinsics::const_iterator
        iterIntrinsicL = sfm_data.GetIntrinsics().find(viewL->id_intrinsic),
        iterIntrinsicR = sfm_data.GetIntrinsics().find(viewR->id_intrinsic);

      const Mat34
        P_L = iterIntrinsicL->second->get_projective_equivalent(poseL),
        P_R = iterIntrinsicR->second->get_projective_equivalent(poseR);

      const Mat3 F_lr = F_from_P(P_L, P_R);
      const double thresholdF = max_reprojection_error_;

      const std::shared_ptr<features::Regions>
        regionsL = regions_provider->get(pair.first),
        regionsR = regions_provider->get(pair.second);

    #if defined(EXHAUSTIVE_MATCHING)
      geometry_aware::GuidedMatching
        <Mat3, openMVG::fundamental::kernel::EpipolarDistanceError>
        (
          F_lr,
          iterIntrinsicL->second.get(),
          *regionsL.get(),
          iterIntrinsicR->second.get(),
          *regionsR.get(),
          Square(thresholdF), Square(0.8),
          vec_corresponding_indexes
        );
    #else
      const Vec3 epipole2  = epipole_from_P(P_R, poseL);

      geometry_aware::GuidedMatching_Fundamental_Fast
        <openMVG::fundamental::kernel::EpipolarDistanceError>
        (
          F_lr,
          epipole2,
          iterIntrinsicL->second.get(),
          *regionsL.get(),
          iterIntrinsicR->second.get(),
          *regionsR.get(),
          iterIntrinsicR->second->w(), iterIntrinsicR->second->h(),
          Square(thresholdF), Square(0.8),
          vec_corresponding_indexes
        );
  #endif

  #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
  #endif // OPENMVG_USE_OPENMP
        {
          putatives_matches[pair].insert(putatives_matches[pair].end(),
            vec_corresponding_indexes.begin(), vec_corresponding_indexes.end());
        }
        ++my_progress_bar;
      }
    }
  }
}

/// Filter inconsistent correspondences by using 3-view correspondences on view triplets
void SfM_Data_Structure_Estimation_From_Known_Poses::filter(
  const SfM_Data & sfm_data,
  const Pair_Set & pairs,
  const std::shared_ptr<Regions_Provider> & regions_provider)
{
  // Compute triplets
  // Triangulate triplet tracks
  //  - keep valid one

  using Triplets = std::vector<graph::Triplet>;
  const Triplets triplets = graph::TripletListing(pairs);

  C_Progress_display my_progress_bar( triplets.size(), std::cout,
    "Per triplet tracks validation (discard spurious correspondences):\n" );
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif // OPENMVG_USE_OPENMP
  for (Triplets::const_iterator it = triplets.begin(); it != triplets.end(); ++it)
  {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif // OPENMVG_USE_OPENMP
    {
      ++my_progress_bar;

      const graph::Triplet & triplet = *it;
      const IndexT I = triplet.i, J = triplet.j , K = triplet.k;

      openMVG::tracks::STLMAPTracks map_tracksCommon;
      openMVG::tracks::TracksBuilder tracksBuilder;
      {
        PairWiseMatches map_matchesIJK;
        if (putatives_matches.count({I,J}))
          map_matchesIJK.insert(*putatives_matches.find({I,J}));

        if (putatives_matches.count({I,K}))
          map_matchesIJK.insert(*putatives_matches.find({I,K}));

        if (putatives_matches.count({J,K}))
          map_matchesIJK.insert(*putatives_matches.find({J,K}));

        if (map_matchesIJK.size() >= 2) {
          tracksBuilder.Build(map_matchesIJK);
          tracksBuilder.Filter(3);
          tracksBuilder.ExportToSTL(map_tracksCommon);
        }

        const std::map<IndexT, std::shared_ptr<openMVG::features::Regions>> regions =
        {{I, regions_provider->get(I)},
         {J, regions_provider->get(J)},
         {K, regions_provider->get(K)},
        };

        // Triangulate the tracks
        for (const auto & track_it : map_tracksCommon)
        {
          const tracks::submapTrack & subTrack = track_it.second;
          std::vector<Vec3> bearing;
          std::vector<Mat34> poses;
          bearing.reserve(subTrack.size());
          poses.reserve(subTrack.size());
          for (const auto & observation_it : subTrack) {
            const size_t imaIndex = observation_it.first;
            const size_t featIndex = observation_it.second;
            const View * view = sfm_data.GetViews().at(imaIndex).get();
            const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
            const Pose3 pose = sfm_data.GetPoseOrDie(view);
            const Vec2 pt = regions.at(imaIndex)->GetRegionPosition(featIndex);
            bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
            poses.emplace_back(pose.asMatrix());
          }
          const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
          Vec4 Xhomogeneous;
          TriangulateNViewAlgebraic(bearing_matrix, poses, &Xhomogeneous);
          const Vec3 X = Xhomogeneous.hnormalized();

          // Test validity of the hypothesis:
          // - residual error
          // - chierality
          bool bChierality = true;
          bool bReprojection_error = true;
          int i(0);
          for (tracks::submapTrack::const_iterator obs_it = subTrack.begin();
            obs_it != subTrack.end() && bChierality && bReprojection_error; ++obs_it, ++i)
          {
            const View * view = sfm_data.views.at(obs_it->first).get();

            const Pose3 pose = sfm_data.GetPoseOrDie(view);
            bChierality &= CheiralityTest(bearing[i], pose, X);

            const size_t imaIndex = obs_it->first;
            const size_t featIndex = obs_it->second;
            const Vec2 pt = regions.at(imaIndex)->GetRegionPosition(featIndex);
            const IntrinsicBase * cam = sfm_data.intrinsics.at(view->id_intrinsic).get();
            const Vec2 residual = cam->residual(pose(X), pt);
            bReprojection_error &= residual.squaredNorm() < max_reprojection_error_;
          }
          if (bChierality && bReprojection_error)
          // TODO: Add an angular check ?
          {
            #ifdef OPENMVG_USE_OPENMP
            #pragma omp critical
            #endif // OPENMVG_USE_OPENMP
            {
              openMVG::tracks::submapTrack::const_iterator iterI, iterJ, iterK;
              iterI = iterJ = iterK = subTrack.begin();
              std::advance(iterJ,1);
              std::advance(iterK,2);

              triplets_matches[{I,J}].emplace_back(iterI->second, iterJ->second);
              triplets_matches[{J,K}].emplace_back(iterJ->second, iterK->second);
            }
          }
        }
      }
    }
  }
  // Clear putatives matches since they are no longer required
  matching::PairWiseMatches().swap(putatives_matches);
}

/// Init & triangulate landmark observations from validated 3-view correspondences
void SfM_Data_Structure_Estimation_From_Known_Poses::triangulate(
  SfM_Data & sfm_data,
  const std::shared_ptr<Regions_Provider> & regions_provider)
{
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  openMVG::tracks::TracksBuilder tracksBuilder;
  tracksBuilder.Build(triplets_matches);
  tracksBuilder.Filter(3);
  tracksBuilder.ExportToSTL(map_tracksCommon);
  matching::PairWiseMatches().swap(triplets_matches);

  // Generate new Structure tracks
  sfm_data.structure.clear();

  SfM_Data_Structure_Computation_Robust structure_estimator(max_reprojection_error_);
  C_Progress_display my_progress_bar( map_tracksCommon.size(), std::cout,
    "Tracks to structure conversion:\n" );
  // Fill sfm_data with the computed tracks (no 3D yet)
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif // OPENMVG_USE_OPENMP
  for (int i = 0; i < map_tracksCommon.size(); ++i)
  {
    ++my_progress_bar;

    tracks::STLMAPTracks::const_iterator itTracks = map_tracksCommon.begin();
    std::advance(itTracks, i);
    {
      const tracks::submapTrack & track = itTracks->second;

      Observations obs;
      for (const auto & track_obs : track)
      {
        const IndexT imaIndex = track_obs.first;
        const IndexT featIndex = track_obs.second;
        const std::shared_ptr<features::Regions> regions = regions_provider->get(imaIndex);
        const Vec2 pt = regions->GetRegionPosition(featIndex);
        obs[imaIndex] = Observation(pt, featIndex);
      }
      Landmark landmark;
      if (structure_estimator.robust_triangulation(sfm_data, obs, landmark))
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
      #endif // OPENMVG_USE_OPENMP
      {
        sfm_data.structure[itTracks->first] = landmark;
      }
    }
  }
}

} // namespace sfm
} // namespace openMVG
