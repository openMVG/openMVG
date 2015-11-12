// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/structure_from_known_poses/structure_estimator.hpp"

#include "openMVG/matching/metric.hpp"
#include "openMVG/robust_estimation/guided_matching.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::features;
using namespace openMVG::geometry;


/// Camera pair epipole (Projection of camera center 2 in the image plane 1)
static Vec3 epipole_from_P(const Mat34& P1, const Pose3& P2)
{
  const Vec3 c = P2.center();
  Vec4 center;
  center << c(0), c(1), c(2), 1.0;
  return P1*center;
}

/// Export point feature based vector to a matrix [(x,y)'T, (x,y)'T]
/// Use the camera intrinsics in order to get undistorted pixel coordinates
template<typename MatT >
static void PointsToMat(
  const IntrinsicBase * cam,
  const PointFeatures & vec_feats,
  MatT & m)
{
  m.resize(2, vec_feats.size());
  typedef typename MatT::Scalar Scalar; // Output matrix type

  size_t i = 0;
  for( PointFeatures::const_iterator iter = vec_feats.begin();
    iter != vec_feats.end(); ++iter, ++i)
  {
    if (cam)
      m.col(i) = cam->get_ud_pixel(Vec2(iter->x(), iter->y()));
    else
      m.col(i) << iter->x(), iter->y();
  }
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
  for (Pair_Set::const_iterator it = pairs.begin(); it != pairs.end(); ++it)
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

    const View * viewL = sfm_data.GetViews().at(it->first).get();
    const Pose3 poseL = sfm_data.GetPoseOrDie(viewL);
    const Intrinsics::const_iterator iterIntrinsicL = sfm_data.GetIntrinsics().find(viewL->id_intrinsic);
    const View * viewR = sfm_data.GetViews().at(it->second).get();
    const Pose3 poseR = sfm_data.GetPoseOrDie(viewR);
    const Intrinsics::const_iterator iterIntrinsicR = sfm_data.GetIntrinsics().find(viewR->id_intrinsic);

    if (sfm_data.GetIntrinsics().count(viewL->id_intrinsic) != 0 ||
        sfm_data.GetIntrinsics().count(viewR->id_intrinsic) != 0)
    {
      const Mat34 P_L = iterIntrinsicL->second.get()->get_projective_equivalent(poseL);
      const Mat34 P_R = iterIntrinsicR->second.get()->get_projective_equivalent(poseR);

      const Mat3 F_lr = F_from_P(P_L, P_R);
      const double thresholdF = 4.0;

    #if defined(EXHAUSTIVE_MATCHING)
      geometry_aware::GuidedMatching
        <Mat3, openMVG::fundamental::kernel::EpipolarDistanceError>
        (
          F_lr,
          iterIntrinsicL->second.get(),
          *regions_provider->regions_per_view.at(it->first),
          iterIntrinsicR->second.get(),
          *regions_provider->regions_per_view.at(it->second),
          Square(thresholdF), Square(0.8),
          vec_corresponding_indexes
        );
    #else
      const Vec3 epipole2  = epipole_from_P(P_R, poseL);

      const features::Regions * regions = regions_provider->regions_per_view.at(it->first).get();
      geometry_aware::GuidedMatching_Fundamental_Fast
        <openMVG::fundamental::kernel::EpipolarDistanceError>
        (
          F_lr,
          epipole2,
          iterIntrinsicL->second.get(),
          *regions_provider->regions_per_view.at(it->first),
          iterIntrinsicR->second.get(),
          *regions_provider->regions_per_view.at(it->second),
          iterIntrinsicR->second->w(), iterIntrinsicR->second->h(),
          Square(thresholdF), Square(0.8),
          vec_corresponding_indexes
        );
  #endif

  #ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
  #endif // OPENMVG_USE_OPENMP
        {
          ++my_progress_bar;
          for (size_t i = 0; i < vec_corresponding_indexes.size(); ++i)
            putatives_matches[*it].push_back(vec_corresponding_indexes[i]);
        }
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

  typedef std::vector< graph::Triplet > Triplets;
  const Triplets triplets = graph::tripletListing(pairs);

  C_Progress_display my_progress_bar( triplets.size(), std::cout,
    "Per triplet tracks validation (discard spurious correspondences):\n" );
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel
#endif // OPENMVG_USE_OPENMP
  for( Triplets::const_iterator it = triplets.begin(); it != triplets.end(); ++it)
  {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif // OPENMVG_USE_OPENMP
    {
      #ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
      #endif // OPENMVG_USE_OPENMP
      {++my_progress_bar;}

      const graph::Triplet & triplet = *it;
      const IndexT I = triplet.i, J = triplet.j , K = triplet.k;

      openMVG::tracks::STLMAPTracks map_tracksCommon;
      openMVG::tracks::TracksBuilder tracksBuilder;
      {
        PairWiseMatches map_matchesIJK;
        if(putatives_matches.find(std::make_pair(I,J)) != putatives_matches.end())
          map_matchesIJK.insert(*putatives_matches.find(std::make_pair(I,J)));

        if(putatives_matches.find(std::make_pair(I,K)) != putatives_matches.end())
          map_matchesIJK.insert(*putatives_matches.find(std::make_pair(I,K)));

        if(putatives_matches.find(std::make_pair(J,K)) != putatives_matches.end())
          map_matchesIJK.insert(*putatives_matches.find(std::make_pair(J,K)));

        if (map_matchesIJK.size() >= 2) {
          tracksBuilder.Build(map_matchesIJK);
          tracksBuilder.Filter(3, false);
          tracksBuilder.ExportToSTL(map_tracksCommon);
        }

        // Triangulate the tracks
        for (tracks::STLMAPTracks::const_iterator iterTracks = map_tracksCommon.begin();
          iterTracks != map_tracksCommon.end(); ++iterTracks) {
          {
            const tracks::submapTrack & subTrack = iterTracks->second;
            Triangulation trianObj;
            for (tracks::submapTrack::const_iterator iter = subTrack.begin(); iter != subTrack.end(); ++iter) {
              const size_t imaIndex = iter->first;
              const size_t featIndex = iter->second;
              const View * view = sfm_data.GetViews().at(imaIndex).get();
              const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
              const Pose3 pose = sfm_data.GetPoseOrDie(view);
              const Vec2 pt = regions_provider->regions_per_view.at(imaIndex)->GetRegionPosition(featIndex);
              trianObj.add(cam->get_projective_equivalent(pose), cam->get_ud_pixel(pt));
            }
            const Vec3 Xs = trianObj.compute();
            if (trianObj.minDepth() > 0 && trianObj.error()/(double)trianObj.size() < 4.0)
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

                triplets_matches[std::make_pair(I,J)].push_back(IndMatch(iterI->second, iterJ->second));
                triplets_matches[std::make_pair(J,K)].push_back(IndMatch(iterJ->second, iterK->second));
                triplets_matches[std::make_pair(I,K)].push_back(IndMatch(iterI->second, iterK->second));
              }
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

  // Fill sfm_data with the computed tracks (no 3D yet)
  Landmarks & structure = sfm_data.structure;
  IndexT idx(0);
  for (tracks::STLMAPTracks::const_iterator itTracks = map_tracksCommon.begin();
    itTracks != map_tracksCommon.end();
    ++itTracks, ++idx)
  {
    const tracks::submapTrack & track = itTracks->second;
    structure[idx] = Landmark();
    Observations & obs = structure.at(idx).obs;
    for (tracks::submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
    {
      const size_t imaIndex = it->first;
      const size_t featIndex = it->second;
      const Vec2 pt = regions_provider->regions_per_view.at(imaIndex)->GetRegionPosition(featIndex);
      obs[imaIndex] = Observation(pt, featIndex);
    }
  }

  // Triangulate them using a robust triangulation scheme
  SfM_Data_Structure_Computation_Robust structure_estimator(true);
  structure_estimator.triangulate(sfm_data);
}

} // namespace sfm
} // namespace openMVG

