// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"

#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/stl/stlMap.hpp"

#include "third_party/histogram/histogram.hpp"

namespace openMVG{
namespace sfm{

using namespace openMVG::rotation_averaging;

Pair_Set GlobalSfM_Rotation_AveragingSolver::GetUsedPairs() const
{
  return used_pairs;
}

bool GlobalSfM_Rotation_AveragingSolver::Run(
  ERotationAveragingMethod eRotationAveragingMethod,
  ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod,
  const RelativeRotations & relativeRot_In,
  Hash_Map<IndexT, Mat3> & map_globalR
) const
{
  RelativeRotations relativeRotations = relativeRot_In;
  // We work on a copy, since inference can remove some relative motions

  switch (eRelativeRotationInferenceMethod)
  {
    case TRIPLET_ROTATION_INFERENCE_NONE:
    break;
    case TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR:
    {
      //-------------------
      // Triplet inference (test over the composition error)
      //-------------------
      Pair_Set pairs = getPairs(relativeRotations);
      std::vector< graph::Triplet > vec_triplets = graph::tripletListing(pairs);
      //-- Rejection triplet that are 'not' identity rotation (error to identity > 5°)
      TripletRotationRejection(5.0f, vec_triplets, relativeRotations);

      pairs = getPairs(relativeRotations);
      const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
      if (set_remainingIds.empty())
        return false;
      KeepOnlyReferencedElement(set_remainingIds, relativeRotations);
    }
  break;
    default:
    std::cerr
      << "Unknown relative rotation inference method: "
      << (int) eRelativeRotationInferenceMethod << std::endl;
  }

  // Compute contiguous index (mapping between sparse index and contiguous index)
  //  from ranging in [min(Id), max(Id)] to  [0, nbCam]

  const Pair_Set pairs = getPairs(relativeRotations);
  Hash_Map<IndexT, IndexT> reindexForward, reindexBackward;
  reindex(pairs, reindexForward, reindexBackward);

  for (RelativeRotations::iterator iter = relativeRotations.begin();  iter != relativeRotations.end(); ++iter)
  {
    RelativeRotation & rel = *iter;
    rel.i = reindexForward[rel.i];
    rel.j = reindexForward[rel.j];
  }

  //- B. solve global rotation computation
  bool bSuccess = false;
  std::vector<Mat3> vec_globalR(reindexForward.size());
  switch (eRotationAveragingMethod)
  {
    case ROTATION_AVERAGING_L2:
    {
      //- Solve the global rotation estimation problem:
      bSuccess = rotation_averaging::l2::L2RotationAveraging(
        reindexForward.size(),
        relativeRotations,
        vec_globalR);
      //- Non linear refinement of the global rotations
      if (bSuccess)
        bSuccess = rotation_averaging::l2::L2RotationAveraging_Refine(
          relativeRotations,
          vec_globalR);

      // save kept pairs (restore original pose indices using the backward reindexing)
      for (RelativeRotations::iterator iter = relativeRotations.begin();  iter != relativeRotations.end(); ++iter)
      {
        RelativeRotation & rel = *iter;
        rel.i = reindexBackward[rel.i];
        rel.j = reindexBackward[rel.j];
      }
      used_pairs = getPairs(relativeRotations);
    }
    break;
    case ROTATION_AVERAGING_L1:
    {
      using namespace openMVG::rotation_averaging::l1;

      //- Solve the global rotation estimation problem:
      const size_t nMainViewID = 0; //arbitrary choice
      std::vector<bool> vec_inliers;
      bSuccess = rotation_averaging::l1::GlobalRotationsRobust(
        relativeRotations, vec_globalR, nMainViewID, 0.0f, &vec_inliers);

      std::cout << "\ninliers: " << std::endl;
      std::copy(vec_inliers.begin(), vec_inliers.end(), std::ostream_iterator<bool>(std::cout, " "));
      std::cout << std::endl;

      // save kept pairs (restore original pose indices using the backward reindexing)
      for (size_t i = 0; i < vec_inliers.size(); ++i)
      {
        if (vec_inliers[i])
        {
          used_pairs.insert(
            Pair(reindexBackward[relativeRotations[i].i],
                 reindexBackward[relativeRotations[i].j]));
        }
      }
    }
    break;
    default:
    std::cerr
      << "Unknown rotation averaging method: "
      << (int) eRotationAveragingMethod << std::endl;
  }

  if (bSuccess)
  {
    //-- Setup the averaged rotations
    for (size_t i = 0; i < vec_globalR.size(); ++i)  {
      map_globalR[reindexBackward[i]] = vec_globalR[i];
    }
  }
  else {
    std::cerr << "Global rotation solving failed." << std::endl;
  }

  return bSuccess;
}

/// Reject edges of the view graph that do not produce triplets with tiny
///  angular error once rotation composition have been computed.
void GlobalSfM_Rotation_AveragingSolver::TripletRotationRejection(
  const double max_angular_error,
  std::vector< graph::Triplet > & vec_triplets,
  RelativeRotations & relativeRotations) const
{
  const size_t edges_start_count = relativeRotations.size();

  RelativeRotations_map map_relatives = getMap(relativeRotations);
  RelativeRotations_map map_relatives_validated;

  //--
  // ROTATION OUTLIERS DETECTION
  //--

  std::vector< graph::Triplet > vec_triplets_validated;
  vec_triplets_validated.reserve(vec_triplets.size());

  std::vector<float> vec_errToIdentityPerTriplet;
  vec_errToIdentityPerTriplet.reserve(vec_triplets.size());
  // Compute the composition error for each length 3 cycles
  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graph::Triplet & triplet = vec_triplets[i];
    const IndexT I = triplet.i, J = triplet.j , K = triplet.k;

    //-- Find the three relative rotations
    const Pair ij(I,J), ji(J,I);
    const Mat3 RIJ = (map_relatives.count(ij)) ?
      map_relatives.at(ij).Rij : Mat3(map_relatives.at(ji).Rij.transpose());

    const Pair jk(J,K), kj(K,J);
    const Mat3 RJK = (map_relatives.count(jk)) ?
      map_relatives.at(jk).Rij : Mat3(map_relatives.at(kj).Rij.transpose());

    const Pair ki(K,I), ik(I,K);
    const Mat3 RKI = (map_relatives.count(ki)) ?
      map_relatives.at(ki).Rij : Mat3(map_relatives.at(ik).Rij.transpose());

    const Mat3 Rot_To_Identity = RIJ * RJK * RKI; // motion composition
    const float angularErrorDegree = static_cast<float>(R2D(getRotationMagnitude(Rot_To_Identity)));
    vec_errToIdentityPerTriplet.push_back(angularErrorDegree);

    if (angularErrorDegree < max_angular_error)
    {
      vec_triplets_validated.push_back(triplet);

      if (map_relatives.count(ij))
        map_relatives_validated[ij] = map_relatives.at(ij);
      else
        map_relatives_validated[ji] = map_relatives.at(ji);

      if (map_relatives.count(jk))
        map_relatives_validated[jk] = map_relatives.at(jk);
      else
        map_relatives_validated[kj] = map_relatives.at(kj);

      if (map_relatives.count(ki))
        map_relatives_validated[ki] = map_relatives.at(ki);
      else
        map_relatives_validated[ik] = map_relatives.at(ik);
    }
  }
  map_relatives = std::move(map_relatives_validated);

  // update to keep only useful triplets
  relativeRotations.clear();
  relativeRotations.reserve(map_relatives.size());
  std::transform(map_relatives.begin(), map_relatives.end(), std::back_inserter(relativeRotations), stl::RetrieveValue());
  std::transform(map_relatives.begin(), map_relatives.end(), std::inserter(used_pairs, used_pairs.begin()), stl::RetrieveKey());

  // Display statistics about rotation triplets error:
  std::cout << "\nStatistics about rotation triplets:" << std::endl;
  minMaxMeanMedian<float>(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  std::sort(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  if (!vec_errToIdentityPerTriplet.empty())
  {
    Histogram<float> histo(0.0f, *max_element(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end()), 20);
    histo.Add(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());
    std::cout << histo.ToString() << std::endl;
  }

  {
    std::cout << "\nTriplets filtering based on composition error on unit cycles\n";
    std::cout << "#Triplets before: " << vec_triplets.size() << "\n"
    << "#Triplets after: " << vec_triplets_validated.size() << std::endl;
  }

  vec_triplets = std::move(vec_triplets_validated);

  const size_t edges_end_count = relativeRotations.size();
  std::cout << "\n #Edges removed by triplet inference: " << edges_start_count - edges_end_count << std::endl;
}

} // namespace sfm
} // namespace openMVG
