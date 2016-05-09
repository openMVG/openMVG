
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP

namespace openMVG{
namespace sfm{

enum ETranslationAveragingMethod
{
  TRANSLATION_AVERAGING_L1 = 1,
  TRANSLATION_AVERAGING_L2_DISTANCE_CHORDAL = 2,
  TRANSLATION_AVERAGING_SOFTL1 = 3
};

} // namespace sfm
} // namespace openMVG

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/graph/graph.hpp"

namespace openMVG{
namespace sfm{

class GlobalSfM_Translation_AveragingSolver
{
  std::vector<RelativeInfo_Vec> vec_relative_motion_;

public:

  /// Use features in normalized camera frames
  bool Run(
    ETranslationAveragingMethod eTranslationAveragingMethod,
    SfM_Data & sfm_data,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches & tripletWise_matches
  );

private:
  bool Translation_averaging(
    ETranslationAveragingMethod eTranslationAveragingMethod,
    SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR);

  void Compute_translations(
    const SfM_Data & sfm_data,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches &tripletWise_matches);

  //-- Compute the relative translations on the rotations graph.
  // Compute relative translations by using triplets of poses.
  // Use an edge coverage algorithm to reduce the graph covering complexity
  // Complexity: sub-linear in term of edges count.
  void ComputePutativeTranslation_EdgesCoverage(
    const SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    std::vector<RelativeInfo_Vec> & vec_triplet_relative_motion,
    matching::PairWiseMatches & newpairMatches);

  // Robust estimation and refinement of triplet of translations
  bool Estimate_T_triplet(
    const SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const Features_Provider * normalized_features_provider,
    const Matches_Provider * matches_provider,
    const graph::Triplet & poses_id,
    std::vector<Vec3> & vec_tis,
    double & dPrecision, // UpperBound of the precision found by the AContrario estimator
    std::vector<size_t> & vec_inliers,
    openMVG::tracks::STLMAPTracks & rig_tracks,
    const std::string & sOutDirectory) const;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
