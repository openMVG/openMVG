// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP

#include <string>
#include <vector>

#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/tracks/tracks.hpp"

namespace openMVG { namespace graph { struct Triplet; } }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

enum ETranslationAveragingMethod
{
  TRANSLATION_AVERAGING_L1 = 1,
  TRANSLATION_AVERAGING_L2_DISTANCE_CHORDAL = 2,
  TRANSLATION_AVERAGING_SOFTL1 = 3
};

struct SfM_Data;
struct Matches_Provider;
struct Features_Provider;

class GlobalSfM_Translation_AveragingSolver
{
  std::vector<RelativeInfo_Vec> vec_relative_motion_;

public:

  bool Run(
    ETranslationAveragingMethod eTranslationAveragingMethod,
    openMVG::sfm::SfM_Data & sfm_data,
    const openMVG::sfm::Features_Provider * features_provider,
    const openMVG::sfm::Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches & tripletWise_matches
  );

private:
  bool Translation_averaging(
    ETranslationAveragingMethod eTranslationAveragingMethod,
    sfm::SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR);

  void Compute_translations(
    const sfm::SfM_Data & sfm_data,
    const sfm::Features_Provider * features_provider,
    const sfm::Matches_Provider * matches_provider,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches &tripletWise_matches);

  //-- Compute the relative translations on the rotations graph.
  // Compute relative translations by using triplets of poses.
  // Use an edge coverage algorithm to reduce the graph covering complexity
  // Complexity: sub-linear in term of edges count.
  void ComputePutativeTranslation_EdgesCoverage(
    const sfm::SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const sfm::Features_Provider * features_provider,
    const sfm::Matches_Provider * matches_provider,
    std::vector<RelativeInfo_Vec> & vec_triplet_relative_motion,
    matching::PairWiseMatches & newpairMatches);

  // Robust estimation and refinement of triplet of translations
  bool Estimate_T_triplet(
    const sfm::SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const sfm::Features_Provider * features_provider,
    const sfm::Matches_Provider * matches_provider,
    const graph::Triplet & poses_id,
    std::vector<Vec3> & vec_tis,
    double & dPrecision, // UpperBound of the precision found by the AContrario estimator
    std::vector<uint32_t> & vec_inliers,
    openMVG::tracks::STLMAPTracks & rig_tracks,
    const std::string & sOutDirectory) const;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
