
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP

namespace openMVG{
namespace globalSfM{

enum ETranslationAveragingMethod
{
  TRANSLATION_AVERAGING_L1 = 1,
  TRANSLATION_AVERAGING_L2 = 2
};

}
}

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/graph/graph.hpp"

namespace openMVG{
namespace globalSfM{

class GlobalSfM_Translation_AveragingSolver
{
  RelativeInfo_Vec vec_initialRijTijEstimates;

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
    const Features_Provider * normalized_features_provider,
    const matching::PairWiseMatches & map_Matches_E,
    const Hash_Map<IndexT, Mat3> & map_globalR);

  void Compute_translations(
    const SfM_Data & sfm_data,
    const Features_Provider * normalized_features_provider,
    const matching::PairWiseMatches & map_Matches_E,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    matching::PairWiseMatches &tripletWise_matches);

  //-- Perform a trifocal estimation of the graph contain in vec_triplets with an
  // edge coverage algorithm. It's complexity is sub-linear in term of edges count.
  void ComputePutativeTranslation_EdgesCoverage(
    const SfM_Data & sfm_data,
    const Hash_Map<IndexT, Mat3> & map_globalR,
    const Features_Provider * normalized_features_provider,
    const matching::PairWiseMatches & map_Matches_E,
    const std::vector< graphUtils::Triplet > & vec_triplets,
    RelativeInfo_Vec & vec_initialEstimates,
    matching::PairWiseMatches & newpairMatches);

  // Robust estimation and refinement of a translation and 3D points of an image triplets.
  bool Estimate_T_triplet(
    const openMVG::tracks::STLMAPTracks & map_tracksCommon,
    const Features_Provider * features_provider,
    const std::vector<Mat3> & vec_global_R_Triplet,
    std::vector<Vec3> & vec_tis,
    double & dPrecision, // UpperBound of the precision found by the AContrario estimator
    std::vector<size_t> & vec_inliers,
    const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
    const size_t nI,
    const size_t nJ,
    const size_t nK,
    const std::string & sOutDirectory) const;
};

} // namespace globalSfM
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_PIPELINES_GLOBAL_TRANSLATION_AVERAGING_HPP
