// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_PIPELINES_SFKP_STRUCTURE_ESTIMATOR_HPP
#define OPENMVG_SFM_PIPELINES_SFKP_STRUCTURE_ESTIMATOR_HPP

#include <memory>

#include "openMVG/matching/indMatch.hpp"

namespace openMVG { namespace sfm { struct Regions_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

class SfM_Data_Structure_Estimation_From_Known_Poses
{
public:

  explicit SfM_Data_Structure_Estimation_From_Known_Poses
  (
    double max_reprojection_error // pixels
  );


  /// Use geometry of the views to compute a putative structure from features and descriptors.
  void run(
    SfM_Data & sfm_data,
    const Pair_Set & pairs,
    const std::shared_ptr<Regions_Provider> & regions_provider);

private:

  /// Use guided matching to find corresponding 2-view correspondences
  void match(
    const SfM_Data & sfm_data,
    const Pair_Set & pairs,
    const std::shared_ptr<Regions_Provider> & regions_provider);

  /// Filter inconsistent correspondences by using 3-view correspondences on view triplets
  void filter(
    const SfM_Data & sfm_data,
    const Pair_Set & pairs,
    const std::shared_ptr<Regions_Provider> & regions_provider);

  /// Init & triangulate landmark observations from validated 3-view correspondences
  void triangulate(
    SfM_Data & sfm_data,
    const std::shared_ptr<Regions_Provider> & regions_provider);

private:
  //--
  // DATA (temporary)
  //--
  matching::PairWiseMatches putatives_matches;
  matching::PairWiseMatches triplets_matches;
  double max_reprojection_error_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_PIPELINES_SFKP_STRUCTURE_ESTIMATOR_HPP
