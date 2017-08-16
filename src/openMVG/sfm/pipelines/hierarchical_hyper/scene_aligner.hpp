// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

namespace ceres {
  class Problem;
}

namespace openMVG {
namespace sfm {

struct SfM_Data;

class SceneAligner
{
public:
  bool computeTransformAndDestinationSeparators
  (SfM_Data & destination_sfm_data,
    const SfM_Data & sfm_data_first, // first submap scene
    const SfM_Data & sfm_data_second, // second submap scene
    std::vector<double> & second_base_node_pose, // second base node pose with respect to the first base node (what we're looking for)
    double & scaling_factor, // scaling factor to complete the transformation
    const std::vector<size_t> & common_track_ids // which tracks id are commonly reconstructed in both submaps
  );

  SceneAligner(Bundle_Adjustment_Ceres::BA_Ceres_options options = Bundle_Adjustment_Ceres::BA_Ceres_options());

protected:
  Bundle_Adjustment_Ceres::BA_Ceres_options ceres_options_;

  virtual void configureProblem(ceres::Problem & problem,
    SfM_Data &destination_sfm_data,
    const SfM_Data &sfm_data_first,
    const SfM_Data & sfm_data_second,
    std::vector<double> & second_base_node_pose,
    double & scaling_factor,
    const std::vector<size_t> & common_track_ids
    );
};

bool MergeScenesUsingCommonTracks(SfM_Data & destination_sfm_data,
    const SfM_Data & sfm_data_first, // first submap scene
    const SfM_Data & sfm_data_second, // second submap scene
    const std::vector<size_t> & common_track_ids, // which tracks id are commonly reconstructed in both submaps
    SceneAligner *smap_aligner);

void transformSfMDataScene(SfM_Data & destination_sfm_data,
    const SfM_Data & original_sfm_data,
    const Mat3 & rotation,
    const Vec3 & translation,
    const double scaling_factor);

} // namespace sfm
} // namespace openMVG

