// This file is part of a pose-only algorithm of Linear Global Translation (LiGT)

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the license terms of
// Creative Commons Attribution Share Alike 4.0 International.
// Details are available at https://choosealicense.com/licenses/cc-by-sa-4.0/

#ifndef LIGT_ALGORITHM_CONVERTER
#define LIGT_ALGORITHM_CONVERTER

#pragma once

#include "LiGT_algorithm.hpp"

#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/tracks/tracks.hpp"

using namespace openMVG;

namespace LiGT {

//------------------
//-- Bibliography --
//------------------
//- [1] "A Pose-only Solution to Visual Reconstruction and Navigation".
//- Authors: Qi Cai, Lilian Zhang, Yuanxin Wu, Wenxian Yu and Dewen Hu.
//- Date: December 2022.
//- Journal: IEEE T-PAMI.
//
//- [2] "Equivalent constraints for two-view geometry: pose solution/pure rotation identification and 3D reconstruction".
//- Authors: Qi Cai; Yuanxin Wu; Lilian Zhang; Peike Zhang.
//- Date: December 2019.
//- Journal: IJCV.


class LiGTBuilder : public LiGTProblem{
  // Readme. For this version, we provide following converter functions
  // to transform the data type in openMVG into LiGT.
  // In the future, we might directly use openMVG's data type in the LiGT algorithm.
public:
  LiGTBuilder(const sfm::Features_Provider* features_provider,
        const matching::PairWiseMatches& pairWise_matches,
        const sfm::SfM_Data& sfm_data,
        const Hash_Map<IndexT, Mat3>& map_globalR,
        const int min_track_length = 2,
        const int fixed_id = 0);

  virtual ~LiGTBuilder() = default;

  // transform openMVG's tracks info into the LiGT's form
  void BuildTracks(const sfm::Features_Provider* features_provider,
           const matching::PairWiseMatches& pairWise_matches,
           const sfm::SfM_Data& sfm_data);


  // Initialize Global rotations
  void MapR2Rotations(const Hash_Map<IndexT, Mat3>& map_globalR);
};
}

#endif
