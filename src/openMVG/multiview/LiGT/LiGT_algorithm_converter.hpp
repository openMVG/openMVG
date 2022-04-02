
#ifndef LIGT_ALGORITHM_CONVERTER
#define LIGT_ALGORITHM_CONVERTER

#pragma once

#include "LiGT_algorithm.hpp"

namespace LiGT {

using IndexT = uint32_t;
typedef Hash_Map<IndexT, features::PointFeatures> FeatsPerView;

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
    LiGTBuilder(const FeatsPerView& feats_per_view,
                const PairWiseMatches& pairWise_matches,
                const SfM_Data& sfm_data,
                const Hash_Map<IndexT, Mat3>& map_globalR,
                const int fixed_id = 0);

    virtual ~LiGTBuilder() = default;

    // transform openMVG's [Landmarks] into [Tracks] in LiGT
    void Landmarks2Tracks(const Landmarks& landmarks);

    // build [Landmarks] in openMVG
    void BuildLandmarks(const FeatsPerView& feats_per_view,
                        const STLMAPTracks& map_tracks,
                        const SfM_Data& sfm_data,
                        Landmarks& landmarks);

    // build [MapTracks] in openMVG
    void BuildMapTracks(const PairWiseMatches& pairWise_matches,
                        STLMAPTracks& map_tracks);

    // transform openMVG's tracks info into the LiGT's form
    void TracksOpenMVG2LiGT(const FeatsPerView& feats_per_view,
                             const PairWiseMatches& pairWise_matches,
                             const SfM_Data& sfm_data);


    // transform openMVG's map global rotations into the LiGT's attitude form
    void MapR2Attitudes(const Hash_Map<IndexT, Mat3>& map_globalR);
};
}

#endif
