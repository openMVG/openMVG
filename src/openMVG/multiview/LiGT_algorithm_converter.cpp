#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>

#include "LiGT_algorithm_converter.hpp"

using namespace std;
using namespace chrono;
using namespace LiGT;

namespace LiGT {

LiGTBuilder::LiGTBuilder(const FeatsPerView& feats_per_view,
                         const PairWiseMatches& pairWise_matches,
                         const SfM_Data& sfm_data,
                         const Hash_Map<IndexT, Mat3>& map_globalR,
                         const int fixed_id){
    //normalized coordinate
    TracksOpenMVG2LiGT(feats_per_view,pairWise_matches,sfm_data);

    // load global rotations
    MapR2Attitudes(map_globalR);

    time_use_ = 0;
    fixed_id_ = fixed_id;

    output_file_ = "";
    time_file_ = "";

    // check tracks and build estimated information [EstInfo]
    CheckTracks();

}

void LiGTBuilder::MapR2Attitudes(const Hash_Map<IndexT, Mat3>& map_globalR){

    map<IndexT,Matrix3d> tmp_map_R;

    // sort by the key value (ViewId)
    for (auto& R: map_globalR){
        tmp_map_R.insert({R.first, R.second});
    }

    // generate global rotations
    global_rotations_.clear();

    for (auto& R: tmp_map_R){
        global_rotations_.emplace_back(R.second);
    }

    tmp_map_R.clear();
}

void LiGTBuilder::Landmarks2Tracks(const Landmarks& landmarks){
    PtsId pts_id = 0;
    for ( const auto& iter_landmark : landmarks) {
        const Landmark& landmark = iter_landmark.second;
        const Observations& observations = landmark.obs;
        TrackInfo track_info;
        for ( const auto& iter_obs : observations) {
            ObsInfo obs_info;
            obs_info.view_id = iter_obs.first;
            obs_info.pts_id = pts_id;

            auto& pt = iter_obs.second.x;
            obs_info.coord << pt.x(), pt.y(), 1;
            track_info.track.emplace_back(obs_info);
        }
        tracks_.emplace_back(track_info);
        pts_id++;
    }
}

void LiGTBuilder::BuildLandmarks(const FeatsPerView& feats_per_view,
                    const STLMAPTracks& map_tracks,
                    const SfM_Data& sfm_data,
                    Landmarks& landmarks){
    // Init the putative landmarks
    {
        // For every track add the obervations:
        // - views and feature positions that see this landmark
        // Fill sfm_data with the computed tracks (no 3D yet)

        IndexT idx(0);
        for (STLMAPTracks::const_iterator itTracks = map_tracks.begin();
             itTracks != map_tracks.end();
             ++itTracks, ++idx)
        {
            const submapTrack& track = itTracks->second;
            landmarks[idx] = Landmark();
            Observations& obs = landmarks.at(idx).obs;
            for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
            {
                const size_t imaIndex = it->first;
                const size_t featIndex = it->second;
                auto& pt = feats_per_view.at(imaIndex)[featIndex];

                const Vec2 x = pt.coords().cast<double>();

                const View* view = sfm_data.views.at(imaIndex).get();
                const IntrinsicBase* cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();

                const Matrix3d& Kinv = dynamic_cast<const cameras::Pinhole_Intrinsic*>(cam)->Kinv();

                // Note: LiGT needs to use the normalized feature coordinate
                Vec2 normalized_coord;
                normalized_coord(0) = Kinv(0,0) * x(0) + Kinv(0,1) * x(1) + Kinv(0,2) * 1;
                normalized_coord(1) = Kinv(1,0) * x(0) + Kinv(1,1) * x(1) + Kinv(1,2) * 1;

                obs[imaIndex] = Observation(normalized_coord, featIndex);
            }
        }
    }
}

void LiGTBuilder::BuildMapTracks(const PairWiseMatches& pairWise_matches,
                    STLMAPTracks& map_tracks){

    const openMVG::matching::PairWiseMatches& map_Matches = pairWise_matches;
    tracks::TracksBuilder tracksBuilder;

    tracksBuilder.Build(map_Matches);
    tracksBuilder.Filter();
    tracksBuilder.ExportToSTL(map_tracks);

}

void LiGTBuilder::TracksOpenMVG2LiGT(const FeatsPerView& feats_per_view,
                         const PairWiseMatches& pairWise_matches,
                         const SfM_Data& sfm_data){
    tracks::STLMAPTracks map_tracks;
    Landmarks landmarks;

    tracks_.clear();

    BuildMapTracks(pairWise_matches,map_tracks);
    BuildLandmarks(feats_per_view, map_tracks, sfm_data, landmarks);
    Landmarks2Tracks(landmarks);

    landmarks.clear();
    map_tracks.clear();
};

}
