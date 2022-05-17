// This file is part of a pose-only algorithm of Linear Global Translation (LiGT)

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the license terms of
// Creative Commons Attribution Share Alike 4.0 International.
// Details are available at https://choosealicense.com/licenses/cc-by-sa-4.0/

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <chrono>
#include <map>

#include "LiGT_algorithm_converter.hpp"

#include "openMVG/cameras/cameras.hpp"

using namespace Eigen;
using namespace std::chrono;

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::tracks;
using namespace openMVG::cameras;

namespace LiGT {

LiGTBuilder::LiGTBuilder(const Features_Provider* features_provider,
             const PairWiseMatches& pairWise_matches,
             const SfM_Data& sfm_data,
             const Hash_Map<IndexT, Mat3>& map_globalR,
             const int min_track_length,
             const int fixed_id){
  // build observation and bearing vectors
  BuildTracks(features_provider, pairWise_matches, sfm_data);

  // load global rotations
  MapR2Rotations(map_globalR);

  time_use_ = 0;
  fixed_id_ = fixed_id;
  min_track_length_ = min_track_length;

  // check tracks and build estimated information [EstInfo]
  CheckTracks();
}

void LiGTBuilder::MapR2Rotations(const Hash_Map<IndexT, Mat3>& map_globalR){
  std::map<IndexT,Matrix3d> tmp_map_R;

  // sort by the key value (ViewId)
  for (auto& R: map_globalR){
    tmp_map_R.insert({R.first, R.second});
  }

  // generate global rotations
  global_rotations_.clear();

  for (auto& R: tmp_map_R){
    global_rotations_.emplace_back(R.second);
  }
}

void LiGTBuilder::BuildTracks(const Features_Provider* features_provider,
             const PairWiseMatches& pairWise_matches,
             const SfM_Data& sfm_data){
  // Build OpenMVG tracks
  tracks::STLMAPTracks map_tracks;
  {
    tracks::TracksBuilder tracksBuilder;

    tracksBuilder.Build(pairWise_matches);
    tracksBuilder.Filter(min_track_length_);
    tracksBuilder.ExportToSTL(map_tracks);
  }

  tracks_.clear();

  // For every track add the observations to LiGT Tracks
  // - views and feature positions that see this landmark

  tracks_.reserve(map_tracks.size());
  IndexT trackId(0);
  for (const auto& trackIt : map_tracks)
  {
    TrackInfo track_info;
    const submapTrack& track = trackIt.second;
    track_info.track.reserve(track.size());
    for (const auto& it : track)
    {
      const size_t imaIndex = it.first;
      const size_t featIndex = it.second;
      const auto& pt = features_provider->feats_per_view.at(imaIndex)[featIndex];

      const View* view = sfm_data.views.at(imaIndex).get();
      const IntrinsicBase* cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();

      const Vec2 x = cam->get_ud_pixel(pt.coords().cast<double>());

      ObsInfo obs_info;
      obs_info.view_id = imaIndex;
      obs_info.pts_id = trackId;
      obs_info.coord = (*cam)(x);
      track_info.track.emplace_back(obs_info);
    }
    tracks_.emplace_back(track_info);
    ++trackId;
  }
};

}
