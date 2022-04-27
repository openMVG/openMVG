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
    BuildTracks(feats_per_view, pairWise_matches, sfm_data);

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

void LiGTBuilder::BuildTracks(const FeatsPerView& feats_per_view,
                         const PairWiseMatches& pairWise_matches,
                         const SfM_Data& sfm_data){
    // Build OpenMVG tracks
    tracks::STLMAPTracks map_tracks;
    {
      const openMVG::matching::PairWiseMatches& map_Matches = pairWise_matches;
      tracks::TracksBuilder tracksBuilder;

      tracksBuilder.Build(map_Matches);
      tracksBuilder.Filter();
      tracksBuilder.ExportToSTL(map_tracks);
    }

    tracks_.clear();

    // For every track add the observations to LiGT Tracks
    // - views and feature positions that see this landmark

    tracks_.reserve(map_tracks.size());
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
#endif
    for (IndexT trackId = 0; trackId < map_tracks.size(); ++trackId)
    {
      openMVG::tracks::STLMAPTracks::const_iterator itTracks = map_tracks.begin();
      std::advance(itTracks, trackId);
      TrackInfo track_info;
      const submapTrack& track = itTracks->second;
      for (const auto& it : track)
      {
          const size_t imaIndex = it.first;
          const size_t featIndex = it.second;
          const auto& pt = feats_per_view.at(imaIndex)[featIndex];

          const View* view = sfm_data.views.at(imaIndex).get();
          const IntrinsicBase* cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();

          const Vec2 x = cam->get_ud_pixel(pt.coords().cast<double>());

          ObsInfo obs_info;
          obs_info.view_id = imaIndex;
          obs_info.pts_id = trackId;
          obs_info.coord = (*cam)(x).col(0);
          track_info.track.emplace_back(obs_info);
      }
      tracks_.emplace_back(track_info);
    }

    map_tracks.clear();
};

}
