// This file is part of a pose-only algorithm of Linear Global Translation (LiGT)

// Copyright (c) 2022, Qi Cai and Yuanxin Wu

// This Source Code Form is subject to the license terms of
// Creative Commons Attribution Share Alike 4.0 International.
// Details are available at https://choosealicense.com/licenses/cc-by-sa-4.0/


#ifndef LIGT_TYPES
#define LIGT_TYPES

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/types.hpp"

#include <set>
#include <unordered_map>
#include <vector>

using namespace openMVG;

namespace LiGT{

using ViewId = IndexT;
using PtsId = IndexT;
using ObsId = IndexT;

// pose storage
using Rotations = std::vector<Mat3>;
using Translations = std::vector<Vec3>;

using Poses = std::unordered_map<ViewId, openMVG::geometry::Pose3>;

// image observation information
struct ObsInfo {
  ViewId view_id;
  PtsId pts_id;
  Eigen::Vector3d coord;
};

// track information
using Track = std::vector<ObsInfo>;

struct TrackInfo {
  Track track;
  bool is_used = true;
};

using Tracks = std::vector<TrackInfo>;

// estimated information
using EstimatedViewIds = std::vector<ViewId>;
using OriginViewIds = std::vector<ViewId>;
using Origin2EstViewIds = std::unordered_map<ViewId, ViewId>;
using Est2OriginViewIds = std::unordered_map<ViewId, ViewId>;

struct EstInfo {
  EstimatedViewIds estimated_view_ids;
  OriginViewIds origin_view_ids;
  Origin2EstViewIds origin2est_view_ids;
  Est2OriginViewIds est2origin_view_ids;

  void BuildMap(){
    origin2est_view_ids.clear();
    est2origin_view_ids.clear();

    if (estimated_view_ids.size() > 0
        && origin_view_ids.size() > 0
        && estimated_view_ids.size() == origin_view_ids.size())
    {
      for ( ViewId i = 0; i < estimated_view_ids.size(); ++i){
        origin2est_view_ids.insert({origin_view_ids[i],estimated_view_ids[i]});
        est2origin_view_ids.insert({estimated_view_ids[i],origin_view_ids[i]});
      }
    }
  }
};

}

#endif
