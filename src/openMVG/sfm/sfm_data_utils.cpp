// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_utils.hpp"
#include <unordered_map>
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/sfm/sfm_data.hpp"

namespace {

using namespace openMVG;
using namespace openMVG::sfm;

inline bool SortViewFunc(const std::string &s1, const std::string &s2)
{
  return s1.compare(s2) < 0;
}

std::unordered_map<IndexT, IndexT> ComputeReorderMap(const std::set<IndexT> &v)
{
  std::vector<IndexT> w;
  std::copy(v.begin(), v.end(), std::back_inserter(w));
  std::sort(w.begin(), w.end());

  std::unordered_map<IndexT, IndexT> res;
  for (IndexT idx = 0; idx < w.size(); ++idx)
  {
    res[w[idx]] = idx;
  }

  return res;
}

}

namespace openMVG {
namespace sfm {

void GroupSharedIntrinsics(SfM_Data & sfm_data)
{
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  // Build hash & build a set of the hash in order to maintain unique Ids
  std::set<size_t> hash_index;
  std::vector<size_t> hash_value;

  for (const auto & intrinsic_it : intrinsics)
  {
    const cameras::IntrinsicBase * intrinsicData = intrinsic_it.second.get();
    const size_t hashVal = intrinsicData->hashValue();
    hash_index.insert(hashVal);
    hash_value.push_back(hashVal);
  }

  // From hash_value(s) compute the new index (old to new indexing)
  Hash_Map<IndexT, IndexT> old_new_reindex;
  size_t i = 0;
  for (const auto & intrinsic_it : intrinsics)
  {
    old_new_reindex[intrinsic_it.first] = std::distance(hash_index.cbegin(), hash_index.find(hash_value[i]));
    ++i;
  }
  //--> Save only the required Intrinsics (do not need to keep all the copy)
  Intrinsics intrinsic_updated;
  for (const auto & intrinsic_it : intrinsics)
  {
    intrinsic_updated[old_new_reindex[intrinsic_it.first]] = intrinsics[intrinsic_it.first];
  }
  // Update intrinsics (keep only the necessary ones) -> swapping
  intrinsics.swap(intrinsic_updated);

  // Update views intrinsic IDs (since some intrinsic position have changed in the map)
  for (auto & view_it: views)
  {
    View * v = view_it.second.get();
    // Update the Id only if a corresponding index exists
    if (old_new_reindex.count(v->id_intrinsic))
      v->id_intrinsic = old_new_reindex[v->id_intrinsic];
  }
}

void SortAndCleanSfMData(SfM_Data &sfm_data)
{
  std::vector<std::string> view_names;
  for (auto view_iter : sfm_data.GetViews())
  {
    const View * view = view_iter.second.get();
    if (!sfm_data.IsIntrinsicDefined(view))
    {
      continue;
    }
    view_names.push_back(view->s_Img_path);
  }

  SortAndCleanSfMData(sfm_data, view_names);
}

void SortAndCleanSfMData(SfM_Data &sfm_data, std::vector<std::string> &view_names)
{
  std::sort(view_names.begin(), view_names.end(), SortViewFunc);

  std::unordered_map<std::string, IndexT> valid_name_id_map;
  std::set<IndexT> valid_intrinsic_ids;
  std::unordered_map<IndexT, IndexT> intrinsic_map;
  std::unordered_map<IndexT, IndexT> view_map;

  for (auto view_iter : sfm_data.GetViews())
  {
    const View * view = view_iter.second.get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
    {
      continue;
    }

    valid_name_id_map[view->s_Img_path] = view->id_view;
    valid_intrinsic_ids.insert(view->id_intrinsic);
  }
  intrinsic_map = ComputeReorderMap(valid_intrinsic_ids);

  SfM_Data sfm_data_new;
  sfm_data_new.s_root_path = sfm_data.s_root_path;

  for (IndexT idx = 0; idx < view_names.size(); ++idx)
  {
    auto iter = valid_name_id_map.find(view_names[idx]);
    if (iter == valid_name_id_map.end())
    {
      continue;
    }

    auto view = sfm_data.views[iter->second];
    view->id_view = idx;
    view->id_pose = idx;
    view->id_intrinsic = intrinsic_map[view->id_intrinsic];
    sfm_data_new.views[idx] = view;
    sfm_data_new.poses[idx] = sfm_data.poses[view->id_pose];
    view_map[iter->second] = idx;
  }

  for (const auto & map : intrinsic_map)
  {
    sfm_data_new.intrinsics[map.second] = sfm_data.intrinsics[map.first];
  }

  for (const auto & landmark : sfm_data.structure)
  {
    Observations obs_new;
    for (const auto & observation : landmark.second.obs)
    {
      auto view_map_iter = view_map.find(observation.first);
      if (view_map_iter == view_map.end())
      {
        continue;
      }
      auto view_id_new = view_map_iter->second;
      std::pair<IndexT, Observation> ob_new(view_id_new, observation.second);
      obs_new.insert(ob_new);
    }
    if (obs_new.empty())
    {
      continue;
    }
    sfm_data_new.structure[landmark.first].obs = std::move(obs_new);
    sfm_data_new.structure[landmark.first].X = std::move(landmark.second.X);
  }

  std::swap(sfm_data, sfm_data_new);
}

} // namespace sfm
} // namespace openMVG
