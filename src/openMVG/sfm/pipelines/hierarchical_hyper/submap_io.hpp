#pragma once

#include "openMVG/sfm/sfm_data_io_cereal.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"

#include <cereal/types/set.hpp>
#include <cereal/types/memory.hpp>

namespace openMVG {
namespace sfm {

/** Serialization
* @param ar Archive
*/
template <class Archive>
inline void HsfmSubmap::serialize ( Archive & ar )
{
  ar(cereal::make_nvp("track_ids", track_ids),
     cereal::make_nvp("is_parent", is_parent),
     cereal::make_nvp("separator", separator),
     cereal::make_nvp("children_submaps_first", children_submaps.first),
     cereal::make_nvp("children_submaps_second", children_submaps.second),
     cereal::make_nvp("parent_id", parent_id),
     cereal::make_nvp("sfm_data", sfm_data));
}

} // namespace sfm
} // namespace openMVG
