#pragma once

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io_cereal.hpp"

#include <cereal/types/set.hpp>
#include <cereal/types/memory.hpp>

#include <set>

namespace openMVG {
namespace sfm {

// struct used to store the submaps created by clustering.
// A HsfmSubmap can be a parent, if it has been clustered into two
// children submaps.
// It stores a sfm_data, and the tracks ids corresponding to the landmarks contained in the submap.
// If it is a parent it also stores the ids of the separator tracks (which are tracks contained in
// both children submaps), and the submap id of its children.
struct HsfmSubmap
{
  SfM_Data sfm_data;
  std::set<size_t> track_ids;
  bool is_parent = false;
  std::set<size_t> separator;// set of track ids in the separator (if parent)
  std::pair<IndexT, IndexT> children_submaps = {0,0}; // pair of submaps index for children (if any)
  IndexT parent_id = 0;

  /** Serialization
  * @param ar Archive
  */
  template <class Archive>
  void serialize ( Archive & ar )
  {
    ar(cereal::make_nvp("track_ids", track_ids),
       cereal::make_nvp("is_parent", is_parent),
       cereal::make_nvp("separator", separator),
       cereal::make_nvp("children_submaps_first", children_submaps.first),
       cereal::make_nvp("children_submaps_second", children_submaps.second),
       cereal::make_nvp("parent_id", parent_id),
       cereal::make_nvp("sfm_data", sfm_data));
  }
};

using HsfmSubmaps = std::map<openMVG::IndexT, openMVG::sfm::HsfmSubmap>;

// export submap data as sfm_data.json and ply (with separator tracks highlighted)
// in output directory
bool ExportSubmapData(const HsfmSubmaps & submaps, const IndexT submap_id, const std::string & file_path);

} // namespace sfm
} // namespace openMVG
