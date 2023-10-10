#include "openMVG/sfm/pipelines/sequential/sequential_SfM_util.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/stl/stl.hpp"

namespace openMVG {
namespace sfm {

using namespace openMVG::matching;

// TODO(trifocal future) GetTripletWithMostMatches
// Get the PairWiseMatches that have the most support point
std::vector<openMVG::matching::PairWiseMatches::const_iterator>
GetPairWithMostMatches(const SfM_Data& sfm_data, const PairWiseMatches& matches, int clamp_count) {

  std::vector<openMVG::matching::PairWiseMatches::const_iterator> sorted_pairwise_matches_iterators;
  // List Views that supports valid intrinsic
  std::set<IndexT> valid_views;
  for (const auto & view : sfm_data.GetViews())
  {
    const View * v = view.second.get();
    if (sfm_data.GetIntrinsics().find(v->id_intrinsic) != sfm_data.GetIntrinsics().end())
      valid_views.insert(v->id_view);
  }

  if (sfm_data.GetIntrinsics().empty() || valid_views.empty())
  {
    OPENMVG_LOG_ERROR
      << "Unable to choose an initial pair, since there is no defined intrinsic data.";
    return {};
  }

  // Try to list the clamp_count top pairs that have valid intrinsics
  std::vector<uint32_t > vec_NbMatchesPerPair;
  std::vector<openMVG::matching::PairWiseMatches::const_iterator> vec_MatchesIterator;
  const openMVG::matching::PairWiseMatches & map_Matches = matches;
  for (openMVG::matching::PairWiseMatches::const_iterator
    iter = map_Matches.begin();
    iter != map_Matches.end(); ++iter)
  {
    const Pair current_pair = iter->first;
    if (valid_views.count(current_pair.first) &&
      valid_views.count(current_pair.second) )
    {
      vec_NbMatchesPerPair.push_back(iter->second.size());
      vec_MatchesIterator.push_back(iter);
    }
  }
  // sort the Pairs in descending order according their correspondences count
  using namespace stl::indexed_sort;
  std::vector<sort_index_packet_descend<uint32_t, uint32_t>> packet_vec(vec_NbMatchesPerPair.size());
  sort_index_helper(packet_vec, &vec_NbMatchesPerPair[0], std::min((size_t)clamp_count, vec_NbMatchesPerPair.size()));

  for (size_t i = 0; i < std::min((size_t)clamp_count, vec_NbMatchesPerPair.size()); ++i) {
    const uint32_t index = packet_vec[i].index;
    sorted_pairwise_matches_iterators.emplace_back(vec_MatchesIterator[index]);
  }
  return sorted_pairwise_matches_iterators;
}

} // namespace sfm
} // namespace openMVG
