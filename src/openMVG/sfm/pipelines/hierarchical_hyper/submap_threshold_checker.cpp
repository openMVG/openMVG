#include "submap_threshold_checker.hpp"

namespace openMVG{
namespace sfm{


SubmapThresholdChecker::SubmapThresholdChecker(int min_number_views_per_submap)
  : min_number_views_per_submap_(min_number_views_per_submap)
{}

SubmapTracksThresholdChecker::SubmapTracksThresholdChecker(int tracks_threshold)
  : SubmapThresholdChecker(2), tracks_threshold_(tracks_threshold)
{}


bool SubmapTracksThresholdChecker::operator()(const HsfmSubmap &smap) const
{
  return ((!smap.is_parent)
          && smap.sfm_data.GetViews().size() > min_number_views_per_submap_
          && (smap.track_ids.size() > tracks_threshold_));
}

bool SubmapTracksThresholdChecker::operator()(const std::pair<IndexT, HsfmSubmap> &smap_with_id) const
{
  return operator ()(smap_with_id.second);
}

SubmapViewThresholdChecker::SubmapViewThresholdChecker(int views_threshold)
  : SubmapThresholdChecker(2), views_threshold_(views_threshold)
{}

bool SubmapViewThresholdChecker::operator()(const HsfmSubmap &smap) const
{
  return ((!smap.is_parent)
          && smap.sfm_data.GetViews().size() > std::max(min_number_views_per_submap_, views_threshold_));
}

bool SubmapViewThresholdChecker::operator()(const std::pair<IndexT, HsfmSubmap> &smap_with_id) const
{
  return operator ()(smap_with_id.second);
}

} // namespace sfm
} // namespace openMVG
