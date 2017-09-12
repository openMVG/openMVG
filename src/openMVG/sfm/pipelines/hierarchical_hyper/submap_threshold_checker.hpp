#ifndef SUBMAPPARTITIONABLEPREDICATE_HPP
#define SUBMAPPARTITIONABLEPREDICATE_HPP

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"

namespace openMVG{
namespace sfm{

class SubmapThresholdChecker
{
public:
  SubmapThresholdChecker(int min_number_views_per_submap);

  virtual bool operator()(const HsfmSubmap & smap) const = 0;

  virtual bool operator()(const std::pair<openMVG::IndexT, HsfmSubmap> & smap) const = 0;

protected:

  int min_number_views_per_submap_;

};

class SubmapTracksThresholdChecker : public SubmapThresholdChecker
{
public:
  SubmapTracksThresholdChecker(int tracks_threshold);

  bool operator()(const HsfmSubmap & smap) const;
  bool operator()(const std::pair<openMVG::IndexT, HsfmSubmap> & smap_with_id) const;

private:

  int tracks_threshold_;

};

class SubmapViewThresholdChecker : public SubmapThresholdChecker
{
public:
  SubmapViewThresholdChecker(int views_threshold);

  bool operator()(const HsfmSubmap & smap) const;
  bool operator()(const std::pair<openMVG::IndexT, HsfmSubmap> & smap_with_id) const;

private:

  int views_threshold_;

};

}// namespace sfm
}// namespace openMVG

#endif // SUBMAPPARTITIONABLEPREDICATE_HPP
