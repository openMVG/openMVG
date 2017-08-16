#pragma once

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp"

namespace openMVG{
namespace sfm{

class SceneAligner;

class SubmapMerger
{
public:
  SubmapMerger(const HsfmSubmaps &submaps);

  bool Merge(const std::string &debug_out_dir = "");

  HsfmSubmaps getSubmaps() const {return submaps_;}

protected:
  bool MergeSubmapPair(const IndexT parent_id, const std::string &debug_out_dir);
  virtual bool BundleAdjustment_FixedSeparators(const IndexT parent_id);

  std::unique_ptr<SceneAligner> scene_aligner_;
  HsfmSubmaps submaps_;

};

}
}
