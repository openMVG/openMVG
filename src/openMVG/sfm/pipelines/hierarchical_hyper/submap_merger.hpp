// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp"

namespace openMVG{
namespace sfm{

class SceneAligner;

/**
 * @brief This class is used to merge a set of hypersfm submaps back together
 */
class SubmapMerger
{
public:
  explicit SubmapMerger(const HsfmSubmaps &submaps);

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
