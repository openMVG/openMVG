// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SUBMAPPARTITIONABLEPREDICATE_HPP
#define SUBMAPPARTITIONABLEPREDICATE_HPP

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"

namespace openMVG{
namespace sfm{

/**
 * @brief A base class for a functor that checks if a submap is partitionable,
 * given some threshold.
 */
class SubmapThresholdChecker
{
public:
  SubmapThresholdChecker() = default;

  virtual bool operator()(const HsfmSubmap & smap) const = 0;
  virtual bool operator()(const std::pair<openMVG::IndexT, HsfmSubmap> & smap) const = 0;

};

/**
 * @brief A functor used to determine if a submap should be partitioned based on a TRACKS threshold
 */
class SubmapTracksThresholdChecker : public SubmapThresholdChecker
{
public:
  explicit SubmapTracksThresholdChecker(int tracks_threshold);

  bool operator()(const HsfmSubmap & smap) const;
  bool operator()(const std::pair<openMVG::IndexT, HsfmSubmap> & smap_with_id) const;

private:

  int tracks_threshold_;

};

/**
 * @brief A functor used to determine if a submap should be partitioned based on a VIEW threshold
 */
class SubmapViewThresholdChecker : public SubmapThresholdChecker
{
public:
  explicit SubmapViewThresholdChecker(int views_threshold);

  bool operator()(const HsfmSubmap & smap) const;
  bool operator()(const std::pair<openMVG::IndexT, HsfmSubmap> & smap_with_id) const;

private:

  int views_threshold_;

};

}// namespace sfm
}// namespace openMVG

#endif // SUBMAPPARTITIONABLEPREDICATE_HPP
