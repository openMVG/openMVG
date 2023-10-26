// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON, Ricardo Fabbri and Gabriel Andrade.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_MATCH_CONSTRAINT_HPP
#define OPENMVG_MULTIVIEW_MATCH_CONSTRAINT_HPP

namespace openMVG
{
/// Define various enumerate for each available two view triangulation solver
enum class MultiviewMatchConstraint : unsigned char
{
  POSITION,           // only match projected/repreojected poisition for RANSAC inlier etc
  ORIENTATION,        // also match projected/repreojected orientation
  SCALE,              // also match projected/repreojected scale
  DEFAULT = POSITION
};

static inline bool isValid( const MultiviewMatchConstraint method )
{
  return
    method >= MultiviewMatchConstraint::POSITION &&
    method <= MultiviewMatchConstraint::DEFAULT;
}

}  // namespace openMVG

#endif //OPENMVG_MULTIVIEW_MATCH_CONSTRAINT_HPP
