// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_TYPE_HPP
#define OPENMVG_MULTIVIEW_TRIANGULATION_TYPE_HPP

namespace openMVG
{
/// Define various enumerate for each available two view triangulation solver
enum class ETriangulationMethod : unsigned char
{
  DIRECT_LINEAR_TRANSFORM, // DLT
  L1_ANGULAR,
  LINFINITY_ANGULAR,
  INVERSE_DEPTH_WEIGHTED_MIDPOINT,
  DEFAULT = INVERSE_DEPTH_WEIGHTED_MIDPOINT
};

static inline bool isValid( const ETriangulationMethod method )
{
  return
    method >= ETriangulationMethod::DIRECT_LINEAR_TRANSFORM &&
    method <= ETriangulationMethod::DEFAULT;
}

}  // namespace openMVG

#endif //OPENMVG_MULTIVIEW_TRIANGULATION_TYPE_HPP
