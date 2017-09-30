// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ALTERNATIVE_VO_VO_VIEWER_DRAWABLE_ELEMENTS_HPP_
#define ALTERNATIVE_VO_VO_VIEWER_DRAWABLE_ELEMENTS_HPP_

#include "openMVG/numeric/numeric.h"

#include <QColor>

#include <vector>

namespace alternative_vo
{
/**
* @brief Helper structure use to represent a line and it's associated color
* @brief a line is composed of a list of points
*/
struct VOViewerLine
{
  /// Positions of the points
  std::vector< openMVG::Vec2f > m_pts;

  /// Color of the elements
  QColor m_color;

  /// Default color for a line
  static const QColor DEFAULT_LINE_COLOR;
};

/**
* Helper structure used to represent a point and it's associated color
*/
struct VOViewerPoint
{
  /// Position of the point
  openMVG::Vec2f m_pt;

  /// Color of the point
  QColor m_color;

  /// Default color for a tracked point
  static const QColor DEFAULT_TRACKED_POINT_COLOR;

  /// Default color for a newly added point
  static const QColor DEFAULT_NEW_POINT_COLOR;
};
}

#endif
