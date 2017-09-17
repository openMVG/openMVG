// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef CONTROLPOINTTABLEVIEW_HPP
#define CONTROLPOINTTABLEVIEW_HPP

#include <QDialog>
#include <QTableWidget>
#include "openMVG/sfm/sfm_data.hpp"

namespace control_point_GUI {

/// QT Interface to edit Landmarks GCP data:
/// Allow to edit X,Y,Z Ground Control Point coordinates
/// Allow to delete GCP
class ControlPointTableView : public QDialog
{
public:
  ControlPointTableView
  (
    const openMVG::sfm::SfM_Data * sfm_data,
    QWidget *parent = nullptr
  );

  /// Update control points X,Y,Z data (if valid datum is provided)
  void updateControlPoints(openMVG::sfm::Landmarks & control_points);

  /// Delete selected control_points row(s) on Key_Delete event
  void keyReleaseEvent(QKeyEvent* event);

private:
  // Input SfM scene data & Control points data
  const openMVG::sfm::SfM_Data * sfm_data_;
  // Copy of sfm_data_ control_points (work on a copy, to do not break existing ones)
  openMVG::sfm::Landmarks control_points_;
  /// The graphic interface allowing to delete and edit control points data
  QTableWidget * table_;
};

} // namespace control_point_GUI

#endif /* CONTROLPOINTTABLEVIEW_HPP */
