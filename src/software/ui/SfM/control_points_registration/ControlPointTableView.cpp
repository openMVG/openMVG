// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ControlPointTableView.hpp"
#include <QDialog>
#include <QSet>
#include <QVBoxLayout>
#include <QWidget>
#include <QTableWidget>
#include <QMessageBox>
#include <QHeaderView>
#include <QKeyEvent>

namespace control_point_GUI {

using namespace openMVG::sfm;

ControlPointTableView::ControlPointTableView
(
  const SfM_Data * sfm_data,
  QWidget *parent
)
: QDialog(parent),
  sfm_data_(sfm_data)
{
  if (sfm_data_ == nullptr)
    return;
  if (sfm_data_->GetControl_Points().empty())
    return;

  //--
  // Init a TableWidget from control points data
  //--

  const int row_count = sfm_data_->GetControl_Points().size();
  const int col_count = 5; // {Id, #obs, X, Y, Z}

  table_ = new QTableWidget(row_count, col_count);
  table_->verticalHeader()->hide();
  //QDialog * control_point_coords_gui = new QDialog(parent);
  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(table_);
  //control_point_coords_gui->setLayout(layout);
  this->setLayout(layout);

  table_->setHorizontalHeaderItem(0, new QTableWidgetItem("Id"));
  table_->setHorizontalHeaderItem(1, new QTableWidgetItem("#obs"));
  table_->setHorizontalHeaderItem(2, new QTableWidgetItem("X"));
  table_->setHorizontalHeaderItem(3, new QTableWidgetItem("Y"));
  table_->setHorizontalHeaderItem(4, new QTableWidgetItem("Z"));

  int i = 0;
  for (Landmarks::const_iterator iterL = sfm_data_->GetControl_Points().begin();
    iterL != sfm_data_->GetControl_Points().end(); ++iterL, ++i)
  {
    // ID, #obs, X, Y, Z
    QTableWidgetItem * item;

    // ID
    item = new QTableWidgetItem( QString::number(iterL->first) );
    item->setFlags(item->flags() & ~Qt::ItemIsEditable); // not editable item
    table_->setItem(i,0,item);
    // #obs
    item = new QTableWidgetItem( QString::number(iterL->second.obs.size()) );
    item->setFlags(item->flags() & ~Qt::ItemIsEditable); // not editable item
    table_->setItem(i,1,item);
    // X
    item = new QTableWidgetItem( QString::number(iterL->second.X(0)) );
    table_->setItem(i,2,item);
    // Y
    item = new QTableWidgetItem( QString::number(iterL->second.X(1)) );
    table_->setItem(i,3,item);
    // Z
    item = new QTableWidgetItem( QString::number(iterL->second.X(2)) );
    table_->setItem(i,4,item);
  }
  // resize the window according to the table
  table_->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
  this->setMinimumSize(table_->horizontalHeader()->length(), table_->verticalHeader()->length());
  // Work internally on a copy (in order to don't break existing keypoint)
  control_points_ = sfm_data_->GetControl_Points();
}

/// Update control points X,Y,Z data (if valid datum is provided)
void ControlPointTableView::updateControlPoints
(
  Landmarks & control_points
)
{
  // Get back the (updated) control points coordinates
  for (int i = 0; i < table_->rowCount(); ++i)
  {
    bool bOk_X = false, bOk_Y = false, bOk_Z = false;

    const int index = table_->item(i, 0)->text().toInt(&bOk_X);

    const double
      valueX = table_->item(i, 2)->text().toFloat(&bOk_X),
      valueY = table_->item(i, 3)->text().toFloat(&bOk_Y),
      valueZ = table_->item(i, 4)->text().toFloat(&bOk_Z);

    if (bOk_X && bOk_Y && bOk_Z)
    {
      control_points_[index].X << valueX, valueY, valueZ;
    }
    else
    {
      QMessageBox msgBox;
      msgBox.setText("Floating point values are required for the GCP coordinates.");
      msgBox.exec();
      return;
    }
  }
  control_points = control_points_;
}

/// Delete selected control_points row(s) on Key_Delete event
void ControlPointTableView::keyReleaseEvent
(
  QKeyEvent* event
)
{
  if (event->key() == Qt::Key_Delete)
  {
    QItemSelection selection( table_->selectionModel()->selection() );

    QSet<int> selected_rows_ids;
    foreach ( const QModelIndex & index, selection.indexes() )
    {
      selected_rows_ids.insert( index.row() );
    }

    int i = 0;
    for ( const int & index :  selected_rows_ids )
    {
      // Retrieve the control point Ids
      bool bOk_X = false;
      const int cp_index = table_->item(index - i, 0)->text().toInt(&bOk_X);
      // Remove the control points entry and it's GUI representation
      control_points_.erase(cp_index);
      table_->removeRow(index - i);
      ++i;
    }
  }
  QDialog::keyReleaseEvent(event);
}

} // namespace control_point_GUI
