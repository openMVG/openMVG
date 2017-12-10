// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QtGui>
#include <QTreeView>
#include <string>

#include "GraphicsView.hpp"
#include "document.hpp"

class MainWindow : public QMainWindow
{
  Q_OBJECT

private:
  QTabWidget * m_tabWidget;
  QWidget * m_tab_1;
  QTreeView * m_treeView_Images;   // Image list of the project

  QStatusBar *m_statusbar;         // Status bar
  QMenuBar *m_menubar;             // Menu bar

  QMenu *m_menuHelp;               // Menu
  QMenu *m_menuFile;               // Help

  //Actions
  QAction * m_open_action;
  QAction * m_save_action;
  QAction * m_register_action;
  QAction * m_edit_cp_action;
  QAction * m_delete_control_point_action;
  QAction * m_help_action;

  // -- DOCUMENT
  Document m_doc;
  std::string m_sfm_data_filename;
  // -- END DOCUMENT

  // -- VIEW
  control_point_GUI::GraphicsView * m_widget;
  // -- END VIEW

  private slots:

  /// Delete all the SfM_Data scene control points
  void removeAllControlPoints();

  /// Display a brief help
  void help();

  /// handle the event Double Click on the SfM_Data image view list
  /// - Display the chosen image & display view's control point observations
  void doubleClickImageList();

  /// Save the SfM_Data scene (perhaps updated with GCP data)
  void saveProject();

  /// Open a SfM_Data scene
  void openProject();

  /// Display the GUI that allow to edit 3D GCP positions
  void editControlPoints();

  /// Perform the registration of the SfM_Data scene to the GCP
  void registerProject();

public:

  /// Constructor
  MainWindow(QWidget * parent = 0);

  /* Create the Window layout (shape panels & co.)
  |____________________________
  |Tab1     |                  |
  |---------|                  |
  |         |     Drawing      |
  | Img     |      Widget      |
  | list    |                  |
  |         |                  |
  |_________|__________________|
*/
  void createPanel();

  /// Create windows menus
  void createMenus();

  /// Update actions (related to menu)
  void createActions();

  /// Establish action-event connections
  void createConnections();
};

#endif /* MAINWINDOW_HPP */
