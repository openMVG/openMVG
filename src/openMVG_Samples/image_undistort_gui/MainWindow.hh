// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef _OPENMVG_SAMPLE_IMAGE_UNDISTORT_GUI_MAIN_WINDOW_HH_
#define _OPENMVG_SAMPLE_IMAGE_UNDISTORT_GUI_MAIN_WINDOW_HH_

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"

#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGraphicsView>
#include <QLabel>
#include <QMainWindow>
#include <QMenu>
#include <QPushButton>

namespace image_undistort_gui
{

/**
  * @brief Main class of the GUI interface 
  */
class MainWindow : public QMainWindow
{
public:
  /**
  * @brief Constructor
  * @param parent Parent widget of this window 
  */
  MainWindow( QWidget *parent = nullptr );

public slots:

  /**
  * @brief Action to be executed when user want to open an image 
  */
  void onOpenImage( void );

  /**
  * @brief Action to be executed when user want to undistort an image 
  */
  void onProcessImage( void );

  /**
  * @brief Action to be executed when user want to export an undistorted image 
  */
  void onExportImage( void );

  /**
  * @brief Action to be executed when user select an intrinsic type 
  */
  void onSelectIntrinsicType( void );

  /**
  * @brief Action to be executed when user check/uncheck the grid checkbox
  */
  void onChangeCheckBox( void );

  /**
  * @brief action to be executed when user select no distorsion menu item
  */
  void onMenuSelectNoDist( void );

  /**
  * @brief action to be executed when user select radial 1 menu item
  */
  void onMenuSelectRadial1( void );

  /**
  * @brief action to be executed when user select radial 3 menu item
  */
  void onMenuSelectRadial3( void );

  /**
  * @brief action to be executed when user select brown menu item
  */
  void onMenuSelectBrown( void );

  /**
  * @brief action to be executed when user select fisheye menu item
  */
  void onMenuSelectFisheye( void );

  /**
  * @brief action to be executed when user select quit menu item
  */
  void onQuit( void );

private:
  static const int STATE_INITIAL;
  static const int STATE_HAS_INPUT;
  static const int STATE_HAS_COMPUTED_DISTORTION;

  int m_current_state;

  /**
  * @brief Fill intrinsic values in the combobox 
  */
  void FillIntrinsicTypes( void );

  /**
  * @brief Build virtual image representing a regular grid 
  */
  void BuildGrid( void );

  /**
  * @brief update interface based on internal state 
  */
  void UpdateInterface( void );

  /**
  * @brief Build Interface widgets 
  */
  void BuildInterface();

  /**
  * @brief Create menus items 
  */
  void BuildMenus();

  /**
  * @brief Make connections between interface items and slots 
  */
  void MakeConnections();

  QPushButton *m_open_image;
  QCheckBox *m_check_use_grid; // use a grid instead of an image

  QLabel *m_intrinsic_label;
  QComboBox *m_intrinsic_type; // Pinhole | pinhole radial 1 | pinhole radial 2 | pinhole radial 3 | pinhole brown | pinhole fisheye

  // Focal
  QLabel *m_label_fov;
  QDoubleSpinBox *m_fov_spin;

  // Center point
  QLabel *m_label_ppx;
  QDoubleSpinBox *m_ppx;
  QLabel *m_label_ppy;
  QDoubleSpinBox *m_ppy;

  // Radial/fisheye
  QLabel *m_label_k1;
  QDoubleSpinBox *m_k1;
  QLabel *m_label_k2;
  QDoubleSpinBox *m_k2;
  QLabel *m_label_k3;
  QDoubleSpinBox *m_k3;
  QLabel *m_label_k4;
  QDoubleSpinBox *m_k4;

  // Tangential
  QLabel *m_label_t1;
  QDoubleSpinBox *m_t1;
  QLabel *m_label_t2;
  QDoubleSpinBox *m_t2;

  // Process
  QPushButton *m_process;
  QPushButton *m_exportResult;

  // Display
  QGraphicsView *m_view_old;
  QGraphicsView *m_view_new;
  QGraphicsScene *m_scene_old;
  QGraphicsScene *m_scene_new;

  // in/out images
  openMVG::image::Image<openMVG::image::RGBColor> m_inputImage;
  openMVG::image::Image<openMVG::image::RGBColor> m_outputImage;

  // Menu items
  QMenu *m_fileMenu;
  QAction *m_quit_act;
  QMenu *m_inputMenu;
  QAction *m_open_image_act;
  QAction *m_use_grid_act;
  QMenu *m_intrinsicMenu;
  QAction *m_no_dist_act;
  QAction *m_radial_k1_act;
  QAction *m_radial_k3_act;
  QAction *m_brown_act;
  QAction *m_fisheye_act;
  QMenu *m_processMenu;
  QAction *m_process_act;
  QAction *m_export_act;

  Q_OBJECT
};
} // namespace image_undistort_gui

#endif
