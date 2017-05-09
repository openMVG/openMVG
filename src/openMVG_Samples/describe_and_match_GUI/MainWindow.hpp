// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SAMPLES_FEATURES_PAIR_DEMO_MAIN_WINDOW_HPP
#define OPENMVG_SAMPLES_FEATURES_PAIR_DEMO_MAIN_WINDOW_HPP

#include "openMVG/features/feature.hpp"
#include "openMVG/features/image_describer.hpp"
#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"

#include "ImageView.hpp"

#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGraphicsEllipseItem>
#include <QGraphicsItem>
#include <QGraphicsLineItem>
#include <QGraphicsScene>
#include <QLabel>
#include <QMainWindow>
#include <QMenu>
#include <QPixmap>
#include <QPushButton>

#include <memory>

namespace features_pair_demo
{

/**
  * @brief Main class
  */
class MainWindow : public QMainWindow
{
  Q_OBJECT
public:
  /**
  * @brief Constructor
  */
  MainWindow();

public slots:

  /**
  * @brief Action to be done when user wants to load the first image
  */
  void onOpenImage1( void );

  /**
  * @brief Action to be done when user wants to load the second image
  */
  void onOpenImage2( void );

  /**
  * @brief Action to be done when user wants to compute features and matching bewteen image pairs
  */
  void onComputeMatching( void );

  /**
  * @brief Action to be executed when user move image on the interface
  */
  void onMoveSomething( void );

  /**
  * @brief Action to be executed when user wants to export the final computation to an image
  */
  void onExportImage( void );

  /**
  * @brief Action to be executed when user request closing of image 1
  */
  void onCloseImage1( void );

  /**
  * @brief Action to be executed when user request closing of image 2
  */
  void onCloseImage2( void );

  /**
  * @brief Action to be executed when user request closing both images
  */
  void onCloseAll( void );

  /**
  * @brief Action to be executed when user request to quit the application
  */
  void onQuit( void );

  /**
  * @brief Action to be executed when user request clearing the computation
  */
  void onClearComputation( void );

private:
  /**
  * @brief Set elements enabled/disabled based on currently loaded elements
  */
  void UpdateActivation( void );

  /**
  * @brief Build all interface widgets
  */
  void BuildInterface();

  /**
  * @brief Build all menu items
  */
  void BuildMenus();

  /**
  * @brief Make connections between interface elements
  */
  void MakeConnections();

  /**
  * @brief Fill comboboxes with all feature types and feature presets
  */
  void PopulateFeatureType();

  /**
  * @brief Fill comboboxes with all matching types
  */
  void PopulateMatchingType();

  /**
  * @brief remove previous computation (internal data)
  */
  void ClearFeaturesandMatch();

  /**
  * @brief remove previous computation on the view (only view data, not internal ones)
  */
  void ClearFeaturesAndMatchItems();

  /**
  * @brief Redraw match lines at the correct position (i.e. :after image move)
  */
  void MoveMatchLines();

  /**
  * @brief Get feature preset from interface
  * @return descriptor preset given the interface choice
  */
  openMVG::features::EDESCRIBER_PRESET GetFeaturePreset( void );

  /**
  * @brief Get describer instance based on the choices made in the interface
  * @param preset Current preset
  * @return An image describer
  */
  std::unique_ptr<openMVG::features::Image_describer> GetFeatureDescriber( const openMVG::features::EDESCRIBER_PRESET preset );

  /**
  * @brief Get Matcher type
  * @return current matcher type based on the choices made in the interface
  */
  openMVG::matching::EMatcherType GetMatcherType( void );

  QLabel *m_labelFeatureName;
  QComboBox *m_comboFeatureName;
  QLabel *m_labelFeatureMode;
  QComboBox *m_comboFeatureMode;
  QCheckBox *m_checkUpright;

  QLabel *m_matchingType;
  QComboBox *m_comboMatchingName;
  QLabel *m_labelDistanceThreshold;
  QDoubleSpinBox *m_spinDistThreshold;

  QPushButton *m_loadImage1;
  QPushButton *m_loadImage2;
  QPushButton *m_computeMatching;

  QPushButton *m_exportResult;

  QGraphicsScene *m_scn;
  ImageView *m_view;

  std::string m_image1Path;
  QPixmap *m_pixmap1;
  QGraphicsItem *m_pixmap1Item;

  std::string m_image2Path;
  QPixmap *m_pixmap2;
  QGraphicsItem *m_pixmap2Item;

  std::vector<openMVG::features::PointFeature> m_feats1;
  std::vector<QGraphicsEllipseItem *> m_ellipsesFeat1;
  std::vector<openMVG::features::PointFeature> m_feats2;
  std::vector<QGraphicsEllipseItem *> m_ellipsesFeat2;
  openMVG::matching::IndMatches m_matches;
  std::vector<QGraphicsLineItem *> m_lineMatch;

  // menus items
  QMenu *m_fileMenu;
  QMenu *m_processingMenu;

  QAction *m_fileOpenImage1;
  QAction *m_fileOpenImage2;
  QAction *m_fileCloseImage1;
  QAction *m_fileCloseImage2;
  QAction *m_fileCloseAll;
  QAction *m_fileExportImage;
  QAction *m_fileQuit;

  QAction *m_processingComputeMatch;
  QAction *m_processingClearMatch;
};

} // namespace features_pair_demo

#endif
