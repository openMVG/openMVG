// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SAMPLES_FEATURES_PAIR_DEMO_MAIN_WINDOW_HPP
#define OPENMVG_SAMPLES_FEATURES_PAIR_DEMO_MAIN_WINDOW_HPP

#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp"

#include "ImageView.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGraphicsEllipseItem>
#include <QGraphicsItem>
#include <QGraphicsLineItem>
#include <QGraphicsScene>
//#include <QGraphicsView>
#include <QLabel>
#include <QMainWindow>
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

private:
  /**
  * @brief Build all interface widgets 
  */
  void BuildInterface();

  /**
  * @brief Build all menus items 
  */
  void BuildMenus();

  /**
  * @brief Make connections between interface elements
  */
  void MakeConnections();

  /** 
  * @brief Get list of all feature type (and mode)
  */
  void PopulateFeatureType();

  /**
  * @brief Get list of all matching type
  */
  void PopulateMatchingType();

  /**
  * @brief remove previous computation on the view (only view data, not internal ones)
  */
  void ClearFeaturesAndMatchItems();

  /**
  * @brief Redraw match lines at correct position (ie:after image move)
  */
  void MoveMatcheLines();

  /**
  * @brief Get feature preset from interface
  * @return descriptor preset given the interface choice 
  */
  openMVG::features::EDESCRIBER_PRESET GetFeaturePreset( void );

  /**
  * @brief Get describer instance based on interface choice 
  * @param preset Current preset 
  * @return An image describer 
  */
  std::unique_ptr<openMVG::features::Image_describer> GetFeatureDescriber( const openMVG::features::EDESCRIBER_PRESET preset );

  /**
  * @brief Get Matcher type 
  * @return current matcher type based on interface choices 
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

  ImageView *m_view;
};

} // namespace features_pair_demo

#endif