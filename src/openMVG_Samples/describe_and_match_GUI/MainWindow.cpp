// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "MainWindow.hpp"

#include "nonFree/sift/SIFT_describer.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/matching/regions_matcher.hpp"

#include <QApplication>
#include <QDir>
#include <QFileDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QMenuBar>
#include <QProgressDialog>
#include <QString>

#include <memory>
#include <random>

using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::matching;
using namespace openMVG::matching_image_collection;

namespace features_pair_demo
{

/**
  * @brief Constructor
  */
MainWindow::MainWindow()
{
  BuildInterface();
  BuildMenus();
  MakeConnections();

  UpdateActivation();

  setWindowTitle( "Features pair demo" );
  resize( 1024, 768 );
}

/**
  * @brief Action to be done when user wants to load the first image
  */
void MainWindow::onOpenImage1( void )
{
  QDir home    = QDir::home();
  QString file = QFileDialog::getOpenFileName( this, "Open Image 1", home.path(), tr( "Image Files (*.png *.jpg *.pbm *.pgm *.ppm *.tiff)" ) );
  if ( !file.isEmpty() && !file.isNull() )
  {
    QImage img( file );
    *m_pixmap1 = QPixmap::fromImage( img );
    m_scn->removeItem( m_pixmap1Item );
    m_pixmap1Item = m_scn->addPixmap( *m_pixmap1 );
    m_pixmap1Item->setFlag( QGraphicsItem::ItemIsMovable, true );

    // Move the first image relative to the second one
    QPointF pos = m_pixmap2Item->pos();
    pos.setX( pos.x() - m_pixmap1Item->boundingRect().width() );
    m_pixmap1Item->setPos( m_pixmap1Item->mapToScene( pos ) );

    m_view->fitInView( m_scn->itemsBoundingRect(), Qt::KeepAspectRatio );

    m_image1Path = file.toStdString();

    ClearFeaturesandMatch();
    ClearFeaturesAndMatchItems();

    UpdateActivation();
  }
}

/**
  * @brief Action to be done when user wants to load the second image
  */
void MainWindow::onOpenImage2( void )
{
  QDir home    = QDir::home();
  QString file = QFileDialog::getOpenFileName( this, "Open Image 2", home.path(), tr( "Image Files (*.png *.jpg *.pbm *.pgm *.ppm *.tiff)" ) );
  if ( !file.isEmpty() && !file.isNull() )
  {
    QImage img( file );
    *m_pixmap2 = QPixmap::fromImage( img );
    m_scn->removeItem( m_pixmap2Item );
    m_pixmap2Item = m_scn->addPixmap( *m_pixmap2 );
    m_pixmap2Item->setFlag( QGraphicsItem::ItemIsMovable, true );

    // Move the second image relative to the first one
    QPointF pos = m_pixmap1Item->pos();
    pos.setX( pos.x() + m_pixmap1Item->boundingRect().width() );
    m_pixmap2Item->setPos( m_pixmap2Item->mapToScene( pos ) );

    m_view->fitInView( m_scn->itemsBoundingRect(), Qt::KeepAspectRatio );

    m_image2Path = file.toStdString();

    ClearFeaturesandMatch();
    ClearFeaturesAndMatchItems();

    UpdateActivation();
  }
}

/**
  * @brief Action to be done when user wants to compute features and matching bewteen image pairs
  */
void MainWindow::onComputeMatching( void )
{
  if ( !( m_image1Path.size() > 0 && m_image2Path.size() > 0 ) )
  {
    // TODO : add a message box to tell that we want two image before processing
    return;
  }

  // 0 clear previously computed values
  ClearFeaturesAndMatchItems();
  ClearFeaturesandMatch();

  QProgressDialog progress( "Computation in progress ...", "Abort", 0, 4, this );
  progress.setWindowModality( Qt::WindowModal );
  progress.setMinimumDuration( 0 );

  // 1.0 -> Load images from names
  Image<unsigned char> image1, image2;
  ReadImage( m_image1Path.c_str(), &image1 );
  if ( progress.wasCanceled() )
    return;
  ReadImage( m_image2Path.c_str(), &image2 );

  if ( progress.wasCanceled() )
    return;
  progress.setValue( 1 );

  // 1.1 Get the describer
  EDESCRIBER_PRESET preset                         = GetFeaturePreset();
  std::unique_ptr<Image_describer> image_describer = GetFeatureDescriber( preset );

  // 1.2 Compute description of the two images
  Image<unsigned char> *mask = nullptr; // The mask is null by default

  std::unique_ptr<Regions> regions1, regions2;
  image_describer->Describe( image1, regions1, mask );
  if ( progress.wasCanceled() )
    return;
  image_describer->Describe( image2, regions2, mask );
  if ( progress.wasCanceled() )
    return;

  progress.setValue( 2 );

  // 2.1 Get matching params
  const EMatcherType matcherType = GetMatcherType();
  const float fDistRatio         = static_cast<float>( m_spinDistThreshold->value() );

  m_matches.clear();

  // 2.1 Matches features
  DistanceRatioMatch(
      fDistRatio, matcherType,
      *regions1.get(),
      *regions2.get(),
      m_matches );

  if ( progress.wasCanceled() )
  {
    m_matches.clear();
    return;
  }
  progress.setValue( 3 );

  // 3 draw result

  m_feats1 = regions1->GetRegionsPositions();
  m_feats2 = regions2->GetRegionsPositions();

  // Position anchor of the images (in global space and local space)
  const QPointF gpos1 = m_pixmap1Item->pos();
  const QPointF gpos2 = m_pixmap2Item->pos();

  const QPointF pos1 = m_pixmap1Item->mapFromScene( gpos1 );
  const QPointF pos2 = m_pixmap2Item->mapFromScene( gpos2 );

  // Default color for unmatched points
  QPen penFeat( QColor( 0, 0, 0 ) );

  // Generate random color per pair
  std::map<int, QPen>
   feat1Pen,
   feat2Pen,
   matchesPen;

  std::random_device rd;
  std::mt19937 rng( rd() );

  std::uniform_int_distribution<int>
    distrib_h( 0, 359 ),
    distrib_s( 200, 255 ),
    distrib_v( 200, 255 );

  for ( size_t i = 0; i < m_matches.size(); ++i )
  {
    const int id1 = m_matches[ i ].i_;
    const int id2 = m_matches[ i ].j_;

    const QColor col = QColor::fromHsv( distrib_h( rng ), distrib_s( rng ), distrib_v( rng ) );

    QPen curPen( col );
    curPen.setWidth( 2 );
    feat1Pen[ id1 ] = curPen;
    feat2Pen[ id2 ] = curPen;
    matchesPen[ i ] = curPen;
  }

  for ( size_t i = 0; i < m_feats1.size(); ++i )
  {
    QPen &curPen              = feat1Pen.count( i ) ? feat1Pen[ i ] : penFeat;
    const PointFeature &f     = m_feats1[ i ];
    QGraphicsEllipseItem *ell = m_scn->addEllipse( f.x() + pos1.x(), f.y() + pos1.y(), 10.0, 10.0, curPen );
    ell->setParentItem( m_pixmap1Item );

    m_ellipsesFeat1.emplace_back( ell );
  }

  for ( size_t i = 0; i < m_feats2.size(); ++i )
  {
    QPen &curPen              = feat2Pen.count( i ) ? feat2Pen[ i ] : penFeat;
    const PointFeature &f     = m_feats2[ i ];
    QGraphicsEllipseItem *ell = m_scn->addEllipse( f.x() + pos2.x(), f.y() + pos2.y(), 10.0, 10.0, curPen );
    ell->setParentItem( m_pixmap2Item );

    m_ellipsesFeat2.emplace_back( ell );
  }

  // Correspondances
  for ( size_t i = 0; i < m_matches.size(); ++i )
  {
    QPen curPen           = matchesPen[ i ];
    const PointFeature &L = m_feats1[ m_matches[ i ].i_ ];
    const PointFeature &R = m_feats2[ m_matches[ i ].j_ ];

    m_lineMatch.emplace_back( m_scn->addLine( L.x() + gpos1.x(), L.y() + gpos1.y(), R.x() + gpos2.x(), R.y() + gpos2.y(), curPen ) );
  }

  progress.setValue( 4 );

  UpdateActivation();
}

/**
  * @brief Action to be executed when user move image on the interface
  */
void MainWindow::onMoveSomething( void )
{
  QGraphicsItem *item = m_scn->mouseGrabberItem();
  if ( item == m_pixmap1Item ||
       item == m_pixmap2Item )
  {
    MoveMatchLines();
  }
}

/**
  * @brief Action to be executed when user wants to export the final computation to an image
  */
void MainWindow::onExportImage( void )
{
  QDir home    = QDir::home();
  QString file = QFileDialog::getSaveFileName( this, "Save File",
                                               home.path(),
                                               tr( "Images (*.png *.xpm *.jpg)" ) );
  if ( !file.isEmpty() && !file.isNull() )
  {
    QPixmap pixMap = m_view->grab();
    pixMap.save( file );
  }
}

/**
  * @brief Action to be executed when user request closing of image 1
  */
void MainWindow::onCloseImage1( void )
{
  m_feats1.clear();
  m_matches.clear();
  m_image1Path.clear();

  for ( auto item : m_ellipsesFeat1 )
  {
    m_scn->removeItem( item );
  }
  for ( auto item : m_lineMatch )
  {
    m_scn->removeItem( item );
  }
  m_ellipsesFeat1.clear();
  m_lineMatch.clear();

  m_scn->removeItem( m_pixmap1Item );

  UpdateActivation();
}

/**
  * @brief Action to be executed when user request closing of image 2
  */
void MainWindow::onCloseImage2( void )
{
  m_feats2.clear();
  m_matches.clear();
  m_image2Path.clear();

  for ( auto item : m_ellipsesFeat2 )
  {
    m_scn->removeItem( item );
  }
  for ( auto item : m_lineMatch )
  {
    m_scn->removeItem( item );
  }
  m_ellipsesFeat2.clear();
  m_lineMatch.clear();

  m_scn->removeItem( m_pixmap2Item );

  UpdateActivation();
}

/**
  * @brief Action to be executed when user request closing both images
  */
void MainWindow::onCloseAll( void )
{
  onCloseImage1();
  onCloseImage2();

  UpdateActivation();
}

/**
  * @brief Action to be executed when user request to quit the application
  */
void MainWindow::onQuit( void )
{
  QApplication::quit();
}

/**
  * @brief Action to be executed when user request clearing the computation
  */
void MainWindow::onClearComputation( void )
{
  ClearFeaturesandMatch();
  ClearFeaturesAndMatchItems();

  UpdateActivation();
}

/**
  * @brief Set elements enabled/disabled based on currently loaded elements
  */
void MainWindow::UpdateActivation( void )
{
  // has loaded image 1 so we can close it
  if ( m_image1Path.size() > 0 )
  {
    m_fileCloseImage1->setEnabled( true );
    m_fileCloseAll->setEnabled( true );
  }
  else
  {
    m_fileCloseImage1->setEnabled( false );
  }

  // has loaded image 2 so we can close it
  if ( m_image2Path.size() > 0 )
  {
    m_fileCloseImage2->setEnabled( true );
    m_fileCloseAll->setEnabled( true );
  }
  else
  {
    m_fileCloseImage2->setEnabled( false );
  }

  // No image opened , close all should be removed
  if ( !( m_image1Path.size() > 0 || m_image2Path.size() > 0 ) )
  {
    m_fileCloseAll->setEnabled( false );
  }

  // Both loaded, could do computation and exportation
  if ( m_image1Path.size() > 0 && m_image2Path.size() > 0 )
  {
    m_computeMatching->setEnabled( true );
    m_processingClearMatch->setEnabled( true );
    m_processingComputeMatch->setEnabled( true );
    m_fileExportImage->setEnabled( true );
  }
  else
  {

    m_computeMatching->setEnabled( false );
    m_processingClearMatch->setEnabled( false );
    m_processingComputeMatch->setEnabled( false );
    m_fileExportImage->setEnabled( false );
  }
}

/**
  * @brief Build all interface widgets
  */
void MainWindow::BuildInterface()
{
  // Features
  QGridLayout *featureLayout = new QGridLayout;
  QGroupBox *groupFeatures   = new QGroupBox( "Feature" );
  m_labelFeatureName         = new QLabel( "Feature" );
  m_labelFeatureMode         = new QLabel( "Mode" );
  m_checkUpright             = new QCheckBox( "Upright" );
  m_checkUpright->setTristate( false );
  m_checkUpright->setCheckState( Qt::Unchecked );

  m_comboFeatureName = new QComboBox;
  m_comboFeatureMode = new QComboBox;

  featureLayout->addWidget( m_labelFeatureName, 0, 0 );
  featureLayout->addWidget( m_comboFeatureName, 0, 1 );
  featureLayout->addWidget( m_labelFeatureMode, 0, 2 );
  featureLayout->addWidget( m_comboFeatureMode, 0, 3 );
  featureLayout->addWidget( m_checkUpright, 0, 4 );

  groupFeatures->setLayout( featureLayout );

  // Matching
  QGroupBox *groupMatching    = new QGroupBox( "Matching" );
  QGridLayout *matchingLayout = new QGridLayout;

  m_matchingType           = new QLabel( "Matching" );
  m_comboMatchingName      = new QComboBox;
  m_labelDistanceThreshold = new QLabel( "Distance" );
  m_spinDistThreshold      = new QDoubleSpinBox();
  m_spinDistThreshold->setDecimals( 3 );
  m_spinDistThreshold->setValue( 0.8 );

  matchingLayout->addWidget( m_matchingType, 0, 0 );
  matchingLayout->addWidget( m_comboMatchingName, 0, 1 );
  matchingLayout->addWidget( m_labelDistanceThreshold, 0, 2 );
  matchingLayout->addWidget( m_spinDistThreshold, 0, 3 );

  groupMatching->setLayout( matchingLayout );

  // Processing
  QGroupBox *groupProcessing    = new QGroupBox( "Processing" );
  QGridLayout *processingLayout = new QGridLayout;

  m_loadImage1      = new QPushButton( "Load Image 1" );
  m_loadImage2      = new QPushButton( "Load Image 2" );
  m_computeMatching = new QPushButton( "Match" );

  processingLayout->addWidget( m_loadImage1, 0, 0 );
  processingLayout->addWidget( m_computeMatching, 0, 1 );
  processingLayout->addWidget( m_loadImage2, 0, 2 );

  groupProcessing->setLayout( processingLayout );

  // Result
  QGroupBox *groupResult    = new QGroupBox( "Result" );
  QGridLayout *resultLayout = new QGridLayout;

  m_scn  = new QGraphicsScene;
  m_view = new ImageView( m_scn );
  m_view->setDragMode( QGraphicsView::ScrollHandDrag );
  m_pixmap1     = new QPixmap;
  m_pixmap2     = new QPixmap;
  m_pixmap1Item = m_scn->addPixmap( *m_pixmap1 );
  m_pixmap2Item = m_scn->addPixmap( *m_pixmap2 );

  resultLayout->addWidget( m_view );

  groupResult->setLayout( resultLayout );

  m_exportResult = new QPushButton( "Export image" );

  // Add everything to the window
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget( groupFeatures );
  mainLayout->addWidget( groupMatching );
  mainLayout->addWidget( groupProcessing );
  mainLayout->addWidget( groupResult );
  mainLayout->addWidget( m_exportResult );

  QWidget *dummyWidget = new QWidget;
  dummyWidget->setLayout( mainLayout );
  setCentralWidget( dummyWidget );

  PopulateFeatureType();
  PopulateMatchingType();
}

/**
  * @brief Build all menus items
  */
void MainWindow::BuildMenus()
{
  // File
  m_fileMenu = new QMenu( "File" );

  m_fileOpenImage1 = m_fileMenu->addAction( "Open image 1" );
  m_fileOpenImage2 = m_fileMenu->addAction( "Open image 2" );
  m_fileMenu->addSeparator();
  m_fileCloseImage1 = m_fileMenu->addAction( "Close image 1" );
  m_fileCloseImage2 = m_fileMenu->addAction( "Close image 2" );
  m_fileCloseAll    = m_fileMenu->addAction( "Close all" );
  m_fileCloseAll->setShortcut( QKeySequence::Close );
  m_fileMenu->addSeparator();
  m_fileExportImage = m_fileMenu->addAction( "Save result image" );
  m_fileExportImage->setShortcut( QKeySequence::Save );
  m_fileMenu->addSeparator();
  m_fileQuit = m_fileMenu->addAction( "Quit" );
  m_fileQuit->setShortcut( QKeySequence::Quit );

  // Processing
  m_processingMenu = new QMenu( "Processing" );

  m_processingComputeMatch = m_processingMenu->addAction( "Compute match" );
  m_processingClearMatch   = m_processingMenu->addAction( "Clear computation" );

  QMenuBar *menuBar = this->menuBar();

  menuBar->addMenu( m_fileMenu );
  menuBar->addMenu( m_processingMenu );
}

/**
  * @brief Make connections between interface elements
  */
void MainWindow::MakeConnections()
{
  // Buttons
  connect( m_loadImage1, SIGNAL( clicked() ), this, SLOT( onOpenImage1() ) );
  connect( m_loadImage2, SIGNAL( clicked() ), this, SLOT( onOpenImage2() ) );
  connect( m_computeMatching, SIGNAL( clicked() ), this, SLOT( onComputeMatching() ) );
  connect( m_exportResult, SIGNAL( clicked() ), this, SLOT( onExportImage() ) );

  // Moving an image
  connect( m_scn, SIGNAL( changed( const QList<QRectF> & ) ), this, SLOT( onMoveSomething() ) );

  // menus
  connect( m_fileQuit, SIGNAL( triggered() ), this, SLOT( onQuit() ) );
  connect( m_fileOpenImage1, SIGNAL( triggered() ), this, SLOT( onOpenImage1() ) );
  connect( m_fileOpenImage2, SIGNAL( triggered() ), this, SLOT( onOpenImage2() ) );
  connect( m_fileCloseImage1, SIGNAL( triggered() ), this, SLOT( onCloseImage1() ) );
  connect( m_fileCloseImage2, SIGNAL( triggered() ), this, SLOT( onCloseImage2() ) );
  connect( m_fileCloseAll, SIGNAL( triggered() ), this, SLOT( onCloseAll() ) );
  connect( m_fileExportImage, SIGNAL( triggered() ), this, SLOT( onExportImage() ) );
  connect( m_processingComputeMatch, SIGNAL( triggered() ), this, SLOT( onComputeMatching() ) );
  connect( m_processingClearMatch, SIGNAL( triggered() ), this, SLOT( onClearComputation() ) );
}

/**
  * @brief Get list of all feature type (and mode)
  */
void MainWindow::PopulateFeatureType()
{
  m_comboFeatureName->addItem( "SIFT" );
  m_comboFeatureName->addItem( "SIFT_ANATOMY" );
  m_comboFeatureName->addItem( "AKAZE_FLOAT" );
  m_comboFeatureName->addItem( "AKAZE_MLDB" );

  m_comboFeatureMode->addItem( "NORMAL" );
  m_comboFeatureMode->addItem( "HIGH" );
  m_comboFeatureMode->addItem( "ULTRA" );
}

/**
  * @brief Get list of all matching type
  */
void MainWindow::PopulateMatchingType()
{
  m_comboMatchingName->addItem( "BRUTEFORCEL2" );
  m_comboMatchingName->addItem( "ANNL2" );
  m_comboMatchingName->addItem( "CASCADEHASHINGL2" );
  m_comboMatchingName->insertSeparator( 4 );
  m_comboMatchingName->addItem( "BRUTEFORCEHAMMING" );

  m_comboMatchingName->setCurrentIndex( 1 );
}

/**
  * @brief remove previous computation (internal data)
  */
void MainWindow::ClearFeaturesandMatch()
{
  m_feats1.clear();
  m_feats2.clear();
  m_matches.clear();
}

/**
  * @brief remove previous computation on the view (only view data, not internal ones)
  */
void MainWindow::ClearFeaturesAndMatchItems()
{
  for (auto item : m_ellipsesFeat1 )
  {
    m_scn->removeItem( item );
  }
  for (auto item : m_ellipsesFeat2 )
  {
    m_scn->removeItem( item );
  }
  for (auto item : m_lineMatch )
  {
    m_scn->removeItem( item );
  }
  m_ellipsesFeat1.clear();
  m_ellipsesFeat2.clear();
  m_lineMatch.clear();
}

/**
  * @brief Redraw match lines at correct position (ie:after image move)
  */
void MainWindow::MoveMatchLines()
{
  // Position anchor of the images (in global space and local space)
  QPointF gpos1 = m_pixmap1Item->pos();
  QPointF gpos2 = m_pixmap2Item->pos();

  for ( size_t i = 0; i < m_matches.size(); ++i )
  {
    const PointFeature &L = m_feats1[ m_matches[ i ].i_ ];
    const PointFeature &R = m_feats2[ m_matches[ i ].j_ ];

    QGraphicsLineItem *cur_line = m_lineMatch[ i ];

    cur_line->setLine( L.x() + gpos1.x(), L.y() + gpos1.y(), R.x() + gpos2.x(), R.y() + gpos2.y() );
  }
}

/**
  * @brief Get feature preset from interface
  * @return descriptor preset given the interface choice
  */
EDESCRIBER_PRESET MainWindow::GetFeaturePreset()
{
  // get selected item
  const std::string sPreset = m_comboFeatureMode->currentText().toStdString();
  EDESCRIBER_PRESET preset;
  if ( sPreset == "NORMAL" )
  {
    preset = NORMAL_PRESET;
  }
  else if ( sPreset == "HIGH" )
  {
    preset = HIGH_PRESET;
  }
  else if ( sPreset == "ULTRA" )
  {
    preset = ULTRA_PRESET;
  }
  else
  {
    std::cerr << "Unknown preset value" << std::endl;
  }

  return preset;
}

/**
  * @brief Get describer instance based on interface choice
  * @param preset Current preset
  * @return An image describer
  */
std::unique_ptr<Image_describer> MainWindow::GetFeatureDescriber( const EDESCRIBER_PRESET preset )
{
  const std::string sImage_Describer_Method = m_comboFeatureName->currentText().toStdString();
  const bool bUpRight                       = m_checkUpright->checkState() == Qt::Checked;

  std::unique_ptr<Image_describer> image_describer;

  if ( sImage_Describer_Method == "SIFT" )
  {
    image_describer.reset( new SIFT_Image_describer( SIFT_Image_describer::Params(), !bUpRight ) );
  }
  else if ( sImage_Describer_Method == "SIFT_ANATOMY" )
  {
    image_describer.reset( new SIFT_Anatomy_Image_describer( SIFT_Anatomy_Image_describer::Params() ) );
  }
  else if ( sImage_Describer_Method == "AKAZE_FLOAT" )
  {
    image_describer = AKAZE_Image_describer::create( AKAZE_Image_describer::Params( AKAZE::Params(), AKAZE_MSURF ), !bUpRight );
  }
  else if ( sImage_Describer_Method == "AKAZE_MLDB" )
  {
    image_describer = AKAZE_Image_describer::create( AKAZE_Image_describer::Params( AKAZE::Params(), AKAZE_MLDB ), !bUpRight );
  }
  else
  {
    std::cerr << "Unknown feature value" << std::endl;
  }

  image_describer->Set_configuration_preset( preset );

  return image_describer;
}

/**
  * @brief Get Matcher type
  * @return current matcher type based on interface choices
  */
EMatcherType MainWindow::GetMatcherType( void )
{
  const std::string sNearestMatchingMethod = m_comboMatchingName->currentText().toStdString();
  if ( sNearestMatchingMethod == "BRUTEFORCEL2" )
  {
    return BRUTE_FORCE_L2;
  }
  else if ( sNearestMatchingMethod == "BRUTEFORCEHAMMING" )
  {
    return BRUTE_FORCE_HAMMING;
  }
  else if ( sNearestMatchingMethod == "ANNL2" )
  {
    return ANN_L2;
  }
  else if ( sNearestMatchingMethod == "CASCADEHASHINGL2" )
  {
    return CASCADE_HASHING_L2;
  }

  return BRUTE_FORCE_L2;
}

} // namespace features_pair_demo
