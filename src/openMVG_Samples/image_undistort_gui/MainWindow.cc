// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "MainWindow.hh"

#include <iostream>

#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"
#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/image/image_io.hpp"

#include "QImageInterface.hh"

#include <QCoreApplication>
#include <QDir>
#include <QFileDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QMenuBar>

namespace image_undistort_gui
{

const int MainWindow::STATE_INITIAL                 = 0;
const int MainWindow::STATE_HAS_INPUT               = 1;
const int MainWindow::STATE_HAS_COMPUTED_DISTORTION = 2;

/**
  * @brief Constructor
  * @param parent Parent widget of this window 
  */
MainWindow::MainWindow( QWidget *parent )
    : QMainWindow( parent ),
      m_current_state( STATE_INITIAL )
{
  BuildInterface();
  BuildMenus();
  MakeConnections();

  resize( 1024, 768 );
  setWindowTitle( "ImageUndistort" );

  UpdateInterface();
}

/**
  * @brief Action to be executed when user want to open an image 
  */
void MainWindow::onOpenImage( void )
{
  QString fileName = QFileDialog::getOpenFileName( this, tr( "Open File" ),
                                                   QDir::homePath(),
                                                   tr( "Images (*.png *.jpg *.tiff)" ) );
  if ( !( fileName.isNull() || fileName.isEmpty() ) )
  {
    const std::string imagePath = fileName.toStdString();
    QImage tmp( fileName );
    m_inputImage = QImageToOpenMVGImage( tmp );

    QPixmap pix = QPixmap::fromImage( tmp );
    m_scene_old->clear();
    m_scene_old->addPixmap( pix );
    m_view_old->fitInView( m_scene_old->itemsBoundingRect(), Qt::KeepAspectRatio );

    m_scene_new->clear();

    const double w = m_inputImage.Width();
    const double h = m_inputImage.Height();

    m_ppx->setValue( w / 2.0 );
    m_ppy->setValue( h / 2.0 );
    m_fov_spin->setValue( 1.1 * std::max( w, h ) );

    m_current_state = STATE_HAS_INPUT;
    UpdateInterface();
  }
}

/**
  * @brief Action to be executed when user want to undistort an image 
  */
void MainWindow::onProcessImage( void )
{
  if ( m_inputImage.Width() > 0 && m_inputImage.Height() > 0 )
  {
    openMVG::cameras::IntrinsicBase *intrin = nullptr;

    const int w      = m_inputImage.Width();
    const int h      = m_inputImage.Height();
    const double fov = m_fov_spin->value();
    const double ppx = m_ppx->value();
    const double ppy = m_ppy->value();

    const double k1 = m_k1->value();
    const double k2 = m_k2->value();
    const double k3 = m_k3->value();
    const double k4 = m_k4->value();
    const double t1 = m_t1->value();
    const double t2 = m_t2->value();

    switch ( m_intrinsic_type->currentIndex() )
    {
    case 0:
    {
      intrin = new openMVG::cameras::Pinhole_Intrinsic( w, h, fov, ppx, ppy );
      break;
    }
    case 1:
    {
      intrin = new openMVG::cameras::Pinhole_Intrinsic_Radial_K1( w, h, fov, ppx, ppy, k1 );
      break;
    }
    case 2:
    {
      intrin = new openMVG::cameras::Pinhole_Intrinsic_Radial_K3( w, h, fov, ppx, ppy, k1, k2, k3 );
      break;
    }
    case 3:
    {
      intrin = new openMVG::cameras::Pinhole_Intrinsic_Brown_T2( w, h, fov, ppx, ppy, k1, k2, k3, t1, t2 );
      break;
    }
    case 4:
    {
      intrin = new openMVG::cameras::Pinhole_Intrinsic_Fisheye( w, h, fov, ppx, ppy, k1, k2, k3, k4 );
      break;
    }
    }

    openMVG::cameras::UndistortImageResized( m_inputImage, intrin, m_outputImage );

    QImage tmp = openMVGImageToQImage( m_outputImage );

    QPixmap pix = QPixmap::fromImage( tmp );
    m_scene_new->clear();
    m_scene_new->addPixmap( pix );
    m_view_new->fitInView( m_scene_new->itemsBoundingRect(), Qt::KeepAspectRatio );

    m_current_state = STATE_HAS_COMPUTED_DISTORTION;
    UpdateInterface();
  }
}

/**
  * @brief Action to be executed when user want to export an undistorted image 
  */
void MainWindow::onExportImage( void )
{
  QString fileName = QFileDialog::getSaveFileName( this, "Save image",
                                                   QDir::homePath(),
                                                   tr( "Images (*.png *.jpg *.tiff)" ) );
  if ( !( fileName.isEmpty() || fileName.isNull() ) )
  {
    const std::string name = fileName.toStdString();
    openMVG::image::WriteImage( name.c_str(), m_outputImage );
  }
}

/**
  * @brief Action to be executed when user select an intrinsic type 
  */
void MainWindow::onSelectIntrinsicType( void )
{
  const int idx = m_intrinsic_type->currentIndex();
  switch ( idx )
  {
  case 0:
  { // No dist
    m_k1->setEnabled( false );
    m_k2->setEnabled( false );
    m_k3->setEnabled( false );
    m_k4->setEnabled( false );
    m_t1->setEnabled( false );
    m_t2->setEnabled( false );

    m_label_k1->setVisible( false );
    m_k1->setVisible( false );
    m_label_k2->setVisible( false );
    m_k2->setVisible( false );
    m_label_k3->setVisible( false );
    m_k3->setVisible( false );
    m_label_k4->setVisible( false );
    m_k4->setVisible( false );
    m_label_t1->setVisible( false );
    m_t1->setVisible( false );
    m_label_t2->setVisible( false );
    m_t2->setVisible( false );

    break;
  }
  case 1:
  { // R1
    m_k1->setEnabled( true );
    m_k2->setEnabled( false );
    m_k3->setEnabled( false );
    m_k4->setEnabled( false );
    m_t1->setEnabled( false );
    m_t2->setEnabled( false );

    m_label_k1->setVisible( true );
    m_k1->setVisible( true );
    m_label_k2->setVisible( false );
    m_k2->setVisible( false );
    m_label_k3->setVisible( false );
    m_k3->setVisible( false );
    m_label_k4->setVisible( false );
    m_k4->setVisible( false );
    m_label_t1->setVisible( false );
    m_t1->setVisible( false );
    m_label_t2->setVisible( false );
    m_t2->setVisible( false );

    break;
  }
  case 2:
  { // R3
    m_k1->setEnabled( true );
    m_k2->setEnabled( true );
    m_k3->setEnabled( true );
    m_k4->setEnabled( false );
    m_t1->setEnabled( false );
    m_t2->setEnabled( false );

    m_label_k1->setVisible( true );
    m_k1->setVisible( true );
    m_label_k2->setVisible( true );
    m_k2->setVisible( true );
    m_label_k3->setVisible( true );
    m_k3->setVisible( true );
    m_label_k4->setVisible( false );
    m_k4->setVisible( false );
    m_label_t1->setVisible( false );
    m_t1->setVisible( false );
    m_label_t2->setVisible( false );
    m_t2->setVisible( false );

    break;
  }
  case 3:
  {
    // Brown R3 - T2
    m_k1->setEnabled( true );
    m_k2->setEnabled( true );
    m_k3->setEnabled( true );
    m_k4->setEnabled( false );
    m_t1->setEnabled( true );
    m_t2->setEnabled( true );

    m_label_k1->setVisible( true );
    m_k1->setVisible( true );
    m_label_k2->setVisible( true );
    m_k2->setVisible( true );
    m_label_k3->setVisible( true );
    m_k3->setVisible( true );
    m_label_k4->setVisible( false );
    m_k4->setVisible( false );
    m_label_t1->setVisible( true );
    m_t1->setVisible( true );
    m_label_t2->setVisible( true );
    m_t2->setVisible( true );

    break;
  }
  case 4:
  {
    // Fisheye R4
    m_k1->setEnabled( true );
    m_k2->setEnabled( true );
    m_k3->setEnabled( true );
    m_k4->setEnabled( true );
    m_t1->setEnabled( false );
    m_t2->setEnabled( false );

    m_label_k1->setVisible( true );
    m_k1->setVisible( true );
    m_label_k2->setVisible( true );
    m_k2->setVisible( true );
    m_label_k3->setVisible( true );
    m_k3->setVisible( true );
    m_label_k4->setVisible( true );
    m_k4->setVisible( true );
    m_label_t1->setVisible( false );
    m_t1->setVisible( false );
    m_label_t2->setVisible( false );
    m_t2->setVisible( false );

    break;
  }
  }
}

/**
  * @brief Action to be executed when user check/uncheck the grid checkbox
  */
void MainWindow::onChangeCheckBox( void )
{
  if ( m_check_use_grid->isChecked() )
  {
    BuildGrid();
    m_scene_new->clear();

    m_current_state = STATE_HAS_INPUT;
  }
  else
  {
    m_inputImage = openMVG::image::Image<openMVG::image::RGBColor>( 10, 10 );
    m_scene_old->clear();
    m_scene_new->clear();

    m_current_state = STATE_INITIAL;
  }

  UpdateInterface();
}

/**
  * @brief action to be executed when user select no distorsion menu item
  */
void MainWindow::onMenuSelectNoDist( void )
{
  m_intrinsic_type->setCurrentIndex( 0 );
  UpdateInterface();
}

/**
  * @brief action to be executed when user select radial 1 menu item
  */
void MainWindow::onMenuSelectRadial1( void )
{
  m_intrinsic_type->setCurrentIndex( 1 );
  UpdateInterface();
}

/**
  * @brief action to be executed when user select radial 3 menu item
  */
void MainWindow::onMenuSelectRadial3( void )
{
  m_intrinsic_type->setCurrentIndex( 2 );
  UpdateInterface();
}

/**
  * @brief action to be executed when user select brown menu item
  */
void MainWindow::onMenuSelectBrown( void )
{
  m_intrinsic_type->setCurrentIndex( 3 );
  UpdateInterface();
}

/**
  * @brief action to be executed when user select fisheye menu item
  */
void MainWindow::onMenuSelectFisheye( void )
{
  m_intrinsic_type->setCurrentIndex( 4 );
  UpdateInterface();
}

/**
  * @brief action to be executed when user select quit menu item
  */
void MainWindow::onQuit( void )
{
  QCoreApplication::quit();
}

/**
  * @brief Fill intrinsic values in the combobox 
  */
void MainWindow::FillIntrinsicTypes( void )
{
  m_intrinsic_type->addItem( "No distorsion" );
  m_intrinsic_type->addItem( "Radial 1" );
  m_intrinsic_type->addItem( "Radial 3" );
  m_intrinsic_type->addItem( "Brown : Radial 3 - tangential 2" );
  m_intrinsic_type->addItem( "Fisheye" );
}

/**
  * @brief Build virtual image representing a regular grid 
  */
void MainWindow::BuildGrid( void )
{
  int w        = 1280;
  int h        = 720;
  m_inputImage = openMVG::image::Image<openMVG::image::RGBColor>( w, h );

  int nb_x = 16;
  int nb_y = 9;

  const double delta_x = w / nb_x;
  const double delta_y = h / nb_y;

  // Vertical lines
  for ( int i = 1; i < nb_x; ++i )
  {
    const int x = i * delta_x;
    for ( int y = 0; y < h; ++y )
    {
      m_inputImage( y, x - 1 ) = openMVG::image::RGBColor( 255, 255, 255 );
      m_inputImage( y, x )     = openMVG::image::RGBColor( 255, 255, 255 );
      m_inputImage( y, x + 1 ) = openMVG::image::RGBColor( 255, 255, 255 );
    }
  }

  // Horizontal lines
  for ( int i = 1; i < nb_y; ++i )
  {
    const int y = i * delta_y;
    for ( int x = 0; x < w; ++x )
    {
      m_inputImage( y - 1, x ) = openMVG::image::RGBColor( 255, 255, 255 );
      m_inputImage( y, x )     = openMVG::image::RGBColor( 255, 255, 255 );
      m_inputImage( y + 1, x ) = openMVG::image::RGBColor( 255, 255, 255 );
    }
  }

  QImage tmp = openMVGImageToQImage( m_inputImage );
  m_scene_old->clear();
  m_scene_old->addPixmap( QPixmap::fromImage( tmp ) );
  m_view_old->fitInView( m_scene_old->itemsBoundingRect(), Qt::KeepAspectRatio );

  m_ppx->setValue( w / 2 );
  m_ppy->setValue( h / 2 );
  m_fov_spin->setValue( 1.1 * std::max( w, h ) );
}

/**
  * @brief update interface based on internal state 
  */
void MainWindow::UpdateInterface( void )
{
  switch ( m_current_state )
  {
  case STATE_INITIAL:
  {
    m_open_image->setEnabled( true );
    m_process->setEnabled( false );
    m_exportResult->setEnabled( false );

    if ( m_process_act )
    {
      m_process_act->setEnabled( false );
      m_export_act->setEnabled( false );
    }

    break;
  }
  case STATE_HAS_INPUT:
  {
    if ( m_check_use_grid->isChecked() )
    {
      m_open_image->setEnabled( false );
      m_open_image_act->setEnabled( false );
    }
    else
    {
      m_open_image->setEnabled( true );
      m_open_image_act->setEnabled( true );
    }
    m_exportResult->setEnabled( false );
    m_export_act->setEnabled( false );

    m_process_act->setEnabled( true );
    m_process->setEnabled( true );

    break;
  }
  case STATE_HAS_COMPUTED_DISTORTION:
  {
    m_process->setEnabled( true );
    m_exportResult->setEnabled( true );

    m_process_act->setEnabled( true );
    m_export_act->setEnabled( true );

    break;
  }
  }

  if ( !m_fileMenu || !m_processMenu || !m_inputMenu )
    return;

  switch ( m_intrinsic_type->currentIndex() )
  {
  case 0:
  {
    m_no_dist_act->setChecked( true );
    m_radial_k1_act->setChecked( false );
    m_radial_k3_act->setChecked( false );
    m_brown_act->setChecked( false );
    m_fisheye_act->setChecked( false );

    break;
  }
  case 1:
  {
    m_no_dist_act->setChecked( false );
    m_radial_k1_act->setChecked( true );
    m_radial_k3_act->setChecked( false );
    m_brown_act->setChecked( false );
    m_fisheye_act->setChecked( false );

    break;
  }
  case 2:
  {
    m_no_dist_act->setChecked( false );
    m_radial_k1_act->setChecked( false );
    m_radial_k3_act->setChecked( true );
    m_brown_act->setChecked( false );
    m_fisheye_act->setChecked( false );

    break;
  }
  case 3:
  {
    m_no_dist_act->setChecked( false );
    m_radial_k1_act->setChecked( false );
    m_radial_k3_act->setChecked( false );
    m_brown_act->setChecked( true );
    m_fisheye_act->setChecked( false );

    break;
  }
  case 4:
  {
    m_no_dist_act->setChecked( false );
    m_radial_k1_act->setChecked( false );
    m_radial_k3_act->setChecked( false );
    m_brown_act->setChecked( false );
    m_fisheye_act->setChecked( true );

    break;
  }
  }
}

/**
  * @brief Build Interface widgets 
  */
void MainWindow::BuildInterface()
{
  QGroupBox *inputGrp = new QGroupBox( "Input" );

  QGridLayout *inputLayout = new QGridLayout;

  m_open_image     = new QPushButton( "Open image" );
  m_check_use_grid = new QCheckBox( "Use grid" );

  inputLayout->addWidget( m_open_image, 0, 0 );
  inputLayout->addWidget( m_check_use_grid, 0, 1 );

  inputGrp->setLayout( inputLayout );

  QGroupBox *intrinsicGrp = new QGroupBox( "Intrinsic" );

  QGridLayout *intrinsicLayout = new QGridLayout;

  m_intrinsic_label = new QLabel( "Type" );
  m_intrinsic_type  = new QComboBox;

  QGroupBox *pinholeGrp = new QGroupBox( "Pinhole" );

  QGridLayout *pinholeLayout = new QGridLayout;

  m_label_fov = new QLabel( "Focal" );
  m_fov_spin  = new QDoubleSpinBox;
  m_fov_spin->setMaximum( 100000.0 );

  m_label_ppx = new QLabel( "Princ. P. X" );
  m_ppx       = new QDoubleSpinBox;
  m_ppx->setMaximum( 50000.0 );
  m_ppx->setMinimum( -50000.0 );
  m_label_ppy = new QLabel( "Princ. P. Y" );
  m_ppy       = new QDoubleSpinBox;
  m_ppy->setMaximum( 50000.0 );
  m_ppy->setMinimum( -50000.0 );

  pinholeLayout->addWidget( m_label_fov, 0, 0, Qt::AlignRight );
  pinholeLayout->addWidget( m_fov_spin, 0, 1 );
  pinholeLayout->addWidget( m_label_ppx, 0, 2, Qt::AlignRight );
  pinholeLayout->addWidget( m_ppx, 0, 3 );
  pinholeLayout->addWidget( m_label_ppy, 0, 4, Qt::AlignRight );
  pinholeLayout->addWidget( m_ppy, 0, 5 );

  pinholeGrp->setLayout( pinholeLayout );

  QGroupBox *distGrp = new QGroupBox( "Distorsion" );

  QGridLayout *distLayout = new QGridLayout;

  m_label_k1 = new QLabel( "K1" );
  m_k1       = new QDoubleSpinBox;
  m_k1->setDecimals( 6 );
  m_k1->setMaximum( 10000.0 );
  m_k1->setMinimum( -10000.0 );
  m_label_k2 = new QLabel( "K2" );
  m_k2       = new QDoubleSpinBox;
  m_k2->setDecimals( 6 );
  m_k2->setMaximum( 10000.0 );
  m_k2->setMinimum( -10000.0 );
  m_label_k3 = new QLabel( "K3" );
  m_k3       = new QDoubleSpinBox;
  m_k3->setDecimals( 6 );
  m_k3->setMaximum( 10000.0 );
  m_k3->setMinimum( -10000.0 );
  m_label_k4 = new QLabel( "K4" );
  m_k4       = new QDoubleSpinBox;
  m_k4->setDecimals( 6 );
  m_k4->setMaximum( 10000.0 );
  m_k4->setMinimum( -10000.0 );
  m_label_t1 = new QLabel( "T1" );
  m_t1       = new QDoubleSpinBox;
  m_t1->setDecimals( 6 );
  m_t1->setMaximum( 10000.0 );
  m_t1->setMinimum( -10000.0 );
  m_label_t2 = new QLabel( "T2" );
  m_t2       = new QDoubleSpinBox;
  m_t2->setDecimals( 6 );
  m_t2->setMaximum( 10000.0 );
  m_t2->setMinimum( -10000.0 );

  distLayout->addWidget( m_label_k1, 0, 0, Qt::AlignRight );
  distLayout->addWidget( m_k1, 0, 1 );
  distLayout->addWidget( m_label_k2, 0, 2, Qt::AlignRight );
  distLayout->addWidget( m_k2, 0, 3 );
  distLayout->addWidget( m_label_k3, 0, 4, Qt::AlignRight );
  distLayout->addWidget( m_k3, 0, 5 );
  distLayout->addWidget( m_label_k4, 0, 6, Qt::AlignRight );
  distLayout->addWidget( m_k4, 0, 7 );
  distLayout->addWidget( m_label_t1, 1, 0, Qt::AlignRight );
  distLayout->addWidget( m_t1, 1, 1 );
  distLayout->addWidget( m_label_t2, 1, 2, Qt::AlignRight );
  distLayout->addWidget( m_t2, 1, 3 );

  distGrp->setLayout( distLayout );

  intrinsicLayout->addWidget( m_intrinsic_label, 0, 0 );
  intrinsicLayout->addWidget( m_intrinsic_type, 0, 1 );

  intrinsicLayout->addWidget( pinholeGrp, 1, 0, 1, 2 );
  intrinsicLayout->addWidget( distGrp, 2, 0, 1, 2 );

  intrinsicGrp->setLayout( intrinsicLayout );

  // Processing
  QGroupBox *processGrp = new QGroupBox( "Processing" );

  QHBoxLayout *layoutProcess = new QHBoxLayout;

  m_process      = new QPushButton( "Undistort" );
  m_exportResult = new QPushButton( "Export" );

  layoutProcess->addWidget( m_process );
  layoutProcess->addWidget( m_exportResult );

  processGrp->setLayout( layoutProcess );

  // Result
  QGroupBox *resultGrp = new QGroupBox( "Images" );

  m_scene_new = new QGraphicsScene;
  m_scene_old = new QGraphicsScene;
  m_view_new  = new QGraphicsView( m_scene_new );
  m_view_old  = new QGraphicsView( m_scene_old );

  QHBoxLayout *resultLayout = new QHBoxLayout;

  resultLayout->addWidget( m_view_old );
  resultLayout->addWidget( m_view_new );

  QVBoxLayout *dummyResultLayout = new QVBoxLayout;

  dummyResultLayout->addLayout( resultLayout );

  resultGrp->setLayout( dummyResultLayout );

  QWidget *dummy           = new QWidget;
  QVBoxLayout *dummyLayout = new QVBoxLayout;

  dummyLayout->addWidget( inputGrp );
  dummyLayout->addWidget( intrinsicGrp );
  dummyLayout->addWidget( processGrp );
  dummyLayout->addWidget( resultGrp );

  dummy->setLayout( dummyLayout );

  setCentralWidget( dummy );

  FillIntrinsicTypes();
  onSelectIntrinsicType();
}

/**
  * @brief Create menus items 
  */
void MainWindow::BuildMenus()
{
  m_fileMenu = new QMenu( "File" );
  m_quit_act = m_fileMenu->addAction( "Quit" );
  m_quit_act->setShortcut( QKeySequence::Quit );
  m_fileMenu->addAction( m_quit_act );

  m_inputMenu      = new QMenu( "Input" );
  m_open_image_act = m_fileMenu->addAction( "Open image" );
  m_open_image_act->setShortcut( QKeySequence::Open );
  m_use_grid_act = m_fileMenu->addAction( "Use grid" );
  m_inputMenu->addAction( m_open_image_act );
  m_inputMenu->addAction( m_use_grid_act );

  m_intrinsicMenu = new QMenu( "Intrinsic Type" );
  m_no_dist_act   = m_fileMenu->addAction( "No distortion" );
  m_no_dist_act->setCheckable( true );
  m_radial_k1_act = m_fileMenu->addAction( "Radial 1" );
  m_radial_k1_act->setCheckable( true );
  m_radial_k3_act = m_fileMenu->addAction( "Radial 3" );
  m_radial_k3_act->setCheckable( true );
  m_brown_act = m_fileMenu->addAction( "Brown (3R-2T)" );
  m_brown_act->setCheckable( true );
  m_fisheye_act = m_fileMenu->addAction( "Fisheye" );
  m_fisheye_act->setCheckable( true );
  m_intrinsicMenu->addAction( m_no_dist_act );
  m_intrinsicMenu->addAction( m_radial_k1_act );
  m_intrinsicMenu->addAction( m_radial_k3_act );
  m_intrinsicMenu->addAction( m_brown_act );
  m_intrinsicMenu->addAction( m_fisheye_act );

  m_processMenu = new QMenu( "Processing" );
  m_process_act = m_fileMenu->addAction( "Undistort" );
  m_export_act  = m_fileMenu->addAction( "Export image" );
  m_processMenu->addAction( m_process_act );
  m_processMenu->addAction( m_export_act );

  QMenuBar *mbar = this->menuBar();
  mbar->addMenu( m_fileMenu );
  mbar->addMenu( m_inputMenu );
  mbar->addMenu( m_intrinsicMenu );
  mbar->addMenu( m_processMenu );
}

/**
  * @brief Make connections between interface items and slots 
  */
void MainWindow::MakeConnections()
{
  connect( m_open_image, SIGNAL( clicked() ), this, SLOT( onOpenImage() ) );
  connect( m_process, SIGNAL( clicked() ), this, SLOT( onProcessImage() ) );
  connect( m_exportResult, SIGNAL( clicked() ), this, SLOT( onExportImage() ) );
  connect( m_intrinsic_type, SIGNAL( currentIndexChanged( int ) ), this, SLOT( onSelectIntrinsicType() ) );
  connect( m_check_use_grid, SIGNAL( stateChanged( int ) ), this, SLOT( onChangeCheckBox() ) );

  connect( m_open_image_act, SIGNAL( triggered() ), this, SLOT( onOpenImage() ) );
  connect( m_process_act, SIGNAL( triggered() ), this, SLOT( onProcessImage() ) );
  connect( m_export_act, SIGNAL( triggered() ), this, SLOT( onExportImage() ) );
  connect( m_quit_act, SIGNAL( triggered() ), this, SLOT( onQuit() ) );
  connect( m_no_dist_act, SIGNAL( triggered() ), this, SLOT( onMenuSelectNoDist() ) );
  connect( m_radial_k1_act, SIGNAL( triggered() ), this, SLOT( onMenuSelectRadial1() ) );
  connect( m_radial_k3_act, SIGNAL( triggered() ), this, SLOT( onMenuSelectRadial3() ) );
  connect( m_brown_act, SIGNAL( triggered() ), this, SLOT( onMenuSelectBrown() ) );
  connect( m_fisheye_act, SIGNAL( triggered() ), this, SLOT( onMenuSelectFisheye() ) );
}

} // namespace image_undistort_gui
