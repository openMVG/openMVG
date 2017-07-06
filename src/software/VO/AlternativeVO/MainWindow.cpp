
// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "MainWindow.hpp"

#include <QVBoxLayout>
#include <QFileDialog>
#include <QMessageBox>
#include <QCoreApplication>
#include <QMenuBar>

namespace alternative_vo
{
/**
* @brief default ctr
*/
MainWindow::MainWindow( )
  : m_vo_processor( nullptr )
{
  BuildInterface();
  BuildMenus();
  MakeConnections();

  m_control_panel->EnableButtons( false );

  resize( 1024 , 768 );
}

/**
* @brief Action to be executed when user open an image folder
*/
void MainWindow::onMenuOpenImageFolder( void )
{
  QString folder = QFileDialog::getExistingDirectory( this , "Select Input image folder" , QDir::homePath() );
  if ( ! ( folder.isNull() || folder.isEmpty() ) )
  {
    m_vo_processor = std::make_shared< VOFolderProcessor >( folder.toStdString() );

    if (m_vo_processor->NbFrame() == 0 )
    {
      QMessageBox::critical( this , "Input Error" , "Provided folder does not contain any valid image file" );
      m_control_panel->EnableButtons( true );
      //      m_control_panel->EnablePlayBtn( false );
      m_vo_processor = nullptr;
    }
    else
    {
      m_control_panel->EnableButtons( true );
    }
  }
}

/**
* @brief Action to be executed when user want to quit the application
* @note : Do nothing, users should love our app and there's no reason to quit us.
*/
void MainWindow::onMenuQuit( void )
{
  QCoreApplication::quit();
}

/**
* Action to be executed when user press play button
*/
void MainWindow::onPlayVOProcess( void )
{
  if (m_timer )
  {
    m_timer->start( 40 ); // 30 ms -> 25 fps
  }
  else
  {
    m_timer = new QTimer;
    connect( m_timer , SIGNAL( timeout() ) , this , SLOT( onTimerTick() ) );
    m_timer->start( 40 );
  }
}

/**
* Action to be executed when user press stop/pause button
*/
void MainWindow::onStopVOProcess( void )
{
  if (m_timer )
  {
    m_timer->stop();
  }
}

/**
* Action to be executed when user press Forward button
*/
void MainWindow::onStepForwardVOProcess( void )
{
  if (m_vo_processor )
  {
    const size_t curr_ID = m_vo_processor->CurrentFrameID( );

    if (m_vo_processor->StepForward() )
    {
      // Has computed a new frame, get, update image and it's corresponding landmarks
      QImage img( m_vo_processor->FullFileName( curr_ID ).c_str() );

      // Get all points
      std::vector< VOViewerPoint > created = m_vo_processor->GetCreatedPoints();
      std::vector< VOViewerPoint > tracked = m_vo_processor->GetCurrentTrackedPoints();

      created.insert( created.end() , tracked.begin() , tracked.end() );

      // Get all lines
      const std::vector< VOViewerLine > trajs = m_vo_processor->GetCurrentTrackTrajectories();

      m_viewer_panel->SetImage( img , trajs , created );
    }
    else
    {
      // Nothing can be done more, stop processing
      m_control_panel->EnablePlayBtn( false );
      m_control_panel->EnableStopBtn( false );
      m_control_panel->EnableForwardBtn( false );
    }

    // Disable further update if we are at the end of the process
    if (m_vo_processor->CurrentFrameID() >= m_vo_processor->NbFrame() )
    {
      m_control_panel->EnablePlayBtn( false );
      m_control_panel->EnableStopBtn( false );
      m_control_panel->EnableForwardBtn( false );
    }
  }
}

/**
* Action to be executed when user press Reset button
* @note Reset is when user want to start again the process
*/
void MainWindow::onResetVOProcess( void )
{
  if (m_vo_processor )
  {
    m_vo_processor->Reset();

    m_control_panel->EnableButtons( true );
  }
}

void MainWindow::BuildInterface()
{
  QWidget * dummy = new QWidget;

  m_control_panel = new ControlButtonsPanel( this );
  m_viewer_panel = new VOViewerPanel( this );

  QVBoxLayout * mainLayout = new QVBoxLayout;

  mainLayout->addWidget( m_viewer_panel );
  mainLayout->addStretch( );
  mainLayout->addWidget( m_control_panel );


  dummy->setLayout( mainLayout );

  setCentralWidget( dummy );
}

void MainWindow::BuildMenus()
{
  m_file_menu = new QMenu( "File" );

  m_folder_open_act = new QAction( "Open image folder" , this );
  m_quit_act = new QAction( "Quit" , this );

  m_file_menu->addAction( m_folder_open_act );
  m_file_menu->addSeparator( );
  m_file_menu->addAction( m_quit_act );

  QMenuBar * menu_bar = menuBar();

  menu_bar->addMenu( m_file_menu );
}

void MainWindow::MakeConnections()
{
  // Menus elts
  connect( m_folder_open_act , SIGNAL( triggered() ) , this , SLOT( onMenuOpenImageFolder() ) );
  connect( m_quit_act , SIGNAL( triggered() ) , this , SLOT( onMenuQuit() ) );

  // Control panel
  connect( m_control_panel , SIGNAL( hasClickedPlay() ) , this , SLOT( onPlayVOProcess() ) );
  connect( m_control_panel , SIGNAL( hasClickedStop() ) , this , SLOT( onStopVOProcess() ) );
  connect( m_control_panel , SIGNAL( hasClickedReset() ) , this , SLOT( onResetVOProcess() ) );
  connect( m_control_panel , SIGNAL( hasClickedForward() ) , this , SLOT( onStepForwardVOProcess() ) );
}

void MainWindow::onTimerTick()
{
  if ( ! m_vo_processor )
  {
    if (m_timer )
    {
      m_timer->stop();
    }
  }
  else
  {
    onStepForwardVOProcess();
  }
}

}
