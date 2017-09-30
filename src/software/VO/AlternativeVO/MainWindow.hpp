// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ALTERNATIVE_VO_MAIN_WINDOW_HPP_
#define ALTERNATIVE_VO_MAIN_WINDOW_HPP_

#include "software/VO/AlternativeVO/ControlButtonsPanel.hpp"
#include "software/VO/AlternativeVO/VOViewerPanel.hpp"
#include "software/VO/AlternativeVO/VOFolderProcessor.hpp"

#include <QMainWindow>
#include <QMenu>
#include <QAction>
#include <QTimer>

#include <memory>

namespace alternative_vo
{

/**
* @brief Entry point of VO application
*/
class MainWindow : public QMainWindow
{
    Q_OBJECT

  public:
    /**
    * @brief default ctr
    */
    MainWindow( );

  public slots:

    /**
    * @brief Action to be executed when user open an image folder
    */
    void onMenuOpenImageFolder( void );

    /**
    * @brief Action to be executed when user want to quit the application
    * @note : Do nothing, users should love our app and there's no reason to quit us.
    */
    void onMenuQuit( void );

    /**
    * Action to be executed when user press play button
    */
    void onPlayVOProcess( void );

    /**
    * Action to be executed when user press stop/pause button
    */
    void onStopVOProcess( void );

    /**
    * Action to be executed when user press Forward button
    */
    void onStepForwardVOProcess( void );

    /**
    * Action to be executed when user press Reset button
    * @note Reset is when user want to start again the process
    */
    void onResetVOProcess( void );

    /**
    * Action to be executed when timer got a tick
    * -> Used to play a video (maybe not the best way to do it)
    */
    void onTimerTick( void );

  private:

    void BuildInterface();
    void BuildMenus();
    void MakeConnections();

    std::shared_ptr< VOFolderProcessor > m_vo_processor;

    QMenu * m_file_menu;
    QAction * m_folder_open_act;
    QAction * m_quit_act;

    QTimer * m_timer;

    VOViewerPanel * m_viewer_panel;
    ControlButtonsPanel * m_control_panel;

};

}

#endif
