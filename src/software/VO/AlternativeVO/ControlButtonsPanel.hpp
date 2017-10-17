// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ALTERNATIVE_VO_CONTROL_BUTTONS_PANEL_HPP_
#define ALTERNATIVE_VO_CONTROL_BUTTONS_PANEL_HPP_

#include <QWidget>
#include <QPushButton>

namespace alternative_vo
{

class ControlButtonsPanel : public QWidget
{
    Q_OBJECT

  public:
    ControlButtonsPanel( QWidget * parent );

    /**
    * @brief Enable/disable all buttons
    * @param enabled True to enable, false to disable
    */
    void EnableButtons( const bool enabled );

    /**
    * @brief Enable/disable play button
    * @param enabled True to enable, false to disable
    */
    void EnablePlayBtn( const bool enabled );

    /**
    * @brief Enable/disable stop button
    * @param enabled True to enable, false to disable
    */
    void EnableStopBtn( const bool enabled );
    /**
    * @brief Enable/disable forward button
    * @param enabled True to enable, false to disable
    */
    void EnableForwardBtn( const bool enabled );
    /**
    * @brief Enable/disable reset button
    * @param enabled True to enable, false to disable
    */
    void EnableResetBtn( const bool enabled );

  public slots:

    // Action executed when user click play button
    void onClickPlayBtn( void );
    // Action executed when user click stop button
    void onClickStopBtn( void );
    // Action executed when user click forward button
    void onClickForwardBtn( void );
    // Action executed when user click reset button
    void onClickResetBtn( void );

  signals:

    // Signal emitted after user click the play button
    void hasClickedPlay( void );
    // Signal emitted after user click the stop button
    void hasClickedStop( void );
    // Signal emitted after user click the forward button
    void hasClickedForward( void );
    // Signal emitted after user click the reset button
    void hasClickedReset( void );

  private:


    void BuildInterface( void );
    void MakeConnections( void );

    QPushButton * m_btn_reset;
    QPushButton * m_btn_stop;
    QPushButton * m_btn_play;
    QPushButton * m_btn_forward;
};

}

#endif
