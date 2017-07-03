
// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/VO/AlternativeVO/ControlButtonsPanel.hpp"

#include <QHBoxLayout>

namespace alternative_vo
{
ControlButtonsPanel::ControlButtonsPanel( QWidget * parent )
  : QWidget( parent )
{
  BuildInterface();
  MakeConnections();
}

/**
* @brief Enable/disable all buttons
* @param enabled True to enable, false to disable
*/
void ControlButtonsPanel::EnableButtons( const bool enabled )
{
  m_btn_forward->setEnabled( enabled );
  m_btn_play->setEnabled( enabled );
  m_btn_stop->setEnabled( enabled );
  m_btn_reset->setEnabled( enabled );
}

/**
* @brief Enable/disable play button
* @param enabled True to enable, false to disable
*/
void ControlButtonsPanel::EnablePlayBtn( const bool enabled )
{
  m_btn_play->setEnabled( enabled );
}

/**
* @brief Enable/disable stop button
* @param enabled True to enable, false to disable
*/
void ControlButtonsPanel::EnableStopBtn( const bool enabled )
{
  m_btn_stop->setEnabled( enabled );
}

/**
* @brief Enable/disable forward button
* @param enabled True to enable, false to disable
*/
void ControlButtonsPanel::EnableForwardBtn( const bool enabled )
{
  m_btn_forward->setEnabled( enabled );
}

/**
* @brief Enable/disable reset button
* @param enabled True to enable, false to disable
*/
void ControlButtonsPanel::EnableResetBtn( const bool enabled )
{
  m_btn_reset->setEnabled( enabled );
}

// Action executed when user click play button
void ControlButtonsPanel::onClickPlayBtn( void )
{
  emit hasClickedPlay();
}

// Action executed when user click stop button
void ControlButtonsPanel::onClickStopBtn( void )
{
  emit hasClickedStop();
}

// Action executed when user click forward button
void ControlButtonsPanel::onClickForwardBtn( void )
{
  emit hasClickedForward();
}

// Action executed when user click reset button
void ControlButtonsPanel::onClickResetBtn( void )
{
  emit hasClickedReset();
}

void ControlButtonsPanel::BuildInterface( void )
{
  m_btn_reset = new QPushButton( "Reset" , this );
  m_btn_stop = new QPushButton( "Stop" , this );
  m_btn_play = new QPushButton( "Play" , this );
  m_btn_forward = new QPushButton( "Forward" , this );

  QHBoxLayout * mainLayout = new QHBoxLayout;

  mainLayout->addWidget( m_btn_reset );
  mainLayout->addWidget( m_btn_stop );
  mainLayout->addWidget( m_btn_play );
  mainLayout->addWidget( m_btn_forward );

  setLayout( mainLayout );
}

void ControlButtonsPanel::MakeConnections( void )
{
  connect( m_btn_forward , SIGNAL( clicked() ) , this , SLOT( onClickForwardBtn() ) );
  connect( m_btn_play , SIGNAL( clicked() ) , this , SLOT( onClickPlayBtn() ) );
  connect( m_btn_reset , SIGNAL( clicked() ) , this , SLOT( onClickResetBtn() ) );
  connect( m_btn_stop , SIGNAL( clicked() ) , this , SLOT( onClickStopBtn() ) );
}
}
