// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/VO/AlternativeVO/VOViewerPanel.hpp"

#include <QHBoxLayout>
#include <QPainter>

namespace alternative_vo
{
VOViewerPanel::VOViewerPanel( QWidget * parent )
{
  BuildInterface();
  MakeConnections();
}

void VOViewerPanel::SetImage( QImage & img )
{
  m_image_view->setPixmap( QPixmap::fromImage( img ) );
}

/**
* @brief Set an image and draw some elements on it
* @param lines List of lines to draw
*/
void VOViewerPanel::SetImage( QImage & img , const std::vector< VOViewerLine > & lines , const std::vector< VOViewerPoint > & points )
{
  QPixmap pmap = QPixmap::fromImage( img );

  QPainter p( & pmap );

  // Draw lines
  for (size_t i = 0; i < lines.size(); ++i )
  {
    p.setPen( lines[i].m_color );

    const VOViewerLine & cur_line = lines[ i ];

    if (cur_line.m_pts.size() > 0 )
    {
      openMVG::Vec2f prev_pt = cur_line.m_pts[ 0 ];
      for (size_t id_point = 1; id_point < cur_line.m_pts.size(); ++id_point )
      {
        openMVG::Vec2f cur_pt = cur_line.m_pts[ id_point ];

        p.drawLine( QPoint( prev_pt( 0 ) , prev_pt( 1 ) ) , QPoint( cur_pt( 0 ) , cur_pt( 1 ) ) );

        prev_pt = cur_pt;
      }
    }


  }

  // Draw points
  for (size_t i = 0; i < points.size(); ++i )
  {
    const VOViewerPoint & cur_pt = points[ i ];
    p.setPen( cur_pt.m_color );
    p.drawPoint( QPoint( cur_pt.m_pt( 0 ) , cur_pt.m_pt( 1 ) ) );
  }

  m_image_view->setPixmap( pmap );

}


void VOViewerPanel::BuildInterface( void )
{
  QHBoxLayout * mainLayout = new QHBoxLayout;

  m_image_view = new QLabel( this );

  mainLayout->addWidget( m_image_view );

  setLayout( mainLayout );
}

void VOViewerPanel::MakeConnections( void )
{

}

}
