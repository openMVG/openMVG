
// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ALTERNATIVE_VO_VO_VIEWER_PANNEL_HPP_
#define ALTERNATIVE_VO_VO_VIEWER_PANNEL_HPP_

#include "software/VO/AlternativeVO/VOViewerDrawableElements.hpp"

#include "openMVG/numeric/numeric.h"

#include <QLabel>
#include <QWidget>
#include <QImage>

namespace alternative_vo
{


class VOViewerPanel : public QWidget
{
  public:
    VOViewerPanel( QWidget * parent );

    void SetImage( QImage & img );

    /**
    * @brief Set an image and draw some elements on it
    * @param lines List of lines to draw
    */
    void SetImage( QImage & img , const std::vector< VOViewerLine > & lines , const std::vector< VOViewerPoint > & points );

  private:

    void BuildInterface( void );
    void MakeConnections( void );

    QLabel * m_image_view;
};
}

#endif