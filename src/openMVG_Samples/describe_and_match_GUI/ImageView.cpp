// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ImageView.hpp"

#include <QWheelEvent>

#include <cmath>

namespace features_pair_demo
{

/**
  * @brief ctr
  * @param scn the scene to render in this widget
  * @param parent The parent widget
  */
ImageView::ImageView( QGraphicsScene *scn, QWidget *parent )
    : QGraphicsView( scn, parent )
{
  setRenderHints( QPainter::Antialiasing | QPainter::SmoothPixmapTransform );
}

/**
  * @brief handling of the zoom effect
  * @param event Container used to answer the mouse wheel informations
  */
void ImageView::wheelEvent( QWheelEvent *event )
{
  // Store current anchor
  const ViewportAnchor anchor = transformationAnchor();
  setTransformationAnchor( QGraphicsView::AnchorUnderMouse );
  const int angle    = event->angleDelta().y();
  const qreal factor = std::pow( 1.01, event->angleDelta().y() );

  scale( factor, factor );
  // Store back current anchor
  setTransformationAnchor( anchor );
}

} // namespace features_pair_demo
