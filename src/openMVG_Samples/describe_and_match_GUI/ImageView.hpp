// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SAMPLES_FEATURES_PAIR_DEMO_IMAGE_VIEW_HPP
#define OPENMVG_SAMPLES_FEATURES_PAIR_DEMO_IMAGE_VIEW_HPP

#include <QGraphicsView>

namespace features_pair_demo
{

// Subclass of qgraphicsview used to provide additional fonctionnalities (zooming with mouse wheel)
class ImageView : public QGraphicsView
{
public:
  /**
  * @brief ctr
  * @param scn the scene to render in this widget
  * @param parent The parent widget
  */
  ImageView( QGraphicsScene *scn, QWidget *parent = nullptr );

  /**
  * @brief handling of the zoom effect
  * @param event Container used to answer the mouse wheel informations
  */
  virtual void wheelEvent( QWheelEvent *event );

private:
};

} // namespace features_pair_demo

#endif
