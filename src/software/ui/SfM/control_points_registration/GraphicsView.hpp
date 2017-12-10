// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef GRAPHICSVIEW_HPP
#define GRAPHICSVIEW_HPP

// Qt4 headers
#include <QMainWindow>
#include <QPointer>
#include <QGraphicsView>

#include "document.hpp"

class QAction;
class QGraphicsPixmapItem;
class QGraphicsScene;
class QLabel;
class QScrollArea;
class QMenu;
class QGraphicsLineItem;
class QGraphicsRectItem;
class QGraphicsItem;

namespace control_point_GUI {

  class GraphicsView : public QGraphicsView
  {
    Q_OBJECT // mandatory for signals and slots

  private: /* static const variables */

  public: /* methods */
    GraphicsView(Document & doc, QWidget * parent = 0);

    void addImage(const QString & qs_filename, float xpos=0.f, float ypos=0.f, bool bClear = false);

    void addNode(QGraphicsItem* it);

    void setCurrentViewId(const openMVG::IndexT index) {_current_view_id = index;}

  protected: /* methods */
    void drawBackground(QPainter * painter, const QRectF &rect);

    void mousePressEvent(QMouseEvent* e );
    void wheelEvent(QWheelEvent * event );

    void zoom(qreal factor, QPointF centerPoint);

  private slots:
    void removeControlPoint();
    void zoomIn();
    void zoomOut();
    void normalSize();

  private: /* methods */
    // Interface construction methods
    void createActions();

  private: /* data members */
    // The graphics view machinery.
    QGraphicsScene * scene;

    // Action
    QAction * open_images_action_;
    QAction * zoomInAct;
    QAction * zoomOutAct;
    QAction * normalSizeAct;
    QAction * removeControlPointAct;

    // Document
    Document & _doc;
    // Current viewed image id
    openMVG::IndexT _current_view_id;
  };

} // namespace control_point_GUI

#endif /* GRAPHICSVIEW_HPP */
