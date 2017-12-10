// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GraphicsView.hpp"
#include "ControlPoint2DNode.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

// Qt4 headers
#include <QtGui>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsLineItem>
#include <QGraphicsItem>
#include <QAction>
#include <QMessageBox>
#include <memory>
#include <QPen>
#include <QInputDialog>

#include <iostream>
#include <fstream>
#include <iterator>

namespace control_point_GUI
{
  using namespace std;
  using namespace openMVG;
  using namespace openMVG::sfm;

  // =========================================================================
  // Public methods
  // =========================================================================
  GraphicsView::GraphicsView(Document & doc, QWidget * parent)
    : QGraphicsView(parent), scene(new QGraphicsScene),
    _doc(doc), _current_view_id(UndefinedIndexT)
  {
    setScene(scene);

    // The OpenGL rendering does not seem to work with too big images.
    //setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));
    setBackgroundRole(QPalette::Dark);
    setAlignment(Qt::AlignCenter);
    setCacheMode(QGraphicsView::CacheBackground);
    setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

    setMouseTracking(true);

    createActions();

    resize(800, 600);
  }

  void GraphicsView::zoomIn()
  {
    scale(1.25, 1.25); update();
  }

  void GraphicsView::zoomOut() { scale(0.8, .8); update(); }
  void GraphicsView::normalSize() { resetTransform();  update(); }

  void GraphicsView::removeControlPoint()
  {
    auto items = scene->selectedItems();

    foreach (auto item, items)
    {
      auto cp = dynamic_cast<ControlPoint2DNode *>(item);

      if (!cp)
        continue;

      auto i = _doc._sfm_data.control_points.find(cp->controlPointId());

      if (i != _doc._sfm_data.control_points.end())
      {
        Landmark & landmark = i->second;

        landmark.obs.erase(_current_view_id);

        if (landmark.obs.empty())
          _doc._sfm_data.control_points.erase(i);
      }

      delete cp;
    }
  }

  // =========================================================================
  // Protected methods
  // =========================================================================
  void GraphicsView::drawBackground(QPainter *painter,
    const QRectF &rect)
  {
  }

  void GraphicsView::mousePressEvent (QMouseEvent* e )
  {
    QGraphicsView::mousePressEvent(e);

    if (e->isAccepted())
      return; // QGraphicsView handled this event

    if (_doc._sfm_data.GetViews().empty() || e->button()!=Qt::LeftButton)
    {
      return;
    }

    e->accept(); // We handled this event

    const QPointF pos =  this->mapToScene(e->pos());

    int index = -1;
    QInputDialog input;
    input.setInputMode(QInputDialog::IntInput);
    input.setWindowTitle("");
    input.setLabelText("Control point ID");
    input.setIntRange(0, 9999);
    input.setIntValue(0);
    input.setIntStep(1);
    if (input.exec() == QDialog::Accepted)
    {
      index = input.intValue();
    }

    if (index != -1)
    {
      // Try to setup a new landmark:
      Landmark & landmark = this->_doc._sfm_data.control_points[index];
      if (landmark.obs.count(_current_view_id) == 0)
      {
        landmark.obs[_current_view_id] = Observation(Vec2(pos.x(), pos.y()), 0);
        // Create a graphical instance that points directly to the control point observation
        Vec2 & ref = landmark.obs[_current_view_id].x;
        this->addNode(new ControlPoint2DNode(pos, ref(0), ref(1), index));
      }
      else
      {
        QMessageBox msgBox;
        msgBox.setText("There is already an observation linked to the depicted control_point in this image. Please select a new Id.");
        msgBox.exec();

        // If the control point was not existing, clear it
        if (landmark.obs.empty())
          this->_doc._sfm_data.control_points.erase(this->_doc._sfm_data.control_points.find(index));
      }
    }
  }

  void GraphicsView::wheelEvent ( QWheelEvent * e)
  {
    if (e->modifiers().testFlag(Qt::ControlModifier))
    // zoom only when CTRL key is pressed
    {
      const int numSteps = e->delta() / 15 / 8;
      if (numSteps == 0) {
        e->ignore();
        return;
      }
      if (numSteps > 0)
        zoomIn();
      else
        zoomOut();

      e->accept();
    }
    else
    {
      QGraphicsView::wheelEvent(e);
    }
  }

  void GraphicsView::zoom(qreal factor, QPointF centerPoint)
  {
    scale(factor, factor);
    //centerOn(centerPoint);
  }

  // =========================================================================
  // Private methods
  // =========================================================================
  void GraphicsView::createActions()
  {
    zoomInAct = new QAction(tr("Zoom &In (25%)"), this);
    zoomInAct->setShortcut(tr("Ctrl++"));
    zoomInAct->setEnabled(false);
    connect(zoomInAct, SIGNAL(triggered()), this, SLOT(zoomIn()));
    addAction(zoomInAct);

    zoomOutAct = new QAction(tr("Zoom &Out (25%)"), this);
    zoomOutAct->setShortcut(tr("Ctrl+-"));
    zoomOutAct->setEnabled(false);
    connect(zoomOutAct, SIGNAL(triggered()), this, SLOT(zoomOut()));
    addAction(zoomOutAct);

    normalSizeAct = new QAction(tr("&Normal Size"), this);
    normalSizeAct->setShortcut(tr("Ctrl+1"));
    normalSizeAct->setEnabled(false);
    connect(normalSizeAct, SIGNAL(triggered()), this, SLOT(normalSize()));
    addAction(normalSizeAct);

    removeControlPointAct = new QAction(tr("&Remove Control Point"), this);
    removeControlPointAct->setShortcut(tr("Del"));
    removeControlPointAct->setEnabled(false);
    connect(removeControlPointAct, SIGNAL(triggered()), this, SLOT(removeControlPoint()));
    addAction(removeControlPointAct);
  }

  void GraphicsView::addImage(const QString & qs_filename, float xpos, float ypos, bool bClear)
  {
    if (bClear) {
      scene->clear();
    }

    QGraphicsPixmapItem * image = new QGraphicsPixmapItem;
    image->setTransformationMode(Qt::FastTransformation);
    scene->addItem(image);

    QPixmap pixmap(qs_filename);
    if (pixmap.isNull()) {
      QMessageBox::information(this, QString::null,
        tr("Cannot load QPixmap %1.").arg(qs_filename));
      return;
    }
    else {
      const QPointF offset = QPointF(xpos, ypos);
      image->setPixmap(pixmap);
      image->setPos(offset);
      zoomInAct->setEnabled(true);
      zoomOutAct->setEnabled(true);
      normalSizeAct->setEnabled(true);
      removeControlPointAct->setEnabled(true);

      const QFileInfo fi(qs_filename);
      const QString name = fi.fileName();
      qDebug() << "Current viewed image: " << name;
    }
  }

  void GraphicsView::addNode(QGraphicsItem* it)
  {
      scene->addItem(it);
  }
} // namespace control_point_GUI
