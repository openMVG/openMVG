// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GraphicsMainWindow.hpp"
#include "node.hpp"

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
  
  // =========================================================================
  // Public methods
  // =========================================================================
  GraphicsView::GraphicsView(QGraphicsScene * scene)
    : QGraphicsView(scene)
  {
  }


  GraphicsMainWindow::GraphicsMainWindow(Document & doc, QWidget * parent)
    : QMainWindow(parent), scene(new QGraphicsScene),
    view(new GraphicsView(scene)), _doc(doc), _current_view_id(UndefinedIndexT)
  {
    // The OpenGL rendering does not seem to work with too big images.
    //view->setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));
    view->setBackgroundRole(QPalette::Dark);
    view->setAlignment(Qt::AlignCenter);
    view->setCacheMode(QGraphicsView::CacheBackground);
    view->setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);
    view->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

    setCentralWidget(view);

    view->setMouseTracking(true);

    createActions();

    resize(800, 600);
  }

  void GraphicsView::zoomIn()
  {
    scale(1.25, 1.25); update();
  }

  void GraphicsView::zoomOut() { scale(0.8, .8); update(); }
  void GraphicsView::normalSize() { resetTransform();  update(); }

  void GraphicsMainWindow::removeControlPoint()
  {
    auto items = scene->selectedItems();

    foreach (auto item, items)
    {
      auto cp = dynamic_cast<control_point_2DNode *>(item);

      if (!cp)
        continue;

      auto i = _doc._sfm_data.control_points.find(cp->id_control_point());

      if (i!=_doc._sfm_data.control_points.end())
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
  void GraphicsMainWindow::drawBackground(QPainter *painter,
    const QRectF &rect)
  {
  }

  void GraphicsMainWindow::mousePressEvent (QMouseEvent* e )
  {
    if (_doc._sfm_data.GetViews().empty())
    {
      return;
    }

    const QPointF pos =  this->view->mapToScene(e->pos());

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
        this->AddNode(new control_point_2DNode(pos, ref(0), ref(1), index));
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
  void GraphicsMainWindow::createActions()
  {
    zoomInAct = new QAction(tr("Zoom &In (25%)"), this);
    zoomInAct->setShortcut(tr("Ctrl++"));
    zoomInAct->setEnabled(false);
    connect(zoomInAct, SIGNAL(triggered()), view, SLOT(zoomIn()));
    addAction(zoomInAct);

    zoomOutAct = new QAction(tr("Zoom &Out (25%)"), this);
    zoomOutAct->setShortcut(tr("Ctrl+-"));
    zoomOutAct->setEnabled(false);
    connect(zoomOutAct, SIGNAL(triggered()), view, SLOT(zoomOut()));
    addAction(zoomOutAct);

    normalSizeAct = new QAction(tr("&Normal Size"), this);
    normalSizeAct->setShortcut(tr("Ctrl+S"));
    normalSizeAct->setEnabled(false);
    connect(normalSizeAct, SIGNAL(triggered()), view, SLOT(normalSize()));
    addAction(normalSizeAct);

    removeControlPointAct = new QAction(tr("&Remove Control Point"), this);
    removeControlPointAct->setShortcut(tr("Del"));
    connect(removeControlPointAct, SIGNAL(triggered()), this, SLOT(removeControlPoint()));
    addAction(removeControlPointAct);
  }

  void GraphicsMainWindow::AddImage(const QString & qs_filename, float xpos, float ypos, bool bClear)
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

      const QFileInfo fi(qs_filename);
      const QString name = fi.fileName();
      qDebug() << "Current viewed image: " << name;
    }
  }

  void GraphicsMainWindow::AddNode(QGraphicsItem* it)
  {
      scene->addItem(it);
  }
} // namespace control_point_GUI

