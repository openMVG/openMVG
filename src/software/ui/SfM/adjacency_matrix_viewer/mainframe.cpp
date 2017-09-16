// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "mainframe.h"
#include "pairgraphicsitem.h"
#include "matchingpairgraphicsview.h"

#ifndef QT_NO_OPENGL
#include <QtOpenGL>
#else
#include <QtWidgets>
#endif
#include <QSvgWidget>

QGraphicsView *MainFrame::view() const
{
  return static_cast<QGraphicsView *>(matching_pair_view);
}

MainFrame::MainFrame(const QString &name, const Document & doc, QWidget *parent)
  : QFrame(parent), doc(doc)
{
  setFrameStyle(Sunken | StyledPanel);
  matching_pair_view = new MatchingPairGraphicsView(this, doc);
  matching_pair_view->setRenderHint(QPainter::Antialiasing, false);
  matching_pair_view->setDragMode(QGraphicsView::RubberBandDrag);
  matching_pair_view->setOptimizationFlags(QGraphicsView::DontSavePainterState);
  matching_pair_view->setViewportUpdateMode(QGraphicsView::SmartViewportUpdate);
  matching_pair_view->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

  // configure the icon size
  const int size = style()->pixelMetric(QStyle::PM_ToolBarIconSize);
  const QSize iconSize(size, size);

  QToolButton *zoom_in_icon = new QToolButton;
  zoom_in_icon->setText(tr("-"));
  zoom_in_icon->setIconSize(iconSize);
  QToolButton *zoom_out_icon = new QToolButton;
  zoom_out_icon->setText(tr("+"));
  zoom_out_icon->setIconSize(iconSize);
  zoom_slider = new QSlider;
  zoom_slider->setMinimum(-500);
  zoom_slider->setMaximum(500);
  zoom_slider->setValue(0);
  zoom_slider->setTickPosition(QSlider::TicksRight);

  // Zoom slider layout
  QVBoxLayout *zoom_slider_layout = new QVBoxLayout;
  zoom_slider_layout->addWidget(zoom_in_icon);
  zoom_slider_layout->addWidget(zoom_slider);
  zoom_slider_layout->addWidget(zoom_out_icon);

  reset_button = new QToolButton;
  reset_button->setText(tr("0"));
  reset_button->setEnabled(false);

  // Label layout
  QHBoxLayout *label_layout = new QHBoxLayout;
  QLabel *label = new QLabel(name);
  QLabel *label2 = new QLabel(tr("Pointer Mode"));
  QToolButton *select_mode_button = new QToolButton;
  select_mode_button->setText(tr("Select"));
  select_mode_button->setCheckable(true);
  select_mode_button->setChecked(true);
  QToolButton *drag_mode_button = new QToolButton;
  drag_mode_button->setText(tr("Drag"));
  drag_mode_button->setCheckable(true);
  drag_mode_button->setChecked(false);
  QToolButton * antialias_button = new QToolButton;
  antialias_button->setText(tr("Antialiasing"));
  antialias_button->setCheckable(true);
  antialias_button->setChecked(false);
  opengl_button = new QToolButton;
  opengl_button->setText(tr("OpenGL"));
  opengl_button->setCheckable(true);
#ifndef QT_NO_OPENGL
  opengl_button->setEnabled(QGLFormat::hasOpenGL());
#else
  opengl_button->setEnabled(false);
#endif

  //
  // Build the window layout
  //

  QButtonGroup *pointer_mode_group = new QButtonGroup(this);
  pointer_mode_group->setExclusive(true);
  pointer_mode_group->addButton(select_mode_button);
  pointer_mode_group->addButton(drag_mode_button);

  label_layout->addWidget(label);
  label_layout->addStretch();
  label_layout->addWidget(label2);
  label_layout->addWidget(select_mode_button);
  label_layout->addWidget(drag_mode_button);
  label_layout->addStretch();
  label_layout->addWidget(antialias_button);
  label_layout->addWidget(opengl_button);

  QGridLayout *top_layout = new QGridLayout;
  top_layout->addLayout(label_layout, 0, 0);
  top_layout->addWidget(matching_pair_view, 1, 0);
  top_layout->addLayout(zoom_slider_layout, 1, 1);
  top_layout->addWidget(reset_button, 2, 1);
  setLayout(top_layout);

  //
  // Confiugre connection between SIGNALs and SLOTs
  //

  connect(reset_button, SIGNAL(clicked()), this, SLOT(resetView()));
  connect(zoom_slider, SIGNAL(valueChanged(int)), this, SLOT(setupMatrix()));
  connect(matching_pair_view->verticalScrollBar(), SIGNAL(valueChanged(int)),
          this, SLOT(setResetButtonEnabled()));
  connect(matching_pair_view->horizontalScrollBar(), SIGNAL(valueChanged(int)),
          this, SLOT(setResetButtonEnabled()));
  connect(select_mode_button, SIGNAL(toggled(bool)), this, SLOT(togglePointerMode()));
  connect(drag_mode_button, SIGNAL(toggled(bool)), this, SLOT(togglePointerMode()));
  connect(antialias_button, SIGNAL(toggled(bool)), this, SLOT(toggleAntialiasing()));
  connect(opengl_button, SIGNAL(toggled(bool)), this, SLOT(toggleOpenGL()));
  connect(zoom_in_icon, SIGNAL(clicked()), this, SLOT(zoomIn()));
  connect(zoom_out_icon, SIGNAL(clicked()), this, SLOT(zoomOut()));

  setupMatrix();
}

void MainFrame::resetView()
{
  zoom_slider->setValue(0);
  matching_pair_view->ensureVisible(QRectF(0, 0, 0, 0));

  matching_pair_view->fitInView(
    matching_pair_view->scene()->itemsBoundingRect(),
    Qt::KeepAspectRatio);

  reset_button->setEnabled(false);
}

void MainFrame::setResetButtonEnabled()
{
  reset_button->setEnabled(true);
}

void MainFrame::setupMatrix()
{
  const qreal scale = qPow(qreal(2), (zoom_slider->value() - 250) / qreal(50));

  QMatrix matrix;
  matrix.scale(scale, scale);

  matching_pair_view->setMatrix(matrix);
  setResetButtonEnabled();
}

void MainFrame::togglePointerMode()
{
  matching_pair_view->setDragMode(select_mode_button->isChecked()
                                  ? QGraphicsView::RubberBandDrag
                                  : QGraphicsView::ScrollHandDrag);
  matching_pair_view->setInteractive(select_mode_button->isChecked());
}

void MainFrame::toggleOpenGL()
{
#ifndef QT_NO_OPENGL
  matching_pair_view->setViewport(
    opengl_button->isChecked() ?
      new QGLWidget(QGLFormat(QGL::SampleBuffers)) :
      new QWidget);
#endif
}

void MainFrame::toggleAntialiasing()
{
  matching_pair_view->setRenderHint(QPainter::Antialiasing, antialias_button->isChecked());
}

void MainFrame::zoomIn(int level)
{
  zoom_slider->setValue(zoom_slider->value() + level);
}

void MainFrame::zoomOut(int level)
{
  zoom_slider->setValue(zoom_slider->value() - level);
}
