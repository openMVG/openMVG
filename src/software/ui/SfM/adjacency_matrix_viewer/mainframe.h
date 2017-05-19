// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <QFrame>
#include <QGraphicsView>
#include "mainwindow.h"

class MatchingPairGraphicsView;
class QSlider;
class QToolButton;

class MainFrame : public QFrame
{
  Q_OBJECT
public:
  explicit MainFrame
  (
    const QString &name,
    const Document & doc,
    QWidget *parent = nullptr
  );

  QGraphicsView *view() const;

public slots:
  void zoomIn(int level = 1);
  void zoomOut(int level = 1);
  void resetView();

private slots:
  void setResetButtonEnabled();
  void setupMatrix();
  void togglePointerMode();
  void toggleOpenGL();
  void toggleAntialiasing();

private:
  MatchingPairGraphicsView *matching_pair_view;

  QToolButton
    *select_mode_button,
    *opengl_button,
    *antialias_button,
    *reset_button;
  QSlider *zoom_slider;

  const Document & doc;
};
