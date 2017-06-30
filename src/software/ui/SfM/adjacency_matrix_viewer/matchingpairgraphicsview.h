// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <QGraphicsView>
//class MainFrame;
#include "mainframe.h"
class Document;

class MatchingPairGraphicsView : public QGraphicsView
{
  Q_OBJECT
public:
  MatchingPairGraphicsView
  (
    MainFrame *v,
    const Document & doc
  );

protected:
#ifndef QT_NO_WHEELEVENT
  void wheelEvent(QWheelEvent *) Q_DECL_OVERRIDE;
#endif

  void mousePressEvent(QMouseEvent * event) Q_DECL_OVERRIDE;

private:
  MainFrame *main_frame;
  const Document &doc;
};
