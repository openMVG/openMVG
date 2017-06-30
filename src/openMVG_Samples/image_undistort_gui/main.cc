// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "MainWindow.hh"

#include <QApplication>

using namespace image_undistort_gui;

int main( int argc , char ** argv )
{
  QApplication app( argc , argv );

  MainWindow win;
  win.show();

  return app.exec();
}
