
// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "MainWindow.hpp"

#include <QApplication>

using namespace alternative_vo;

int main( int argc , char ** argv )
{
  QApplication app( argc , argv );

  MainWindow win;
  win.show();

  return app.exec();
}
