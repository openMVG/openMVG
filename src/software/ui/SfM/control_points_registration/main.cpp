// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <QApplication>
#include "mainLayout.hpp"

int main(int argc, char ** argv)
{
  QApplication app(argc, argv);

  MainWindow * mainWindow = new MainWindow;
  mainWindow->show();
  
  return app.exec();
}
