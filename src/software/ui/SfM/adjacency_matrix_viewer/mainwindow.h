// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <QWidget>
#include <document.hpp>

class QGraphicsScene;
class MainFrame;

class MainWindow : public QWidget
{
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);

private slots:
  void open();

private:

  // Open the openMVG data and create the QT item to visualize the pair
  void populateScene
  (
    const std::string & sSfM_Data_Filename,
    const std::string & sMatchesDir,
    const std::string & sMatchFile
  );

  QGraphicsScene * scene;
  MainFrame * view;

  Document doc;
};
