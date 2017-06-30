// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "pairgraphicsitem.h"
#include "mainwindow.h"
#include "mainframe.h"

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"

#include <QHBoxLayout>
#include <QSplitter>
#include <QtWidgets>

MainWindow::MainWindow(QWidget *parent)
  : QWidget(parent)
{
  view = new MainFrame("Pairwise matches viewer", doc);
  QMenuBar * menu_bar = new QMenuBar;
  QMenu * fileMenu = new QMenu(tr("&File"));
  menu_bar->addMenu(fileMenu);
  QAction * openAct = new QAction(tr("&Open..."), this);
  openAct->setShortcuts(QKeySequence::Open);
  openAct->setStatusTip(tr("Open an existing project"));
  connect(openAct, &QAction::triggered, this, &MainWindow::open);
  fileMenu->addAction(openAct);

  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(menu_bar);
  layout->addWidget(view);
  setLayout(layout);

  setWindowTitle(tr("Pairwise matches viewer"));
}

void MainWindow::open()
{
  const QString sfm_data_fileName = QFileDialog::getOpenFileName(
    this, tr("Choose a sfm_data project file"),
    QString::null, tr("sfm_data files (*.json *.xml *.bin)"));
  if (sfm_data_fileName.isEmpty())
    return;

  const std::string m_sfm_data_filename = sfm_data_fileName.toStdString();

  const QString matches_fileName = QFileDialog::getOpenFileName(
    this, tr("Choose a matches.X file"),
    QFileInfo(sfm_data_fileName).path(), tr("matches files (*.bin *.txt)"));
  if (matches_fileName.isEmpty())
    return;

  const std::string m_matches_data_filename = matches_fileName.toStdString();

  populateScene(m_sfm_data_filename,
    QFileInfo(sfm_data_fileName).path().toStdString(),
    m_matches_data_filename);
  view->view()->setScene(scene);

  emit view->resetView();
}

void MainWindow::populateScene
(
  const std::string & sSfM_Data_Filename,
  const std::string & sMatchesDir,
  const std::string & sMatchFile
)
{
  scene = new QGraphicsScene(this);

  using namespace openMVG;
  using namespace openMVG::features;
  using namespace openMVG::matching;
  using namespace openMVG::sfm;

  if (!Load(doc.sfm_data, sSfM_Data_Filename,
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename
      << "\" cannot be read." << std::endl;
    return;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type used for this scene from the image describer file
  const std::string sImage_describer =
    stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type =
    Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return;
  }

  // Read the features
  doc.feats_provider = std::make_shared<Features_Provider>();
  if (!doc.feats_provider->load(doc.sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return;
  }

  // Read the matches
  doc.matches_provider = std::make_shared<Matches_Provider>();
  if (!doc.matches_provider->load(doc.sfm_data, sMatchFile)) {
    std::cerr << "\nInvalid matches file." << std::endl;
    return;
  }

  std::cout << "Read #Pair: " << doc.matches_provider->getPairs().size()
    << std::endl;

  const Pair_Set pairs = doc.matches_provider->getPairs();

  for (const auto & pair_iter : pairs)
  {
    const auto
      I = pair_iter.first,
      J = pair_iter.second;

    const QColor color(0, 0, 255, 127);
    QGraphicsItem *item =
      new PairGraphicsItem(
        color, I , J,
        doc.matches_provider->pairWise_matches_.at(pair_iter).size());

    item->setPos(
      QPointF(
        J * item->boundingRect().width() * 1.1,
        I * item->boundingRect().height() * 1.1));
    scene->addItem(item);
  }
}
