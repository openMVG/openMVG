// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "mainLayout.hpp"
#include <QtGui>
#include <QDebug>

#include <algorithm>
#include <clocale>

#include "node.hpp"
#include "ControlPointTableView.hpp"

#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"

#include "openMVG/stl/stl.hpp"

void MainWindow::removeAllControlPoints()
{
  m_doc._sfm_data.control_points = Landmarks();
  doubleClickImageList();
}

void MainWindow::help()
{
   QMessageBox::about(this, tr("About"),
   tr(
   "This application allow to perform a rigid registration of a SfM scene to know GCP (Ground Control Points)<br>"
   "<b>Important:</b><br>"
   "Each GCP must have an unique numbered ID.<br>"
   "<b>1-</b> Add GCP image observations:<br>"
   " - double-click on one image on the left (in which you want add a GCP observation)<br>"
   " - click on the displayed image on the right, and choose your ID and move your GCP observation to the right place.<br>"
   " - to move a GCP observation, click on it and move your mouse with left-click pressed.<br>"
   "<b>2-</b> Edit GCP 3D positions:<br>"
   "  - Go to Control Point Edition menu (Ctrl+E) and edit the 3D position<br>"
   "<b>3-</b> Registration:<br>"
   "- Go to Control Point registration menu (Ctrl+R) and save your scene.<br>"
   "<br>"
   "Tip:<br>"
   "CTRL+WHEEL: Image zoom in/out"
   ));
}

void MainWindow::doubleClickImageList()
{
  const QModelIndex modelIndex = m_treeView_Images->currentIndex();
  const QModelIndex res = modelIndex.sibling ( modelIndex.row(), 0 );
  const QString current_string_clicked = res.data().toString();

  std::istringstream is(current_string_clicked.toStdString());
  size_t id;
  is >> id;
  if (m_doc._sfm_data.GetViews().count(id) == 1)
  {
    const View * view = m_doc._sfm_data.GetViews().at(id).get();
    const std::string sView = stlplus::create_filespec(m_doc._sfm_data.s_root_path, view->s_Img_path);
    m_widget->AddImage(QString::fromStdString(sView), 0.f, 0.f, true);
    m_widget->SetCurrentViewId(view->id_view);

    setWindowTitle(QString("Control_point_editor: " + QString::fromStdString(sView)));

    // Display control points information related to this view
    for (Landmarks::iterator iterL = m_doc._sfm_data.control_points.begin();
      iterL != this->m_doc._sfm_data.control_points.end(); ++iterL)
    {
      Landmark & landmark = iterL->second;
      Observations & obs = landmark.obs;
      if (obs.count(view->id_view) != 0)
      {
        // Create a graphical instance that points directly to the control point 2d observation
        // It will allow dynamic update of the control_point observation when the user will move the Node
        Vec2 & ref = obs[view->id_view].x;
        const size_t index = iterL->first;
        m_widget->AddNode(new control_point_2DNode(QPointF(ref(0), ref(1)), ref(0), ref(1), index));
      }
    }
  }
}

void MainWindow::saveProject()
{
  const QString dir = QString::fromStdString(stlplus::folder_part(m_sfm_data_filename));
  const QString sfm_data_filename = QFileDialog::getSaveFileName(this, tr("Choose a sfm_data file (sfm_data.json)"),
    dir, tr("sfm_data (*.json *.xml *.bin)"));
  if (sfm_data_filename.isEmpty())
    return;

  std::setlocale(LC_ALL, "C");
  if (!m_doc.saveData(sfm_data_filename.toStdString()))
  {
    QMessageBox msgBox;
    msgBox.setText("Cannot save the sfm_data file.");
    msgBox.exec();
  }
}

void MainWindow::openProject()
{
  const QString sfm_data_fileName = QFileDialog::getOpenFileName(this, tr("Choose a sfm_data project file"),
    QString::null, tr("sfm_data (*.json *.xml *.bin)"));
  if (sfm_data_fileName.isEmpty())
    return;

  m_sfm_data_filename = sfm_data_fileName.toStdString();

  if (m_doc.loadData(sfm_data_fileName.toStdString()))
  {
    //Add image names in the QT tree view
    {
      QStandardItemModel * model = new QStandardItemModel(0,1, this);
      model->setHeaderData(0, Qt::Horizontal, QObject::tr("Views"));
      m_treeView_Images->setModel(model);

      std::vector<IndexT> view_ids;
      view_ids.reserve(m_doc._sfm_data.GetViews().size());
      std::transform(m_doc._sfm_data.GetViews().begin(), m_doc._sfm_data.GetViews().end(),
        std::back_inserter(view_ids), stl::RetrieveKey());
      std::sort(view_ids.begin(), view_ids.end());

      // Add view in reverse order to have them ordered by ID
      for (std::vector<IndexT>::const_reverse_iterator iter = view_ids.rbegin();
        iter != view_ids.rend(); ++iter)
      {
        Views::const_iterator iterV = m_doc._sfm_data.GetViews().find(*iter);
        const View * view = iterV->second.get();
        if (m_doc._sfm_data.IsPoseAndIntrinsicDefined(view))
        {
          std::ostringstream os;
          os << view->id_view << " " << view->s_Img_path;
          model->insertRow(0);
          model->setData(model->index(0, 0), QString::fromStdString(os.str()));
        }
      }
    }
  }
  else
  {
    QMessageBox msgBox;
    msgBox.setText("Cannot open the provided sfm_data file.");
    msgBox.exec();
  }
}

void MainWindow::editControlPoints()
{
  // Graphical widget to configure the control point position
  if (!m_doc._sfm_data.control_points.empty())
  {
    using namespace control_point_GUI;
    ControlPointTableView control_point_editor(&m_doc._sfm_data, this);
    control_point_editor.exec();
    Landmarks control_points_cpy = m_doc._sfm_data.control_points;
    control_point_editor.update_control_points(control_points_cpy);
    m_doc._sfm_data.control_points = control_points_cpy;
  }
  else
  {
    QMessageBox msgBox;
    msgBox.setText("No defined control points.");
    msgBox.exec();
    return;
  }
}

// Find pattern occurences and replace it by a new one
static std::string string_pattern_replace
(
  std::string subject,
  const std::string& search,
  const std::string& replace
)
{
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != std::string::npos)
  {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
  return subject;
}

void MainWindow::registerProject()
{
  if (m_doc._sfm_data.control_points.size() < 3)
  {
    QMessageBox msgBox;
    msgBox.setText("At least 3 control points are required.");
    msgBox.exec();
    return;
  }

  // Assert that control points can be triangulated
  for (Landmarks::const_iterator iterL = m_doc._sfm_data.control_points.begin();
    iterL != m_doc._sfm_data.control_points.end(); ++iterL)
  {
    if (iterL->second.obs.size() < 2)
    {
      QMessageBox msgBox;
      msgBox.setText("Each control point must be defined in at least 2 pictures.");
      msgBox.exec();
      return;
    }
  }

  //---
  // registration (coarse):
  // - compute the 3D points corresponding to the control point observation for the SfM scene
  // - compute a coarse registration between the controls points & the triangulated point
  // - transform the scene according the found transformation
  //---
  std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
  std::map<IndexT, double> vec_triangulation_errors;
  for (Landmarks::iterator iterCP = m_doc._sfm_data.control_points.begin();
    iterCP != m_doc._sfm_data.control_points.end(); ++iterCP)
  {
    Landmark & landmark = iterCP->second;
    //Triangulate the point:
    Triangulation trianObj;
    const Observations & obs = landmark.obs;
    for(Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = m_doc._sfm_data.views.at(itObs->first).get();
      if (!m_doc._sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      const openMVG::cameras::IntrinsicBase * cam = m_doc._sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const openMVG::geometry::Pose3 pose = m_doc._sfm_data.GetPoseOrDie(view);
      trianObj.add(
        cam->get_projective_equivalent(pose),
        cam->get_ud_pixel(itObs->second.x));
    }
    // Compute the 3D point
    const Vec3 X = trianObj.compute();
    if (trianObj.minDepth() > 0) // Keep the point only if it have a positive depth
    {
      vec_triangulated[iterCP->first] = X;
      vec_control_points[iterCP->first] = landmark.X;
      vec_triangulation_errors[iterCP->first] = trianObj.error()/(double)trianObj.size();
    }
    else
    {
      std::cout << "Invalid triangulation" << std::endl;
      return;
    }
  }

  if (vec_control_points.size() < 3)
  {
    QMessageBox msgBox;
    msgBox.setText("Insufficient number of triangulated control points.");
    msgBox.exec();
    return;
  }

  // compute the similarity
  {
    // data conversion to appropriate container
    Mat x1(3, vec_control_points.size()),
        x2(3, vec_control_points.size());
    for (int i=0; i < vec_control_points.size(); ++i)
    {
      x1.col(i) = vec_triangulated[i];
      x2.col(i) = vec_control_points[i];
    }

    std::cout
      << "Control points observation triangulations:\n"
      << x1 << std::endl << std::endl
      << "Control points coords:\n"
      << x2 << std::endl << std::endl;

    Vec3 t;
    Mat3 R;
    double S;
    if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
    {
      openMVG::geometry::Refine_RTS(x1,x2,&S,&t,&R);
      std::cout << "Found transform:\n"
        << " scale: " << S << "\n"
        << " rotation:\n" << R << "\n"
        << " translation: "<< t.transpose() << std::endl;


      //--
      // Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
      //--

      const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t/S), S);
      openMVG::sfm::ApplySimilarity(sim, m_doc._sfm_data);

      // Display some statistics:
      std::stringstream os;
      for (Landmarks::const_iterator iterL = m_doc._sfm_data.control_points.begin();
        iterL != m_doc._sfm_data.control_points.end(); ++iterL)
      {
        const IndexT CPIndex = iterL->first;
        // If the control point has not been used, continue...
        if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
          continue;

        os
          << "CP index: " << CPIndex << "\n"
          << "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
          << "CP registration error: "
          << (sim(vec_triangulated[CPIndex]) - vec_control_points[CPIndex]).norm() << " user unit(s)"<< "\n\n";
      }
      std::cout << os.str();

      QMessageBox msgBox;
      msgBox.setText(QString::fromStdString(string_pattern_replace(os.str(), "\n", "<br>")));
      msgBox.exec();
    }
    else
    {
      QMessageBox msgBox;
      msgBox.setText("Registration failed. Please check your Control Points coordinates.");
      msgBox.exec();
    }
  }

  //---
  // Bundle adjustment with GCP
  //---
  {
    using namespace openMVG::sfm;
    Bundle_Adjustment_Ceres::BA_Ceres_options options;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    Control_Point_Parameter control_point_opt(20.0, true);
    if (!bundle_adjustment_obj.Adjust(m_doc._sfm_data,
        Optimize_Options
        (
          cameras::Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL, // Adjust structure
          control_point_opt // Use GCP and weight more their observation residuals
          )
        )
      )
    {
      QMessageBox msgBox;
      msgBox.setText("BA with GCP failed.");
      msgBox.exec();
    }
  }
}

MainWindow::MainWindow
(
  QWidget * parent
): QMainWindow()
{
  createPanel();
  createActions();
  createMenus();
  createConnections();

  setWindowTitle(tr("Control_point_editor"));
  QMainWindow::statusBar()->showMessage("Welcome in Control_point_editor GUI.");
  resize(640, 480);
}

void MainWindow::createPanel()
{
  QSplitter *splitter = new QSplitter;
  //-- Create left panel
  m_tabWidget = new QTabWidget;
  //-- Create right panel
  m_widget = new control_point_GUI::GraphicsMainWindow(m_doc, this);
  splitter->addWidget(m_tabWidget);
  splitter->addWidget(m_widget);
  splitter->setStretchFactor(0, 0);
  splitter->setStretchFactor(1, 1);
  setCentralWidget(splitter);

  //-- Add tab inside the m_tabWidget
  m_tab_1 = new QWidget;
  m_tab_1->setObjectName(QString::fromUtf8("m_tab_1"));
  m_tabWidget->addTab(m_tab_1, QString());
  m_tabWidget->setTabText(m_tabWidget->indexOf(m_tab_1), "ImageList");

  //-- Configure tab widgets
  m_treeView_Images = new QTreeView(m_tab_1);
  m_treeView_Images->setRootIsDecorated(false);
  m_treeView_Images->setEditTriggers(QAbstractItemView::NoEditTriggers);
  m_treeView_Images->setObjectName(QString::fromUtf8("m_treeView_Images"));
  m_treeView_Images->setSortingEnabled(true);

  QGridLayout * gridLayout1 = new QGridLayout(m_tab_1);
  gridLayout1->addWidget(m_treeView_Images, 0, 0, 1, 1);
}

void MainWindow::createMenus()
{
  m_menuFile = new QMenu(tr("&File"),this);
  m_menuFile->setObjectName(QString::fromUtf8("m_menuFile"));
  menuBar()->addMenu(m_menuFile);
  QMenu * fileMenu = new QMenu(tr("Open Project"), this);
  m_menuFile->addAction(m_open_action);

  QMenu * m_menuSave = new QMenu(tr("&Save Project"),this);
  m_menuFile->setObjectName(QString::fromUtf8("m_menuSave"));
  m_menuFile->addAction(m_save_action);

  QMenu * m_menuEditCp= new QMenu(tr("Control Point &Edition"),this);
  m_menuEditCp->setObjectName(QString::fromUtf8("m_menuCPEdition"));
  m_menuFile->addAction(m_edit_cp_action);

  QMenu * m_menuRegister = new QMenu(tr("Control Point &Registration"),this);
  m_menuRegister->setObjectName(QString::fromUtf8("m_menuRegister"));
  m_menuFile->addAction(m_register_action);

  QMenu * m_menuDelete = new QMenu(tr("&Delete All Control Point"),this);
  m_menuDelete->setObjectName(QString::fromUtf8("m_menuDelete"));
  m_menuFile->addAction(m_delete_control_point_action);

  m_menuHelp = new QMenu(tr("&Help"),this);
  m_menuHelp->setObjectName(QString::fromUtf8("m_menuHelp"));
  menuBar()->addMenu(m_menuHelp);
  QMenu * helpMenu = new QMenu(tr("Help"), this);
  m_menuHelp->addAction(m_help_action);
}

void MainWindow::createActions()
{
  m_open_action = new QAction(tr("&Open Project..."), this);
  m_open_action->setShortcut(tr("Ctrl+O"));
  connect(m_open_action, SIGNAL(triggered()),
    this, SLOT(openProject()));

  m_save_action = new QAction(tr("&Save Project..."), this);
  m_save_action->setShortcut(tr("Ctrl+S"));
  connect(m_save_action, SIGNAL(triggered()),
    this, SLOT(saveProject()));

  m_register_action = new QAction(tr("Control Point &Registration..."), this);
  m_register_action->setShortcut(tr("Ctrl+R"));
  connect(m_register_action, SIGNAL(triggered()),
    this, SLOT(registerProject()));

  m_edit_cp_action = new QAction(tr("Control Point &Edition..."), this);
  m_edit_cp_action->setShortcut(tr("Ctrl+E"));
  connect(m_edit_cp_action, SIGNAL(triggered()),
    this, SLOT(editControlPoints()));

  m_help_action = new QAction(tr("&Help"), this);
  m_help_action->setShortcut(tr("Ctrl+H"));
  connect(m_help_action, SIGNAL(triggered()),
    this, SLOT(help()));

  m_delete_control_point_action = new QAction(tr("&Delete All Control Point..."), this);
  m_delete_control_point_action->setShortcut(tr("Ctrl+D"));
  connect(m_delete_control_point_action, SIGNAL(triggered()),
    this, SLOT(removeAllControlPoints()));
}

void MainWindow::createConnections()
{
  connect (m_treeView_Images
    ,SIGNAL(doubleClicked(const QModelIndex &))
    ,this
    ,SLOT(doubleClickImageList())
    );
}

