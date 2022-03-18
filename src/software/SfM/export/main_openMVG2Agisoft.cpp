// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Etienne Danvoye

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <iomanip>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir] path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--sfmdata " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl;

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create( sOutDir ))
    {
      std::cerr << "Cannot create the destination folder." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Read the SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream outfile(stlplus::create_filespec(sOutDir, "cameras", "xml").c_str());

  if (!outfile)
  {
    std::cerr << "Cannot open the destination file." << std::endl;
    return EXIT_FAILURE;
  }

  outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  outfile << "<document version=\"1.3.0\">\n";
  outfile << "<chunk>\n";

  outfile << "<sensors>\n";
  for (const auto& intrinsic : sfm_data.GetIntrinsics())
  {
    if (isPinhole(intrinsic.second->getType()))
    {
      const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(intrinsic.second.get());

      outfile << std::setprecision(16) <<
        "<sensor id=\"" << intrinsic.first << "\" label=\"sensor_" << intrinsic.first << "\" type=\"frame\">\n" <<
        "<resolution width=\"" << cam->w() << "\" height=\"" << cam->h() << "\"/>\n" <<
        "<property name=\"fixed\" value=\"false\"/>\n" <<
        "<calibration type=\"frame\" class=\"adjusted\">\n" <<
        "<resolution width=\"" << cam->w() << "\" height=\"" << cam->h() << "\"/>\n";

      outfile <<
        "<fx>" << cam->focal() << "</fx>\n" <<
        "<fy>" << cam->focal() << "</fy>\n";

      outfile <<
        "<cx>" << cam->principal_point()[0] << "</cx>\n" <<
        "<cy>" << cam->principal_point()[1] << "</cy>\n";

      if (cam->have_disto())
      {
        auto params = cam->getParams();
        switch (cam->getType())
        {
        case PINHOLE_CAMERA_RADIAL1:
          outfile <<
            "<k1>" << params[3] << "</k1>\n";
          break;
        case PINHOLE_CAMERA_RADIAL3:
          outfile <<
            "<k1>" << params[3] << "</k1>\n"
            "<k2>" << params[4] << "</k2>\n"
            "<k3>" << params[5] << "</k3>\n";
          break;
        case PINHOLE_CAMERA_BROWN:
          outfile <<
            "<k1>" << params[3] << "</k1>\n"
            "<k2>" << params[4] << "</k2>\n"
            "<p1>" << params[6] << "</p1>\n"
            "<p2>" << params[7] << "</p2>\n"
            "<k3>" << params[5] << "</k3>\n";
          break;
        default:
          std::cerr << "Unsupported camera model for export." << std::endl;
          return EXIT_FAILURE;
        }
      }

      outfile << "</calibration>\n";
      outfile << "</sensor>\n";
    }
  }
  outfile << "</sensors>\n";

  outfile << "<cameras>\n";
  for (const auto& view : sfm_data.GetViews())
  {
    if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
    {
      const openMVG::geometry::Pose3 poseMVG(sfm_data.GetPoseOrDie(view.second.get()));
      auto mat34 = poseMVG.inverse().asMatrix();

      outfile
        << "<camera id=\"" << view.first
        << "\" label=\"" << stlplus::basename_part(view.second->s_Img_path)
        << "\" sensor_id=\"" << view.second->id_intrinsic
        << "\" enabled=\"1\">\n"
        << "<transform>" << mat34 << " 0.0 0.0 0.0 1.0</transform>\n"
        << "</camera>\n";
    }
  }
  outfile << "</cameras>\n";

  outfile << "<region>\n"
    "<center>0 0 0 </center>\n"
    "<size>100 100 100 </size>\n"
    "<R>1 0 0 0 1 0 0 0 1 </R>\n"
    "</region>\n"
    "<settings>\n"
    "<property name=\"accuracy_tiepoints\" value=\"1\"/>\n"
    "<property name=\"accuracy_cameras\" value=\"10\" />\n"
    "<property name=\"accuracy_cameras_ypr\" value=\"2\" />\n"
    "<property name=\"accuracy_markers\" value=\"0.005\" />\n"
    "<property name=\"accuracy_scalebars\" value=\"0.001\" />\n"
    "<property name=\"accuracy_projections\" value=\"0.1\" />\n"
    "</settings>\n";

  outfile <<
    "</chunk>\n"
    "</document>\n";

  if (!outfile)
  {
    std::cerr << "Something went wrong during the export." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
