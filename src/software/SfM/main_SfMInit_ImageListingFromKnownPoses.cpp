// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

bool read_Strecha_Camera(const std::string & camName, PinholeCamera & cam, double & width, double & height)
{
  std::ifstream ifs;
  ifs.open( camName.c_str(), std::ifstream::in);
  if (!ifs.is_open()) {
    std::cerr << "Error: failed to open file '" << camName << "' for reading" << std::endl;
    return false;
  }
  std::vector<double> val;
  while (ifs.good() && !ifs.eof())
  {
    double valT;
    ifs >> valT;
    if (!ifs.fail())
      val.push_back(valT);
  }

  if (val.size() == 3*3 +3 +3*3 +3 + 3 || val.size() == 26) //Strecha cam
  {
    Mat3 K, R;
    K << val[0], val[1], val[2],
      val[3], val[4], val[5],
      val[6], val[7], val[8];
    R << val[12], val[13], val[14],
      val[15], val[16], val[17],
      val[18], val[19], val[20];

    Vec3 C (val[21], val[22], val[23]);
    // Strecha model is P = K[R^T|-R^T t];
    // My model is P = K[R|t], t = - RC
    const Vec3 t (-R.transpose() * C);
    R.transposeInPlace();
    cam = cameras::PinholeCamera(K, R, t);

    //image width and height
    width = val[24];
    height = val[25];
  }
  else
  {
    return false;
  }
  return true;
}


//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
  std::string suffix = ".camera";

  CmdLine cmd;

  std::string sImageDir,
    sGroundTruthDir,
    sOutputDir = "";


  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('g', sGroundTruthDir, "groundTruthDirectory") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );


  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-g|--groundTruthDirectory]\n"
      << "[-o|--outputDirectory]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--groundTruthDirectory" << sGroundTruthDir << std::endl
            << "--outputDirectory " << sOutputDir << std::endl;


  if ( !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sGroundTruthDir ) )
  {
    std::cerr << "\nThe input ground truth directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    std::cerr << "\nInvalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      std::cerr << "\nCannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }


  std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  std::sort(vec_image.begin(), vec_image.end());


  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Poses & poses = sfm_data.poses;
  Intrinsics & intrinsics = sfm_data.intrinsics;
  

  C_Progress_display my_progress_bar( vec_image.size(),
      std::cout, "\n- Image listing From Known Poses- \n" );
  std::ostringstream error_report_stream;
  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    ++iter_image, ++my_progress_bar )
  {
   
    const std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream
          << sImFilenamePart << ": Unkown image file format." << "\n";
      continue; // image cannot be opened
    }


    if (sImFilenamePart.find("mask.png") != std::string::npos
       || sImFilenamePart.find("_mask.png") != std::string::npos)
    {
      error_report_stream
          << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }


    // Test if this image can be read
    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

 
    //read ground truth from file
    // Test if ground truth file exists
    std::string sGTFilePath = stlplus::create_filespec(sGroundTruthDir,sImFilenamePart+suffix);
    if (!stlplus::is_file( sGTFilePath ) )
    {
      error_report_stream << std::endl << "There is not valid GT data to read from " << sGTFilePath << std::endl;
      continue;
    }
    else
    {
      // Load GT
      cameras::PinholeCamera cam;
      double width = 0,height = 0;
      if ( read_Strecha_Camera(sGTFilePath, cam, width, height) )
      {
        Mat3 K = cam._K;
        double focal = (K(0,0)+K(1,1))/2.0;//it will be better K(0,0)==K(1,1)
        double pxx   = K(0,2);
        double pyy   = K(1,2);

        std::shared_ptr<View> view                   = std::make_shared<View>(*iter_image, views.size(), views.size(), views.size(), width, height);
        Pose3 pose                                   = Pose3(cam._R,cam._C);
        std::shared_ptr<Pinhole_Intrinsic> intrinsic = std::make_shared<Pinhole_Intrinsic>(width,height,focal,pxx,pyy);

        // Add the view to the sfm_container
        views[view->id_view] = view;
        // Add the pose to the sfm_container
        poses[view->id_pose] = pose;
        // Add the intrinsic to the sfm_container
        intrinsics[view->id_intrinsic] = intrinsic;

      }
      else
      {
        continue;
      }
    }
    
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    std::cerr
      << "\nWarning & Error messages:" << std::endl
      << error_report_stream.str() << std::endl;
  }

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }

  std::cout << std::endl
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Pose(s) listed in sfm_data: " << sfm_data.GetPoses().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;
  return EXIT_SUCCESS;
}
