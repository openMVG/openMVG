// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
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

using namespace std;
using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;


struct ImagePos
{
    std::string image_name;
    std::pair<bool, Vec3> pos;
    std::pair<bool, Vec3> pos_weight;
    
    ImagePos()
    {
        image_name = "";
        pos.first=false;
        pos.second = openMVG::Vec3::Zero();
        pos_weight.first=false;
        pos_weight.second = openMVG::Vec3::Zero();
    }
};

std::vector<ImagePos> readPosFile(const std::string& filename)
{
    std::vector<ImagePos> poss(0);

    std::ifstream pos_file(filename.c_str(), std::ifstream::in);
    if (!pos_file)
    {
      std::cerr << "Error: Failed to open file '" << filename << "' for reading" << std::endl;
      return poss;
    }

    std::string line;
    while (pos_file)
    {
        getline(pos_file,line);
        std::cout << line << std::endl;
        std::vector<std::string> vec_str;
        stl::split(line, ' ', vec_str);
        int vec_size = vec_str.size();
        std::stringstream ss;
        ss.str(line);
        std::string name;
        double x, y, z, x_w, y_w, z_w;
        bool has_pos, has_pos_w;
        x = y = z = x_w = y_w = z_w = 0.0;
        has_pos = has_pos_w = false;
        std::cout << "elements size : "<<vec_size << std::endl;
        if(vec_size == 1)
        {
            ss >> name;
        }
        else if(vec_size == 4) 
        {
            if(ss >> name >> x >> y >> z)
            {
                has_pos = true;
            }
        }
        else if(vec_size == 7) 
        {
            if(ss >> name >> x >> y >> z >> x_w >> y_w >> z_w)
            {
                has_pos = true;
                has_pos_w = true;
            }
        }
        else
        {
            continue;
        }
        
        
        

        ImagePos img_pos;
        img_pos.image_name = name;
        img_pos.pos.first = has_pos;
        img_pos.pos.second = openMVG::Vec3(x, y, z);
        img_pos.pos_weight.first = has_pos_w;
        img_pos.pos_weight.second = openMVG::Vec3(x_w, y_w, z_w);
        poss.emplace_back(img_pos);
    }
    
    pos_file.close();

    return poss;
}


/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ';', vec_str);
  if (vec_str.size() != 9)  {
    std::cerr << "\n Missing ';' character" << std::endl;
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character" << std::endl;
      return false;
    }
    if (i==0) focal = readvalue;
    if (i==2) ppx = readvalue;
    if (i==5) ppy = readvalue;
  }
  return true;
}



//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{

  CmdLine cmd;

  std::string //sImageDir,
    sfileDatabase = "",
    sOutputDir = "",
    sPosFile = "",
    sKmatrix;

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  int i_GPS_XYZ_method = 0;

  double focal_pixels = -1.0;

  cmd.add( make_option('i', sPosFile, "poseFile") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', focal_pixels, "focal") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );
  
  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--poseFile]\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA) << ": Pinhole\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_RADIAL1) << ": Pinhole radial 1\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_RADIAL3) << ": Pinhole radial 3 (default)\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_BROWN) << ": Pinhole brown 2\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_FISHEYE) << ": Pinhole with a simple Fish-eye distortion\n"
      << "\t" << static_cast<int>(CAMERA_SPHERICAL) << ": Spherical camera\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> (default) view can share some camera intrinsic parameters\n"
      << "\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--poseFile " << sPosFile << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << focal_pixels << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl;

  // Expected properties for each image
  double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if (sPosFile.empty())
  {
    std::cerr << "\nInvalid input pose file" << std::endl;
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

  if (sKmatrix.size() > 0 &&
    !checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy) )
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
    return EXIT_FAILURE;
  }

  if (sKmatrix.size() > 0 && focal_pixels != -1.0)
  {
    std::cerr << "\nCannot combine -f and -k options" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      std::cerr
       << "\nInvalid input database: " << sfileDatabase
       << ", please specify a valid file." << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::vector<ImagePos> vec_image_poss = readPosFile(sPosFile);

  if(vec_image_poss.size() == 0)
  {
    std::cerr << "\nImage count is 0, please specify a valid pose file." 
        << std::endl;
    return EXIT_FAILURE;
  }

  //   std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  //   std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  //sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  C_Progress_display my_progress_bar( vec_image_poss.size(),
      std::cout, "\n- Image listing -\n" );
  std::ostringstream error_report_stream;
  for ( std::vector<ImagePos>::const_iterator iter_image = vec_image_poss.begin();
    iter_image != vec_image_poss.end();
    ++iter_image, ++my_progress_bar )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    width = height = ppx = ppy = focal = -1.0;

    const std::string sImageFilename = iter_image->image_name;
    const std::pair<bool, Vec3> gps_info = iter_image->pos;
    const std::pair<bool, Vec3> prior_w_info = iter_image->pos_weight;
    
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported:
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

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

    width = imgHeader.width;
    height = imgHeader.height;
    ppx = width / 2.0;
    ppy = height / 2.0;


    // Consider the case where the focal is provided manually
    if (sKmatrix.size() > 0) // Known user calibration K matrix
    {
      if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
        focal = -1.0;
    }
    else // User provided focal length value
      if (focal_pixels != -1 )
        focal = focal_pixels;

    // If not manually provided or wrongly provided
    if (focal == -1)
    {
      std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
      exifReader->open( sImageFilename );

      const bool bHaveValidExifMetadata =
        exifReader->doesHaveExifInfo()
        && !exifReader->getModel().empty();

      if (bHaveValidExifMetadata) // If image contains meta data
      {
        const std::string sCamModel = exifReader->getModel();

        // Handle case where focal length is equal to 0
        if (exifReader->getFocal() == 0.0f)
        {
          error_report_stream
            << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
          focal = -1.0;
        }
        else
        // Create the image entry in the list file
        {
          Datasheet datasheet;
          if ( getInfo( sCamModel, vec_database, datasheet ))
          {
            // The camera model was found in the database so we can compute it's approximated focal length
            const double ccdw = datasheet.sensorSize_;
            focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
          }
          else
          {
            error_report_stream
              << stlplus::basename_part(sImageFilename)
              << "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
              << "Please consider add your camera model and sensor width in the database." << "\n";
          }
        }
      }
    }
    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic;

    if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
    {
      // Create the desired camera type
      switch (e_User_camera_model)
      {
        case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>
            (width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_RADIAL1:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_RADIAL3:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_BROWN:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_FISHEYE:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case CAMERA_SPHERICAL:
           intrinsic = std::make_shared<Intrinsic_Spherical>
             (width, height);
        break;
        default:
          std::cerr << "Error: unknown camera model: " << (int) e_User_camera_model << std::endl;
          return EXIT_FAILURE;
      }
    }

    // Build the view corresponding to the image
    if (gps_info.first)
    {
      ViewPriors v(sImageFilename, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (intrinsic == nullptr)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      v.b_use_pose_center_ = true;
      v.pose_center_ = gps_info.second;
      // prior weights
      if (prior_w_info.first == true)
      {
        v.center_weight_ = prior_w_info.second;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<ViewPriors>(v);
    }
    else
    {
      View v(sImageFilename, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (intrinsic == nullptr)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<View>(v);
    }
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    std::cerr
      << "\nWarning & Error messages:" << std::endl
      << error_report_stream.str() << std::endl;
  }

  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    GroupSharedIntrinsics(sfm_data);
  }



  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    return EXIT_FAILURE;
  }

  std::cout << std::endl
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image_poss.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;

   return EXIT_SUCCESS;
}

