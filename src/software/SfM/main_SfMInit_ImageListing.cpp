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
#include <experimental/filesystem>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrixs, std::vector<double> & focals,
                                  std::vector<double> & ppxs, std::vector<double> & ppys)
{
  focals.clear();
  ppxs.clear();
  ppys.clear();
  std::vector<std::string> vec_str;
  stl::split(Kmatrixs, ';', vec_str);
  if (vec_str.size() == 0 || vec_str.size() % 9 != 0)  {
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
    size_t component = i % 9;
    if (component==0) focals.push_back(readvalue);
    if (component==2) ppxs.push_back(readvalue);
    if (component==5) ppys.push_back(readvalue);
  }
  return true;
}

bool checkFocalPixelsStringValidity(const std::string &sFocalPixels, std::vector<double> &focals) {
  focals.clear();
  std::vector<std::string> vec_str;
  std::string focal_pixels = sFocalPixels + ";";
  stl::split(focal_pixels, ';', vec_str);
  if (vec_str.size() == 0) {
    std::cerr << "\n Missing ';' character" << std::endl;
    return false;
  }
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue)) {
      std::cerr << "\n Used an invalid not a number character" << std::endl;
      return false;
    }
    focals.push_back(readvalue);
  }

  return true;
}

std::vector<std::string> checkImageDirsStringValidity(const std::string& sImageDirs) {
  std::vector<std::string> imageDirs;
  std::vector<std::string> vec_str;
  stl::split(sImageDirs, ';', vec_str);
  if (vec_str.size() == 0) {
    throw std::runtime_error("\n Missing ';' character");
  }
  // Check that all ImageDir all exist
  for (size_t i = 0; i < vec_str.size(); ++i) {
    if (!stlplus::folder_exists(vec_str[i])) {
      throw std::runtime_error(vec_str[i] + " does not exist.");
    }
    imageDirs.push_back(std::move(vec_str[i]));
  }

  return imageDirs;
}

inline bool image_filename_less(const std::pair<std::string, size_t> &p1, const std::pair<std::string, size_t> &p2) {
  return (p1.first.compare(p2.first) < 0);
}

std::vector<std::pair<std::string, size_t>> aggregateImageDirs(std::vector<std::string> &imageDirs, std::string &sInputFolder) {
  std::vector<std::pair<std::string, size_t>> sortedImages;

  if (imageDirs.size() == 1) {
    std::vector<std::string> vec_image = stlplus::folder_files(imageDirs[0]);
    std::sort(vec_image.begin(), vec_image.end());
    for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin(); iter_image != vec_image.end(); ++iter_image) {
      if (openMVG::image::GetFormat((*iter_image).c_str()) == openMVG::image::Unknown)
        continue;
      sortedImages.push_back(std::make_pair(*iter_image, 0));
    }
  } else if (imageDirs[0] != sInputFolder){
    for (size_t i=0; i<imageDirs.size(); ++i) {
      std::vector<std::string> vec_image = stlplus::folder_files(imageDirs[i]);
      for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin(); iter_image != vec_image.end(); ++iter_image) {
        if (openMVG::image::GetFormat((*iter_image).c_str()) == openMVG::image::Unknown)
          continue;
        std::string imageFilenameSrc = stlplus::create_filespec(imageDirs[i], *iter_image);
        std::string imageFilenameDest = stlplus::create_filespec(sInputFolder, *iter_image);
        if (!std::experimental::filesystem::exists(imageFilenameDest)) {
          std::experimental::filesystem::remove(imageFilenameDest);
        }
        std::experimental::filesystem::create_symlink(imageFilenameSrc, imageFilenameDest);
        sortedImages.push_back(std::make_pair(*iter_image, i));
      }
    }
    std::sort(sortedImages.begin(), sortedImages.end(), image_filename_less);
  } else {
      std::vector<std::string> vec_image = stlplus::folder_files(imageDirs[0]);
      for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin(); iter_image != vec_image.end(); ++iter_image) {
        if (openMVG::image::GetFormat((*iter_image).c_str()) == openMVG::image::Unknown)
          continue;
        if (!std::experimental::filesystem::is_symlink(*iter_image)) {
          sortedImages.push_back(std::make_pair(*iter_image, 0));
        }
      }
      for (size_t i=1; i<imageDirs.size(); ++i) {
        std::vector<std::string> vec_image = stlplus::folder_files(imageDirs[i]);
        for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin(); iter_image != vec_image.end(); ++iter_image) {
          if (openMVG::image::GetFormat((*iter_image).c_str()) == openMVG::image::Unknown)
            continue;
          std::string imageFilenameSrc = stlplus::create_filespec(imageDirs[i], *iter_image);
          std::string imageFilenameDest = stlplus::create_filespec(sInputFolder, *iter_image);
          if (!std::experimental::filesystem::exists(imageFilenameDest)) {
            std::experimental::filesystem::remove(imageFilenameDest);
          }
          std::experimental::filesystem::create_symlink(imageFilenameSrc, imageFilenameDest);
          sortedImages.push_back(std::make_pair(*iter_image, i));
        }
      }
    std::sort(sortedImages.begin(), sortedImages.end(), image_filename_less);
  }

  return sortedImages;
}

std::pair<bool, Vec3> checkGPS
(
  const std::string & filename,
  const int & GPS_to_XYZ_method = 0
)
{
  std::pair<bool, Vec3> val(false, Vec3::Zero());
  std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
  if (exifReader)
  {
    // Try to parse EXIF metada & check existence of EXIF data
    if ( exifReader->open( filename ) && exifReader->doesHaveExifInfo() )
    {
      // Check existence of GPS coordinates
      double latitude, longitude, altitude;
      if ( exifReader->GPSLatitude( &latitude ) &&
           exifReader->GPSLongitude( &longitude ) &&
           exifReader->GPSAltitude( &altitude ) )
      {
        // Add ECEF or UTM XYZ position to the GPS position array
        val.first = true;
        switch (GPS_to_XYZ_method)
        {
          case 1:
            val.second = lla_to_utm( latitude, longitude, altitude );
            break;
          case 0:
          default:
            val.second = lla_to_ecef( latitude, longitude, altitude );
            break;
        }
      }
    }
  }
  return val;
}


/// Check string of prior weights
std::pair<bool, Vec3> checkPriorWeightsString
(
  const std::string &sWeights
)
{
  std::pair<bool, Vec3> val(true, Vec3::Zero());
  std::vector<std::string> vec_str;
  stl::split(sWeights, ';', vec_str);
  if (vec_str.size() != 3)
  {
    std::cerr << "\n Missing ';' character in prior weights" << std::endl;
    val.first = false;
  }
  // Check that all weight values are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i)
  {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character in local frame origin" << std::endl;
      val.first = false;
    }
    val.second[i] = readvalue;
  }
  return val;
}
//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDirs;
  std::string sCombinationDirectory;
  std::string sfileDatabase;
  std::string sOutputDir;
  std::string sKmatrixs;
  std::string sFocalPixels;

  std::string sPriorWeights;
  std::pair<bool, Vec3> prior_w_info(false, Vec3(1.0,1.0,1.0));

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  int i_GPS_XYZ_method = 0;


  cmd.add( make_option('a', sCombinationDirectory, "combinationDirectory"));
  cmd.add( make_option('i', sImageDirs, "imageDirectories") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', sFocalPixels, "focals") );
  cmd.add( make_option('k', sKmatrixs, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );
  cmd.add( make_switch('P', "use_pose_prior") );
  cmd.add( make_option('W', sPriorWeights, "prior_weights"));
  cmd.add( make_option('m', i_GPS_XYZ_method, "gps_to_xyz_method") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-a|--combinationDirectory]: if you have multiple intrinsics, you must offer the argument to accommodate multiple source image.\n"
      << "[-i|--imageDirectories] sImageDirs: \"input1;input2...\"\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focals] (pixels)\n"
      << "[-k|--intrinsics] Kmatrixs: \"f;0;ppx;0;f;ppy;0;0;1...\"\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t 1: Pinhole\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole brown 2\n"
      << "\t 5: Pinhole with a simple Fish-eye distortion\n"
      << "\t 6: Pinhole radial 1 pba\n"
      << "\t 8: Spherical camera\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> (default) view can share some camera intrinsic parameters\n"
      << "\n"
      << "[-P|--use_pose_prior] Use pose prior if GPS EXIF pose is available"
      << "[-W|--prior_weights] \"x;y;z;\" of weights for each dimension of the prior (default: 1.0)\n"
      << "[-m|--gps_to_xyz_method] XZY Coordinate system:\n"
      << "\t 0: ECEF (default)\n"
      << "\t 1: UTM\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--combinationDirectory " << sCombinationDirectory << std::endl
            << "--imageDirectories " << sImageDirs << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << sFocalPixels << std::endl
            << "--intrinsics " << sKmatrixs << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl;

  // Expected properties for each image
  std::vector<double> focals, ppxs, ppys;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  auto imageDirs = checkImageDirsStringValidity(sImageDirs);
  if (imageDirs.empty()) {
    return EXIT_FAILURE;
  }

  if (imageDirs.size() != 1 && sCombinationDirectory.size() > 0) {
    if ( !stlplus::folder_exists( sCombinationDirectory ) ) {
      if ( !stlplus::folder_create( sCombinationDirectory ) ) {
        std::cerr << "\nCannot create " << sCombinationDirectory << " directory" << std::endl;
        return EXIT_FAILURE;
      }
    } else {
        if (!std::experimental::filesystem::is_empty(sCombinationDirectory)) {
          std::experimental::filesystem::remove_all(sCombinationDirectory);
          if ( !stlplus::folder_create( sCombinationDirectory ) ) {
            std::cerr << "\nCannot create " << sCombinationDirectory << " directory" << std::endl;
            return EXIT_FAILURE;
          }
        }
    }
  } else {
    sCombinationDirectory = imageDirs[0];
  }

  auto sortedImages = aggregateImageDirs(imageDirs, sCombinationDirectory);

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

  if (sFocalPixels.size() > 0 && !checkFocalPixelsStringValidity(sFocalPixels, focals)) {
    std::cerr << "\nInvalid f focals input" << std::endl;
    return EXIT_FAILURE;
  }

  bool is_valid_intrinsic = checkIntrinsicStringValidity(sKmatrixs, focals, ppxs, ppys);
  if (sKmatrixs.size() > 0 &&
    !is_valid_intrinsic)
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
    return EXIT_FAILURE;
  }

  if (sKmatrixs.size() > 0 && sFocalPixels.size() > 0)
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

  // Check if prior weights are given
  if (cmd.used('P') && !sPriorWeights.empty())
  {
    prior_w_info = checkPriorWeightsString(sPriorWeights);
  }
  else if (cmd.used('P'))
  {
    prior_w_info.first = true;
  }

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sCombinationDirectory; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  C_Progress_display my_progress_bar( sortedImages.size(),
      std::cout, "\n- Image listing -\n" );
  std::ostringstream error_report_stream;
  for (size_t image_number = 0; image_number<sortedImages.size(); ++image_number, ++my_progress_bar) {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    double width = -1.0, height = -1.0, ppx = -1.0, ppy = -1.0, focal = -1.0;

    const std::string sImageFilename = stlplus::create_filespec( sCombinationDirectory, sortedImages[image_number].first);
    const std::string sImFilenamePart = stlplus::filename_part( sImageFilename);

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
    if (sKmatrixs.size() > 0) // Known user calibration K matrix
    {
      if (is_valid_intrinsic) {
        focal = focals[sortedImages[image_number].second];
        ppx = ppxs[sortedImages[image_number].second];
        ppy = ppys[sortedImages[image_number].second];
      } else {
        focal = -1.0;
      }
    }
    else // User provided focal length value
      if (sFocalPixels.size() > 0)
        focal = focals[sortedImages[image_number].second];

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
        case PINHOLE_CAMERA_RADIAL1_PBA:
           intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1_PBA>
             (width, height, focal, ppx, ppy, 0.0);
           break;
        default:
          std::cerr << "Error: unknown camera model: " << (int) e_User_camera_model << std::endl;
          return EXIT_FAILURE;
      }
    }

    // Build the view corresponding to the image
    const std::pair<bool, Vec3> gps_info = checkGPS(sImageFilename, i_GPS_XYZ_method);
    if (gps_info.first && cmd.used('P'))
    {
      ViewPriors v(sortedImages[image_number].first, views.size(), views.size(), views.size(), width, height);

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
      View v(sortedImages[image_number].first, views.size(), views.size(), views.size(), width, height);

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
    << "listed #File(s): " << sortedImages.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size() << std::endl;

  return EXIT_SUCCESS;
}
