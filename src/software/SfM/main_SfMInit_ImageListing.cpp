// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include "openMVG_Samples/sensorWidthDatabase/ParseDatabase.hpp"

#include "openMVG/image/image.hpp"
#include "openMVG/stl/split.hpp"

#include "openMVG/sfm/sfm.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::image;
using namespace openMVG::sfm;

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ";", vec_str);
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

  std::string sImageDir,
    sfileDatabase = "",
    sOutputDir = "",
    sKmatrix;

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  double focalPixPermm = -1.0;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', focalPixPermm, "focal") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t 1: Pinhole\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> (default) view can share some camera intrinsic parameters\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << focalPixPermm << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl;

  // Expected properties for each image
  double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
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

  if (sKmatrix.size() > 0 && focalPixPermm != -1.0)
  {
    std::cerr << "\nCannot combine -f and -k options" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      std::cerr << "\nInvalid input database" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    iter_image++ )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    width = height = ppx = ppy = focal = -1.0;

    const std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
      continue; // image cannot be opened

    Image<unsigned char> image;
    if (openMVG::image::ReadImage( sImageFilename.c_str(), &image))  {
      width = image.Width();
      height = image.Height();
      ppx = width / 2.0;
      ppy = height / 2.0;
    }
    else
      continue; // image cannot be read

    std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif());
    exifReader->open( sImageFilename );

    // Consider the case where focal is provided manually
    if ( !exifReader->doesHaveExifInfo() || focalPixPermm != -1)
    {
      if (sKmatrix.size() > 0) // Known user calibration K matrix
      {
        if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
          focal = -1.0;
      }
      else // User provided focal lenght value
        if (focalPixPermm != -1 )
        focal = focalPixPermm;
    }
    else // If image contains meta data
    {
      const std::string sCamName = exifReader->getBrand();
      const std::string sCamModel = exifReader->getModel();

      // Handle case where focal length is equal to 0
      if (exifReader->getFocal() == 0.0f)
      {
        std::cerr << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << std::endl;
        continue;
      }

      // Create the image entry in the list file
      {
        Datasheet datasheet;
        if ( getInfo( sCamName, sCamModel, vec_database, datasheet ))
        {
          // The camera model was found in the database so we can compute it's approximated focal length
          const double ccdw = datasheet._sensorSize;
          focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
        }
        else
        {
          std::cout << stlplus::basename_part(sImageFilename) << ": Camera \""
            << sCamName << "\" model \"" << sCamModel << "\" doesn't exist in the database" << std::endl
            << "Please consider add your camera model and sensor width in the database." << std::endl;
        }
      }
    }

    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic (NULL);

    if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
    {
      // Create the desired camera type
      switch(e_User_camera_model)
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
        default:
          std::cout << "Unknown camera model: " << (int) e_User_camera_model << std::endl;
      }
    }

    // Build the view corresponding to the image
    View v(*iter_image, views.size(), views.size(), views.size(), width, height);

    // Add intrinsic related to the image (if any)
    if (intrinsic == NULL)
    {
      //Since the view have invalid intrinsic data
      // (export the view, with an invalid intrinsic field value)
      v.id_intrinsic = UndefinedIndexT;
    }
    else
    {
      // Add the intrinsic to the sfm_container
      intrinsics[v.id_intrinsic] = intrinsic;
    }

    // Add the view to the sfm_container
    views[v.id_view] = std::make_shared<View>(v);
  }

  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    // Group camera model that share common optics camera properties
    // They must share (camera model, image size, & camera parameters)
    // Grouping is simplified by using a hash function over the camera intrinsic.

    // Build hash & build a set of the hash in order to maintain unique Ids
    std::set<size_t> hash_index;
    std::vector<size_t> hash_value;

    for (Intrinsics::const_iterator iterIntrinsic = intrinsics.begin();
      iterIntrinsic != intrinsics.end();
      ++iterIntrinsic)
    {
      const IntrinsicBase * intrinsicData = iterIntrinsic->second.get();
      const size_t hashVal = intrinsicData->hashValue();
      hash_index.insert(hashVal);
      hash_value.push_back(hashVal);
    }

    // From hash_value(s) compute the new index (old to new indexing)
    Hash_Map<IndexT, IndexT> old_new_reindex;
    size_t i = 0;
    for (Intrinsics::const_iterator iterIntrinsic = intrinsics.begin();
      iterIntrinsic != intrinsics.end();
      ++iterIntrinsic, ++i)
    {
      old_new_reindex[iterIntrinsic->first] = std::distance(hash_index.begin(), hash_index.find(hash_value[i]));
      //std::cout << hash_value[i] // hash
      //  << "\t" << iterIntrinsic->first // old reference index
      // << "\t" << old_new_reindex[iterIntrinsic->first] << std::endl; // new index
    }
    //--> Copy & modify Ids & replace
    //     - for the Intrinsic params
    //     - for the View
    Intrinsics intrinsic_updated;
    for (Intrinsics::const_iterator iterIntrinsic = intrinsics.begin();
      iterIntrinsic != intrinsics.end();
      ++iterIntrinsic)
    {
      intrinsic_updated[old_new_reindex[iterIntrinsic->first]] = intrinsics[iterIntrinsic->first];
    }
    intrinsics.swap(intrinsic_updated); // swap camera intrinsics
    // Update intrinsic ids
    for (Views::iterator iterView = views.begin();
      iterView != views.end();
      ++iterView)
    {
      View * v = iterView->second.get();
      v->id_intrinsic = old_new_reindex[v->id_intrinsic];
    }
  }

  // Store SfM_Data views & intrinsic data
  if (Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(VIEWS|INTRINSICS))
    )
  return EXIT_SUCCESS;
    else
  return EXIT_FAILURE;
}
