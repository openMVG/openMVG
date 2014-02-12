// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifdef USE_EXIV2
#include "openMVG/exif_IO/exif_IO_Exiv2.hpp"
#else
#include "openMVG/exif_IO/exif_IO_openExif.hpp"
#endif

#include "openMVG_Samples/sensorWidthDatabase/ParseDatabase.hpp"

#include "openMVG/image/image.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "JsonBox.h"

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir;
  std::string sfileDatabase = "";
  std::string sOutputDir = "";
  double focalPixPermm = -1.0;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', focalPixPermm, "focal") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--imageDirectory] "
      << "[-d|--sensorWidthDatabase] "
      << "[-o|--outputDirectory] "
      << "[-f|--focal] "
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << focalPixPermm << std::endl;

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    stlplus::folder_create( sOutputDir );
  }

  std::vector<std::string> vec_image = stlplus::folder_all( sImageDir );

  // Create the new file
  JsonBox::Object imageParams;
  JsonBox::Array images;
  std::ofstream imageParamsFile( stlplus::create_filespec( sOutputDir,
                                                   "imageParams.json" ).c_str() );
  if ( !imageParamsFile )
  {
    std::cerr << "\nCould not create imageParams.json" << std::endl;
    return EXIT_FAILURE;
  }

  std::sort(vec_image.begin(), vec_image.end());
  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    iter_image++ )
  {
    // Read meta data to fill width height and focalPixPermm
    std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );

    size_t width = -1;
    size_t height = -1;

#ifdef USE_EXIV2
    std::auto_ptr<Exif_IO> exifReader (new Exif_IO_Exiv2() );
#else
    std::auto_ptr<Exif_IO> exifReader (new Exif_IO_OpenExif() );
#endif
    exifReader->open( sImageFilename );

    // Consider the case where focal is provided

    std::ostringstream os;
    //If image do not contains meta data
    if ( !exifReader->doesHaveExifInfo() || focalPixPermm != -1)
    {
      Image<unsigned char> image;
      if (openMVG::ReadImage( sImageFilename.c_str(), &image))
      {
        width = image.Width();
        height = image.Height();
      }
      else
      {
        Image<RGBColor> imageRGB;
        if (openMVG::ReadImage( sImageFilename.c_str(), &imageRGB))
        {
          width = imageRGB.Width();
          height = imageRGB.Height();
        }
        else
        {
          continue; // image is not considered, cannot be read
        }
      }
      if ( focalPixPermm == -1)
      {
        JsonBox::Object image;
        image["filename"] = JsonBox::Value(*iter_image);
        image["width"] = JsonBox::Value(int(width));
        image["height"] = JsonBox::Value(int(height));
        images.push_back(image);
      }
      else
      {
        JsonBox::Object image;
        image["filename"] = JsonBox::Value(*iter_image);
        image["width"] = JsonBox::Value(int(width));
        image["height"] = JsonBox::Value(int(height));
        image["focalLength"] = JsonBox::Value(focalPixPermm);
        JsonBox::Object principalPoint;
        principalPoint["x"] = JsonBox::Value(width/2.0);
        principalPoint["y"] = JsonBox::Value(height/2.0);
        image["principalPoint"] = principalPoint;
        images.push_back(image);
      }

    }
    else // If image contains meta data
    {
      double focal = focalPixPermm;
      width = exifReader->getWidth();
      height = exifReader->getHeight();
      std::string sCamName = exifReader->getBrand();
      std::string sCamModel = exifReader->getModel();

      std::vector<Datasheet> vec_database;
      Datasheet datasheet;
      if ( parseDatabase( sfileDatabase, vec_database ) )
      {
        if ( getInfo( sCamName, sCamModel, vec_database, datasheet ) )
        {
          // The camera model was found in the database so we can compute it's approximated focal length
          double ccdw = datasheet._sensorSize;
          focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
          JsonBox::Object image;
          image["filename"] = JsonBox::Value(*iter_image);
          image["width"] = JsonBox::Value(int(width));
          image["height"] = JsonBox::Value(int(height));
          image["focalLength"] = JsonBox::Value(focal);
          JsonBox::Object camera;
          camera["name"] = JsonBox::Value(sCamName);
          camera["model"] = JsonBox::Value(sCamModel);
          image["camera"] = camera;
          images.push_back(image);
        }
        else
        {
          std::cout << "Camera \"" << sCamName << "\" model \"" << sCamModel << "\" doesn't exist in the database" << std::endl;
          JsonBox::Object image;
          image["filename"] = JsonBox::Value(*iter_image);
          image["width"] = JsonBox::Value(int(width));
          image["height"] = JsonBox::Value(int(height));
          JsonBox::Object camera;
          camera["name"] = JsonBox::Value(sCamName);
          camera["model"] = JsonBox::Value(sCamModel);
          image["camera"] = camera;
          images.push_back(image);
        }
      }
      else
      {
        std::cout << "Sensor width database \"" << sfileDatabase << "\" doesn't exist." << std::endl;
        std::cout << "Please consider add your camera model in the database." << std::endl;
        JsonBox::Object image;
        image["filename"] = JsonBox::Value(*iter_image);
        image["width"] = JsonBox::Value(int(width));
        image["height"] = JsonBox::Value(int(height));
        JsonBox::Object camera;
        camera["name"] = JsonBox::Value(sCamName);
        camera["model"] = JsonBox::Value(sCamModel);
        image["camera"] = camera;
        images.push_back(image);
      }
    }
  }
  imageParams["images"] = images;
  imageParamsFile << imageParams;
  imageParamsFile.close();
  return EXIT_SUCCESS;
}
