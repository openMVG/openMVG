// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/geodesy/geodesy.hpp"

#include "software/SfM/SfMPlyHelper.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/cmdLine/cmdLine.h"

#include <iostream>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace std;

int main(int argc, char **argv)
{
  std::string
    sInputDirectory = "",
    sOutputPLYFile = "GPS_POSITION.ply";

  CmdLine cmd;
  cmd.add(make_option('i', sInputDirectory, "input-directory") );
  cmd.add(make_option('o', sOutputPLYFile, "output-file") );

  try
  {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string& s)
  {
    std::cout
    << "Geodesy demo:\n"
    << " Export as PLY points the parsed image EXIF GPS positions,\n"
    << " -[i|input-directory] Directory that will be parsed.\n"
    << "-- OPTIONAL PARAMETERS --\n"
    << " -[o|output-file] Output PLY file.\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Loop: Parse images
  //  | If valid GPS data information are found, convert them to XYZ & add it to an array
  // Export all XYZ parsed data to a PLY file

  // Init the EXIF reader
  std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
  if (!exifReader)
  {
    std::cerr << "Cannot instantiate the EXIF metadata reader." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Vec3> vec_gps_xyz_position;

  size_t valid_exif_count = 0;

  const std::vector<std::string> vec_image = stlplus::folder_files( sInputDirectory );
  for (const std::string & it_string : vec_image)
  {
    const std::string sImageFilename = stlplus::create_filespec( sInputDirectory, it_string );

    // Try to parse EXIF metada
    exifReader->open( sImageFilename );

    // Check existence of EXIF data
    if (!exifReader->doesHaveExifInfo())
      continue;

    ++valid_exif_count;
    // Check existence of GPS coordinates
    double latitude, longitude, altitude;
    if ( exifReader->GPSLatitude( &latitude ) &&
         exifReader->GPSLongitude( &longitude ) &&
         exifReader->GPSAltitude( &altitude ) )
    {
      // Add ECEF XYZ position to the GPS position array
      vec_gps_xyz_position.push_back( lla_to_ecef( latitude, longitude, altitude ) );
    }
  }

  std::cout << std::endl
    << "Report:\n"
    << " #file listed: " << vec_image.size() << "\n"
    << " #valid exif data: " << valid_exif_count << "\n"
    << " #valid exif gps data: " << vec_gps_xyz_position.size() << std::endl;

  if ( vec_gps_xyz_position.empty() )
  {
    std::cerr << "No valid GPS data found for the image list" << std::endl;
    return EXIT_FAILURE;
  }

  if ( plyHelper::exportToPly( vec_gps_xyz_position, sOutputPLYFile ) )
  {
    std::cout << sOutputPLYFile << " -> successfully exported on disk." << std::endl;
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
