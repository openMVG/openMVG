// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013, 2014 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/colorHarmonize/colorHarmonizeEngineGlobal.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/system/timer.hpp"

#include <cstdlib>
#include <memory>

using namespace openMVG;

int main( int argc, char **argv )
{
  using namespace std;
  std::cout << "Global Color Harmonization" << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir, sMatchesFile;
  std::string sOutDir = "";
  int selectionMethod = -1;
  int imgRef = -1;

  cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  cmd.add( make_option( 'm', sMatchesFile, "matchesFile" ) );
  cmd.add( make_option( 'o', sOutDir, "outdir" ) );
  cmd.add( make_option( 's', selectionMethod, "selectionMethod" ) );
  cmd.add( make_option( 'r', imgRef, "referenceImage" ) );

  try
  {
    if (argc == 1 ) throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch ( const std::string& s )
  {
    OPENMVG_LOG_INFO
      << "Usage: " << argv[ 0 ] << '\n'
      << "[-i|--input_file] path to a SfM_Data scene\n"
      << "[-m|--sMatchesFile path] i.e path/matches.(h/f/e).txt\n"
      << "[-o|--outdir path]\n"
      << "\n[Optional]\n"
      << "[-s|--selectionMethod int]\n"
      << "[-r|--referenceImage int]";

    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  if ( sSfM_Data_Filename.empty() )
  {
    OPENMVG_LOG_ERROR << "\nIt is an invalid file input";
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutDir ) )
    stlplus::folder_create( sOutDir );

  //---------------------------------------
  // Harmonization process
  //---------------------------------------

  openMVG::system::Timer timer;

  sMatchesDir = stlplus::folder_part(sMatchesFile);
  std::unique_ptr<ColorHarmonizationEngineGlobal> m_colorHarmonizeEngine(
    new ColorHarmonizationEngineGlobal(sSfM_Data_Filename,
    sMatchesDir,
    sMatchesFile,
    sOutDir,
    selectionMethod,
    imgRef));

  if ( m_colorHarmonizeEngine->Process() )
  {
    OPENMVG_LOG_INFO << "ColorHarmonization took (s): " << timer.elapsed();

    return EXIT_SUCCESS;
  }

  OPENMVG_LOG_ERROR << "Something goes wrong in the process";
  return EXIT_FAILURE;
}
