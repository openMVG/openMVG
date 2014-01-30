
// Copyright (c) 2013, 2014 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <memory>

#include "software/colorHarmonize/colorHarmonizeEngineGlobal.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;

int main( int argc, char **argv )
{
  using namespace std;
  std::cout << "Global Color Harmonization" << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sImaDirectory;
  std::string sMatchesDir, sMatchesFile;
  std::string sOutDir = "";
  int selectionMethod = -1;
  int imgRef = -1;

  cmd.add( make_option( 'i', sImaDirectory, "imagesDirectory" ) );
  cmd.add( make_option( 'm', sMatchesFile, "matchesFile" ) );
  cmd.add( make_option( 'o', sOutDir, "outdir" ) );
  cmd.add( make_option( 's', selectionMethod, "selectionMethod" ) );
  cmd.add( make_option( 'r', imgRef, "referenceImage" ) );

  try
  {
    if( argc == 1 ) throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch( const std::string& s )
  {
    std::cerr << "Usage: " << argv[ 0 ] << ' '
    << "[-i|--imagesDirectory path] "
    << "[-m|--sMatchesFile path] "
    << "[-o|--outdir path] "
    << "[-s|--selectionMethod int] "
    << "[-r|--referenceImage int]"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if ( sImaDirectory.empty() )
  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutDir ) )
    stlplus::folder_create( sOutDir );

  //---------------------------------------
  // Harmonization process
  //---------------------------------------

  clock_t timeStart = clock();
  
  sMatchesDir = stlplus::folder_part(sMatchesFile);
  std::auto_ptr<ReconstructionEngine> m_colorHarmonizeEngine;
  {
    m_colorHarmonizeEngine = std::auto_ptr<ReconstructionEngine>(
      new ColorHarmonizationEngineGlobal(sImaDirectory,
      sMatchesDir,
      sMatchesFile,
      sOutDir,
      selectionMethod,
      imgRef));
  }

  if ( m_colorHarmonizeEngine->Process() )
  {
    clock_t timeEnd = clock();
    std::cout << std::endl
      << " ColorHarmonization took : "
      << (timeEnd - timeStart) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "\n Something goes wrong in the process" << std::endl;
  }
  return EXIT_FAILURE;
}
