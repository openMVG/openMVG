
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "software/globalSfM/SfMGlobalEngine.hpp"
using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Global Structure from Motion:\n"
    << "-----------------------------------------------------------\n"
    << "Open Source implementation of the paper:\n"
    << "\"Global Fusion of Relative Motions for "
    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
    << " ICCV 2013." << std::endl
    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sImaDirectory;
  std::string sMatchesDir;
  std::string sOutDir = "";
  bool bColoredPointCloud = false;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('c', bColoredPointCloud, "coloredPointCloud") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path]\n"
    << "[-m|--matchdir path]\n"
    << "[-o|--outdir path]\n"
    << "[-c|--coloredPointCloud 0(default) or 1]\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Incremental reconstruction process
  //---------------------------------------

  clock_t timeStart = clock();
  GlobalReconstructionEngine to3DEngine(sImaDirectory,
                                            sMatchesDir,
                                            sOutDir,
                                            true);

  if (to3DEngine.Process())
  {
    clock_t timeEnd = clock();
    std::cout << std::endl << " Total Ac-Global-Sfm took : " << (timeEnd - timeStart) / CLOCKS_PER_SEC << std::endl;
    
    //-- Export computed data to disk
    to3DEngine.ExportToOpenMVGFormat(bColoredPointCloud);

    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
