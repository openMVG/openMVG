
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
  std::cout << "Global Structure from Motion :"
            << " open Source implementation of: \n"
            << "\"Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion.\""
            << " ICCV 2013 paper. Pierre Moulon, Pascal Monasse and Renaud Marlet." << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sImaDirectory;
  std::string sMatchesDir;
  std::string sOutDir = "";

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << ' '
    << "[-i|--imadir path] "
    << "[-m|--matchdir path] "
    << "[-o|--outdir path] "
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
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
