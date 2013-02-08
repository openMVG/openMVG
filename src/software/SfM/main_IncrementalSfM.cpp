
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "software/SfM/SfMIncrementalEngine.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Incremental reconstruction" << std::endl
            << " Perform incremental SfM (Initial Pair Essential + Resection)." << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sImaDirectory;
  std::string sMatchesDir;
  std::string sOutDir = "";
  bool bPmvsExport = false;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('p', bPmvsExport, "pmvs") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << ' '
    << "[-i|--imadir path] "
    << "[-m|--matchdir path] "
    << "[-o|--outdir path] "
    << "[-p|--pmvs 0 or 1] "
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
  IncrementalReconstructionEngine to3DEngine(sImaDirectory,
                                            sMatchesDir,
                                            sOutDir,
                                            true);

  if (to3DEngine.Process())
  {
    clock_t timeEnd = clock();
    std::cout << std::endl << " Ac-Sfm took : " << (timeEnd - timeStart) / CLOCKS_PER_SEC << " seconds." << std::endl;

    std::cout << std::endl << "Compute and export colorized point cloud" << std::endl;
    std::vector<Vec3> colortracks;
    to3DEngine.ColorizeTracks(colortracks);
    const reconstructorHelper & reconstructorHelperRef = to3DEngine.refToReconstructorHelper();
    reconstructorHelperRef.exportToPly(
      stlplus::create_filespec(sOutDir, "FinalColorized", ".ply"),
      &colortracks);

    // Export to openMVG format
    std::cout << std::endl << "Export 3D scene to openMVG format" << std::endl;
      reconstructorHelperRef.ExportToOpenMVGFormat(
        stlplus::folder_append_separator(sOutDir) + "SfM_output",
        to3DEngine.getFilenamesVector(),
        sImaDirectory,
        to3DEngine.getImagesSize(),
        to3DEngine.getTracks());

    // Manage export data to desired format
    // -> PMVS
    if (bPmvsExport)  {
      std::cout << std::endl << "Export 3D scene to PMVS format" << std::endl;
      reconstructorHelperRef.exportToPMVSFormat(
        stlplus::folder_append_separator(sOutDir) + "PMVS",
        to3DEngine.getFilenamesVector(),
        sImaDirectory);
    }

    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "\n Something goes wrong in the Structure from Motion process" << std::endl;
  }
  return EXIT_FAILURE;
}
