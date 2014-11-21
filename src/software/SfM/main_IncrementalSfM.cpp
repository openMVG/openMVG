
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "software/SfM/SfMIncrementalEngine.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/system/timer.hpp"

using namespace openMVG;

/// From 2 given image file-names, find the two corresponding index in the openMVG 'lists.txt' file.
bool computeIndexFromImageNames(
  const std::string& sMatchesDir,
  const std::pair<std::string,std::string>& initialPairName,
  std::pair<size_t, size_t>& initialPairIndex)
{
  const std::string sListsFile = stlplus::create_filespec(sMatchesDir, "lists.txt" );
  if (!stlplus::is_file(sListsFile)) {
    std::cerr << "\nCannot access input file \""<< sListsFile << "\"" << std::endl;
    return false;
  }

  std::vector<std::string> vec_camImageName;
  if (!openMVG::SfMIO::loadImageList(vec_camImageName, sListsFile, false))
  {
    std::cerr << "\nEmpty or invalid image list." << std::endl;
    return false;
  }

  if (initialPairName.first == initialPairName.second)
  {
    std::cerr << "\nInvalid image names. You cannot use the same image to create a pair." << std::endl;
    return false;
  }

  std::vector<std::string>::const_iterator imageName = find(vec_camImageName.begin(), vec_camImageName.end(), initialPairName.first);
  if(imageName == vec_camImageName.end())
  {
      std::cerr << "\nCannot access to the specified image: \""<< initialPairName.first << "\"" << std::endl;
      return false;
  }
  initialPairIndex.first = std::distance<std::vector<std::string>::const_iterator>(vec_camImageName.begin(), imageName);  
  imageName = find(vec_camImageName.begin(), vec_camImageName.end(), initialPairName.second);
  if(imageName == vec_camImageName.end())
  {
      std::cerr << "\nCannot access to the specified image: \""<< initialPairName.second << "\"" << std::endl;
      return false;
  }
  initialPairIndex.second = std::distance<std::vector<std::string>::const_iterator>(vec_camImageName.begin(), imageName);
  return true;
}


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
  bool bRefinePPandDisto = true;
  bool bRefineFocal = true;
  bool bColoredPointCloud = false;
  std::pair<std::string,std::string> initialPair("","");

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('p', bPmvsExport, "pmvs") );
  cmd.add( make_option('a', initialPair.first, "initialPairA") );
  cmd.add( make_option('b', initialPair.second, "initialPairB") );
  cmd.add( make_option('c', bColoredPointCloud, "coloredPointCloud") );
  cmd.add( make_option('d', bRefinePPandDisto, "refinePPandDisto") );
  cmd.add( make_option('f', bRefineFocal, "refineFocal") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir PATH] \n"
    << "[-m|--matchdir PATH] \n"
    << "[-o|--outdir PATH] \n"
    << "[-p|--pmvs 0 or 1] \n"
    << "[-a|--initialPairA NAME] \n"
    << "[-b|--initialPairB NAME] \n"
    << "[-c|--coloredPointCloud 0(default) or 1]\n"
    << "[-d|--refinePPandDisto \n"
    << "\t 0-> refine only the Focal,\n"
    << "\t 1-> refine Focal, Principal point and radial distortion factors.] \n"
    << "[-f|--refineFocal \n"
    << "\t 0-> refine only Principal point and radial distortion,\n"
    << "\t 1-> refine Focal, Principal point and radial distortion ] \n"
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

  openMVG::Timer timer;
  IncrementalReconstructionEngine to3DEngine(sImaDirectory,
                                            sMatchesDir,
                                            sOutDir,
                                            true);

  // Handle Initial pair parameter
  if (!initialPair.first.empty() && !initialPair.second.empty())
  {
    std::pair<size_t, size_t> initialPairIndex(0,0);
    if(!computeIndexFromImageNames(sMatchesDir, initialPair, initialPairIndex))
      return EXIT_FAILURE;
    to3DEngine.setInitialPair(initialPairIndex);
  }  
  to3DEngine.setIfRefinePrincipalPointAndRadialDisto(bRefinePPandDisto);
  to3DEngine.setIfRefineFocal(bRefineFocal);

  if (to3DEngine.Process())
  {
    clock_t timeEnd = clock();
    std::cout << std::endl << " Ac-Sfm took (s): " << timer.elapsed() << "." << std::endl;

    const reconstructorHelper & reconstructorHelperRef = to3DEngine.refToReconstructorHelper();
    std::vector<Vec3> vec_tracksColor;
    if (bColoredPointCloud)
    {
      // Compute the color of each track
      to3DEngine.ColorizeTracks(vec_tracksColor);
    }
    reconstructorHelperRef.exportToPly(
      stlplus::create_filespec(sOutDir, "FinalColorized", ".ply"),
      bColoredPointCloud ? &vec_tracksColor : NULL);

    // Export to openMVG format
    std::cout << std::endl << "Export 3D scene to openMVG format" << std::endl
      << " -- Point cloud color: " << (bColoredPointCloud ? "ON" : "OFF") << std::endl;

    if (!reconstructorHelperRef.ExportToOpenMVGFormat(
      stlplus::folder_append_separator(sOutDir) + "SfM_output",
      to3DEngine.getFilenamesVector(),
      sImaDirectory,
      to3DEngine.getImagesSize(),
      to3DEngine.getTracks(),
      bColoredPointCloud ? &vec_tracksColor : NULL,
      true,
      std::string("generated by the Sequential OpenMVG Calibration Engine")))
    {
      std::cerr << "Error while saving the scene." << std::endl;
    }

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
