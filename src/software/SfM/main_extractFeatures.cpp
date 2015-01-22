
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"

/// Generic Image Collection image matching
#include "openMVG/matching_image_collection/Matcher_AllInMemory.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "software/SfM/pairwiseAdjacencyDisplay.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/indMatch_utils.hpp"

/// Feature detector and descriptor interface
#include "nonFree/sift/SIFT.hpp"

#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace std;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImaDirectory;
  std::string sImageNameFile;
  std::string sOutDir = "";
  bool bOctMinus1 = false;
  float dPeakThreshold = 0.04f;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('f', sImageNameFile, "imaNameFile"));
  cmd.add( make_option('s', bOctMinus1, "octminus1") );
  cmd.add( make_option('p', dPeakThreshold, "peakThreshold") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imadir path] \n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-f|--imaNameFile path" // TODO
      << "[-s|--octminus1 0 or 1] \n"
      << "[-p|--peakThreshold 0.04 -> 0.01] \n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imadir " << sImaDirectory << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--imaNameFile " << sImageNameFile << std::endl
            << "--octminus1 " << bOctMinus1 << std::endl
            << "--peakThreshold " << dPeakThreshold << std::endl;

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  std::vector<std::string> vec_imagePaths;  
  if(!sImageNameFile.empty())
  {    
    std::ifstream in(sImageNameFile.c_str());
    if(!in.is_open())  {
      std::cerr << "Impossible to read the specified file : " << sImageNameFile << std::endl;
      return false;
    }
    std::string sValue;
    while(std::getline(in, sValue))
    {
        std::string imagePath = stlplus::create_filespec(sImaDirectory, sValue);
        vec_imagePaths.push_back(imagePath);
    }
    in.close();
  }
  else  // Get all images in directory
  {
      std::vector<std::string> vec_fileNames = stlplus::folder_files(sImaDirectory);
      for( std::vector<std::string>::const_iterator it = vec_fileNames.begin(); it != vec_fileNames.end(); ++it)
      {
          std::string imagePath = stlplus::create_filespec(sImaDirectory, *it);
          vec_imagePaths.push_back(imagePath);
      }
  }
  //---------------------------------------
  // b. Compute features and descriptor
  //    - extract sift features and descriptor
  //    - if keypoints already computed, re-load them
  //    - else save features and descriptors on disk
  //---------------------------------------

  typedef Descriptor<unsigned char, 128> DescriptorT;
  typedef SIOPointFeature FeatureT;
  typedef std::vector<FeatureT> FeatsT;
  typedef vector<DescriptorT > DescsT;
  typedef KeypointSet<FeatsT, DescsT > KeypointSetT;

  {
    
    std::cout << "\n\n - EXTRACT FEATURES - " << std::endl;
    Image<unsigned char> imageGray;

    for(size_t i=0; i < vec_imagePaths.size(); ++i)  {
      Timer timer;
      std::cout << "Try image : " << vec_imagePaths[i] << std::endl;
      std::string sFeat = stlplus::create_filespec(sOutDir,
      stlplus::basename_part(vec_imagePaths[i]), "feat");
      std::string sDesc = stlplus::create_filespec(sOutDir,
      stlplus::basename_part(vec_imagePaths[i]), "desc");

      //If descriptors or features file are missing, compute them
      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc)) {
        if (!ReadImage(vec_imagePaths[i].c_str(), &imageGray))
        {
          std::cerr << "Cant read image" << std::endl;
          continue;
        }

        // Compute features and descriptors and export them to files
        KeypointSetT kpSet;
        SIFTDetector(imageGray,
          kpSet.features(), kpSet.descriptors(),
          bOctMinus1, true, dPeakThreshold);
        kpSet.saveToBinFile(sFeat, sDesc);
      }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
