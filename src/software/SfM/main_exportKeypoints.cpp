
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"

#include "software/SfM/SfMIOHelper.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <map>

using namespace openMVG;
using namespace openMVG::matching;
using namespace svg;


int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sImaDirectory;
  std::string sMatchesDir;
  std::string sOutDir = "";

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Export pairwise matches.\nUsage: " << argv[0] << ' '
      << "[-i|--imadir path] "
      << "[-d|--matchdir path] "
      << "[-o|--outdir path] "
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }


  //---------------------------------------
  // Read images names
  //---------------------------------------

  std::vector<std::string> vec_fileNames;  
  if (!SfMIO::loadImageList( vec_fileNames,
      stlplus::create_filespec(sMatchesDir, "lists", "txt"),false)) {
    std::cerr << "\nEmpty input image list" << std::endl;
    return EXIT_FAILURE;
  }


  // ------------
  // For each image, export visually the keypoints
  // ------------

  stlplus::folder_create(sOutDir);
  std::cout << "\n Export extracted keypoints for all images" << std::endl;
  C_Progress_display my_progress_bar( vec_fileNames.size() );
  for (std::vector<std::string>::const_iterator iterFilename = vec_fileNames.begin();
    iterFilename != vec_fileNames.end();
    ++iterFilename, ++my_progress_bar)
  {
    const size_t I = std::distance(
      (std::vector<std::string>::const_iterator)vec_fileNames.begin(),
      iterFilename);

    //--
    //- Read the image size
    Image<RGBColor> image;
    ReadImage( stlplus::create_filespec(sImaDirectory,*iterFilename).c_str() , &image);
    const std::pair<size_t, size_t> dimImage (image.Width(),image.Height());

    svgDrawer svgStream( dimImage.first, dimImage.second);
    svgStream.drawImage(stlplus::create_filespec(sImaDirectory,*iterFilename),
      dimImage.first,
      dimImage.second);

    // Load the features from the feature file
    std::vector<SIOPointFeature> vec_feat;
    loadFeatsFromFile(
      stlplus::create_filespec(sMatchesDir, stlplus::basename_part(*iterFilename), ".feat"),
      vec_feat);

    //-- Draw features
    for (size_t i=0; i< vec_feat.size(); ++i)  {
      const SIOPointFeature & feature = vec_feat[i];
      svgStream.drawCircle(feature.x(), feature.y(), feature.scale(),
          svgStyle().stroke("yellow", 2.0));
    }

    // Write the SVG file
    std::ostringstream os;
    os << stlplus::folder_append_separator(sOutDir)
      << stlplus::basename_part(*iterFilename)
      << "_" << vec_feat.size() << "_.svg";
    ofstream svgFile( os.str().c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }
  return EXIT_SUCCESS;
}
