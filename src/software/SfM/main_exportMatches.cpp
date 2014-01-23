
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"

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
  std::string sMatchFile;
  std::string sOutDir = "";

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('m', sMatchFile, "matchfile") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Export pairwise matches.\nUsage: " << argv[0] << ' '
      << "[-i|--imadir path] "
      << "[-d|--matchdir path] "
      << "[-m|--sMatchFile filename] "
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
  {
    std::ifstream in(stlplus::create_filespec(sMatchesDir, "lists", "txt").c_str());
    std::string sValue;
    while(in>>sValue)
    {
      int n = sValue.find_first_of(';');
      sValue = sValue.substr(0,n);
      vec_fileNames.push_back(sValue);
    }
    in.close();
  }
  if (vec_fileNames.empty()) {
    std::cerr << "\nEmpty input image list" << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Read matches
  //---------------------------------------

  typedef std::map< std::pair<size_t, size_t>, std::vector<IndMatch> > map_pairWiseMatches;
  map_pairWiseMatches map_Matches;
  PairedIndMatchImport(sMatchFile, map_Matches);

  // ------------
  // For each pair, export the matches
  // ------------

  //--
  //- Preprocess the images size
  Image<RGBColor> image;
  std::map< std::string, std::pair<size_t, size_t> > map_imageSize;
  for (std::vector<std::string>::const_iterator iterFilename = vec_fileNames.begin();
    iterFilename != vec_fileNames.end();
    ++iterFilename)
  {
    ReadImage( stlplus::create_filespec(sImaDirectory,*iterFilename).c_str() , &image);
    map_imageSize.insert(
      std::make_pair(*iterFilename,
      std::make_pair(image.Width(), image.Height())));
  }

  stlplus::folder_create(sOutDir);
  std::cout << "\n Export pairwise matches" << std::endl;
  C_Progress_display my_progress_bar( map_Matches.size() );
  for (map< std::pair<size_t, size_t>, std::vector<IndMatch> >::const_iterator iter = map_Matches.begin();
    iter != map_Matches.end();
    ++iter, ++my_progress_bar)
  {
    const size_t I = iter->first.first;
    const size_t J = iter->first.second;

    const std::pair<size_t, size_t> dimImage0 = map_imageSize.find(vec_fileNames[I])->second,
      dimImage1 = map_imageSize.find(vec_fileNames[J])->second;

    svgDrawer svgStream( dimImage0.first + dimImage1.first, max(dimImage0.second, dimImage1.second));
    svgStream.drawImage(stlplus::create_filespec(sImaDirectory,vec_fileNames[I]),
      dimImage0.first,
      dimImage0.second);
    svgStream.drawImage(stlplus::create_filespec(sImaDirectory,vec_fileNames[J]),
      dimImage1.first,
      dimImage1.second, dimImage0.first);

    const vector<IndMatch> & vec_FilteredMatches = iter->second;

    if (!vec_FilteredMatches.empty()) {
      // Load the features from the features files
      std::vector<SIOPointFeature> vec_featI, vec_featJ;
      loadFeatsFromFile(
        stlplus::create_filespec(sMatchesDir, stlplus::basename_part(vec_fileNames[I]), ".feat"),
        vec_featI);
      loadFeatsFromFile(
        stlplus::create_filespec(sMatchesDir, stlplus::basename_part(vec_fileNames[J]), ".feat"),
        vec_featJ);

      //-- Draw link between features :
      for (size_t i=0; i< vec_FilteredMatches.size(); ++i)  {
        const SIOPointFeature & imaA = vec_featI[vec_FilteredMatches[i]._i];
        const SIOPointFeature & imaB = vec_featJ[vec_FilteredMatches[i]._j];
        svgStream.drawLine(imaA.x(), imaA.y(),
          imaB.x()+dimImage0.first, imaB.y(), svgStyle().stroke("green", 2.0));
      }

      //-- Draw features (in two loop, in order to have the features upper the link, svg layer order):
      for (size_t i=0; i< vec_FilteredMatches.size(); ++i)  {
        const SIOPointFeature & imaA = vec_featI[vec_FilteredMatches[i]._i];
        const SIOPointFeature & imaB = vec_featJ[vec_FilteredMatches[i]._j];
        svgStream.drawCircle(imaA.x(), imaA.y(), imaA.scale(),
          svgStyle().stroke("yellow", 2.0));
        svgStream.drawCircle(imaB.x() + dimImage0.first, imaB.y(), imaB.scale(),
          svgStyle().stroke("yellow", 2.0));
      }
    }
    std::ostringstream os;
    os << stlplus::folder_append_separator(sOutDir)
      << iter->first.first << "_" << iter->first.second
      << "_" << iter->second.size() << "_.svg";
    ofstream svgFile( os.str().c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }
  return EXIT_SUCCESS;
}
