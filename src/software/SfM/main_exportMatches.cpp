// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/sfm/sfm.hpp"

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
using namespace openMVG::features;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace svg;

// Convert HUE color to RGB
inline float hue2rgb(float p, float q, float t){
  if(t < 0) t += 1;
  if(t > 1) t -= 1;
  if(t < 1.f/6.f) return p + (q - p) * 6.f * t;
  if(t < 1.f/2.f) return q;
  if(t < 2.f/3.f) return p + (q - p) * (2.f/3.f - t) * 6.f;
  return p;
}

 //
 // Converts an HSL color value to RGB. Conversion formula
 // adapted from http://en.wikipedia.org/wiki/HSL_color_space.
 // Assumes h, s, and l are contained in the set [0, 1] and
 // returns r, g, and b in the set [0, 255].
 void hslToRgb(
   float h, float s, float l,
   unsigned char & r, unsigned char & g, unsigned char & b)
{
  if(s == 0){
    r = g = b = static_cast<unsigned char>(l * 255.f); // achromatic
  }else{
    const float q = l < 0.5f ? l * (1 + s) : l + s - l * s;
    const float p = 2.f * l - q;
    r = static_cast<unsigned char>(hue2rgb(p, q, h + 1.f/3.f) * 255.f);
    g = static_cast<unsigned char>(hue2rgb(p, q, h) * 255.f);
    b = static_cast<unsigned char>(hue2rgb(p, q, h - 1.f/3.f) * 255.f);
  }
}

int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sMatchFile;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('m', sMatchFile, "matchfile") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Export pairwise matches.\nUsage: " << argv[0] << "\n"
      << "[-i|--input_file file] path to a SfM_Data scene\n"
      << "[-d|--matchdir path]\n"
      << "[-m|--sMatchFile filename]\n"
      << "[-o|--outdir path]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }


  //---------------------------------------
  // Read SfM Scene (image view names)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the features
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!matches_provider->load(sfm_data, sMatchFile)) {
    std::cerr << "\nInvalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  // ------------
  // For each pair, export the matches
  // ------------

  stlplus::folder_create(sOutDir);
  std::cout << "\n Export pairwise matches" << std::endl;
  const Pair_Set pairs = matches_provider->getPairs();
  C_Progress_display my_progress_bar( pairs.size() );
  for (Pair_Set::const_iterator iter = pairs.begin();
    iter != pairs.end();
    ++iter, ++my_progress_bar)
  {
    const size_t I = iter->first;
    const size_t J = iter->second;

    const View * view_I = sfm_data.GetViews().at(I).get();
    const std::string sView_I= stlplus::create_filespec(sfm_data.s_root_path,
      view_I->s_Img_path);
    const View * view_J = sfm_data.GetViews().at(J).get();
    const std::string sView_J= stlplus::create_filespec(sfm_data.s_root_path,
      view_J->s_Img_path);

    const std::pair<size_t, size_t>
      dimImage_I = std::make_pair(view_I->ui_width, view_I->ui_height),
      dimImage_J = std::make_pair(view_J->ui_width, view_J->ui_height);

    svgDrawer svgStream( dimImage_I.first + dimImage_J.first, max(dimImage_I.second, dimImage_J.second));
    svgStream.drawImage(sView_I,
      dimImage_I.first,
      dimImage_I.second);
    svgStream.drawImage(sView_J,
      dimImage_J.first,
      dimImage_J.second, dimImage_I.first);

    const vector<IndMatch> & vec_FilteredMatches = matches_provider->_pairWise_matches.at(*iter);

    if (!vec_FilteredMatches.empty()) {

      const PointFeatures & vec_feat_I = feats_provider->getFeatures(view_I->id_view);
      const PointFeatures & vec_feat_J = feats_provider->getFeatures(view_J->id_view);

      //-- Draw link between features :
      for (size_t i=0; i< vec_FilteredMatches.size(); ++i)  {
        const PointFeature & imaA = vec_feat_I[vec_FilteredMatches[i]._i];
        const PointFeature & imaB = vec_feat_J[vec_FilteredMatches[i]._j];
        // Compute a flashy colour for the correspondence
        unsigned char r,g,b;
        hslToRgb( (rand()%360) / 360., 1.0, .5, r, g, b);
        std::ostringstream osCol;
        osCol << "rgb(" << (int)r <<',' << (int)g << ',' << (int)b <<")";
        svgStream.drawLine(imaA.x(), imaA.y(),
          imaB.x()+dimImage_I.first, imaB.y(), svgStyle().stroke(osCol.str(), 2.0));
      }

      //-- Draw features (in two loop, in order to have the features upper the link, svg layer order):
      for (size_t i=0; i< vec_FilteredMatches.size(); ++i)  {
        const PointFeature & imaA = vec_feat_I[vec_FilteredMatches[i]._i];
        const PointFeature & imaB = vec_feat_J[vec_FilteredMatches[i]._j];
        svgStream.drawCircle(imaA.x(), imaA.y(), 3.0,
          svgStyle().stroke("yellow", 2.0));
        svgStream.drawCircle(imaB.x() + dimImage_I.first, imaB.y(), 3.0,
          svgStyle().stroke("yellow", 2.0));
      }
    }
    std::ostringstream os;
    os << stlplus::folder_append_separator(sOutDir)
      << iter->first << "_" << iter->second
      << "_" << vec_FilteredMatches.size() << "_.svg";
    ofstream svgFile( os.str().c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }
  return EXIT_SUCCESS;
}
