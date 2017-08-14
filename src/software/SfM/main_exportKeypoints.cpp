// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/svg_features.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_view.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Export pairwise matches.\nUsage: " << argv[0] << "\n"
      << "[-i|--input_file file] path to a SfM_Data scene\n"
      << "[-d|--matchdir path]\n"
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
  using namespace openMVG::features;
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

  // ------------
  // For each image, export visually the keypoints
  // ------------

  stlplus::folder_create(sOutDir);
  std::cout << "\n Export extracted keypoints for all images" << std::endl;
  C_Progress_display my_progress_bar( sfm_data.views.size() );
  for (Views::const_iterator iterViews = sfm_data.views.begin();
      iterViews != sfm_data.views.end();
      ++iterViews, ++my_progress_bar)
  {
    const View * view = iterViews->second.get();
    const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
      view->s_Img_path);

    std::ostringstream os;
    os << stlplus::folder_append_separator(sOutDir)
      << stlplus::basename_part(sView_filename)
      << "_" << feats_provider->getFeatures(view->id_view).size() << "_.svg";

    Features2SVG
    (
      sView_filename,
      {view->ui_width, view->ui_height},
      feats_provider->getFeatures(view->id_view),
      os.str()
    );
  }
  return EXIT_SUCCESS;
}
