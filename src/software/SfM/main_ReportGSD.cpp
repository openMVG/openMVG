// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_gsd.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"


using namespace openMVG;
using namespace openMVG::sfm;

// Convert from a SfM_Data format to another
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string
    sSfM_Data_Filename_In;

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n' 
        << "[-i|--input_file] path to the input SfM_Data scene\n";

      OPENMVG_LOG_ERROR << s;
      return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
  {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  openMVG::sfm::View *view = sfm_data.GetViews().at(0).get();
  std::shared_ptr<cameras::IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic);
  double ccdw = intrinsic->ccdw();

  // last argument can take both ccd dimensions, it'll take the max value from the calculation 
  // then as GSD should, format is h,w
  double gsd = get_ground_sampling_distance(0, sfm_data, Vec2(0.0,ccdw));

  OPENMVG_LOG_INFO << "Scene GSD is: " << gsd;

  return EXIT_FAILURE;
}
