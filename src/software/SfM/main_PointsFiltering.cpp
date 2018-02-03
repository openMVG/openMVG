// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Filter Structure based on statistics computed per view    :\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sfm_data_filename;
  std::string sfm_data_filename_out;
  double factor = 5.2; // 5.2 * median ~= X84

  cmd.add( make_option('i', sfm_data_filename, "input_file") );
  cmd.add( make_option('o', sfm_data_filename_out, "output_file") );
  cmd.add( make_option('f', factor, "factor") );

  try
  {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string& s)
  {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] path to a SfM_Data scene\n"
      << "[-o|--output_file] path where the filtered SfM_data scene will be saved\n"
      << "\n[Optional]\n"
      << "[-f|--factor] factor apply on the median depth per view for thresholding\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sfm_data_filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sfm_data_filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Structure filtering
  //---------------------------------------

  const IndexT k_min_point_per_pose = 12;
  const IndexT k_min_track_length = 2;
  const double k_min_median_depth = DepthCleaning(
    sfm_data,
    factor,
    k_min_point_per_pose,
    k_min_track_length
  );
  std::cout << "MIN MEDIAN DEPTH VALUE = " << k_min_median_depth << std::endl;

  if (!Save(sfm_data, sfm_data_filename_out, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The output SfM_Data file \""<< sfm_data_filename_out << "\" cannot be saved." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
