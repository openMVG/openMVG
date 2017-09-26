// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_graph_utils.hpp"
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
  std::cout << std::endl;

  CmdLine cmd;

  std::string sfm_data_filename;
  std::string match_filename;
  std::string match_component_filename;
  bool is_biedge = false;
  int min_nodes = 3;

  cmd.add(make_option('i', sfm_data_filename, "input_file"));
  cmd.add(make_option('m', match_filename, "match_file"));
  cmd.add(make_option('o', match_component_filename, "match_component_file"));
  cmd.add(make_option('b', is_biedge, "biedge"));
  cmd.add(make_option('n', min_nodes, "min_nodes"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] path to a SfM_Data scene\n"
      << "[-m|--match_file] path to the matches that corresponds to the provided SfM_Data scene\n"
      << "[-o|--match_component_file] path to the matches components that corresponds to the provided SfM_Data scene\n"
      << "\n[Optional]\n"
      << "[-b|--biedge]\n"
      << "[-n|--min_nodes] Note:value of n should larger than 3\n"
      << std::endl;
    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::file_exists(sfm_data_filename))
  {
    return EXIT_FAILURE;
  }

  if (!stlplus::file_exists(match_filename))
  {
    return EXIT_FAILURE;
  }

  const std::string &match_component_dir = stlplus::folder_part(match_component_filename);
  if (!stlplus::folder_exists(match_component_dir))
  {
    const bool folder_create_flag = stlplus::folder_create(match_component_dir);
    if (!folder_create_flag)
    {
      std::cerr << "Cannot create the output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sfm_data_filename, ESfM_Data(VIEWS | INTRINSICS)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sfm_data_filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Matches reading
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!(matches_provider->load(sfm_data, match_filename)))
  {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  // Split match_filename by connected components;
  const bool success_flag = SplitMatchFileIntoMatchFiles(sfm_data, match_filename,
    match_component_filename, is_biedge, min_nodes);
  if (success_flag)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
