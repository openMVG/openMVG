// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::sfm;

// Convert from a SfM_Data format to another
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string
    sSfM_Data_Filename_In,
    sSfM_Data_Filename_Out;

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_switch('V', "VIEWS"));
  cmd.add(make_switch('I', "INTRINSICS"));
  cmd.add(make_switch('E', "EXTRINSICS"));
  cmd.add(make_switch('S', "STRUCTURE"));
  cmd.add(make_switch('C', "CONTROL_POINTS"));
  cmd.add(make_option('o', sSfM_Data_Filename_Out, "output_file"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
        << "[-i|--input_file] path to the input SfM_Data scene\n"
        << "[-o|--output_file] path to the output SfM_Data scene\n"
        << "\t .json, .bin, .xml, .ply, .baf\n"
        << "\n[Options to export partial data (by default all data are exported)]\n"
        << "\nUsable for json/bin/xml format"
        << "[-V|--VIEWS] export views\n"
        << "[-I|--INTRINSICS] export intrinsics (view orientations)\n"
        << "[-E|--EXTRINSICS] export extrinsics (view poses)\n"
        << "[-S|--STRUCTURE] export structure\n"
        << "[-C|--CONTROL_POINTS] export control points\n"
        << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sSfM_Data_Filename_In.empty() || sSfM_Data_Filename_Out.empty())
  {
    std::cerr << "Invalid input or output filename." << std::endl;
    return EXIT_FAILURE;
  }

  // OptionSwitch is cloned in cmd.add(),
  // so we must use cmd.used() instead of testing OptionSwitch.used
  int flags =
    (cmd.used('V') ? VIEWS      : 0)
  | (cmd.used('I') ? INTRINSICS : 0)
  | (cmd.used('E') ? EXTRINSICS : 0)
  | (cmd.used('S') ? STRUCTURE  : 0);

  flags = (flags) ? flags : ALL;

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Export the SfM_Data scene in the expected format
  if (Save(
    sfm_data,
    sSfM_Data_Filename_Out.c_str(),
    ESfM_Data(flags)))
  {
    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << std::endl
      << "An error occured while trying to save \"" << sSfM_Data_Filename_Out << "\"." << std::endl;
    return EXIT_FAILURE;
  }
}
