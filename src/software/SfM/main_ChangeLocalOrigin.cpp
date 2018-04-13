// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/stl/split.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"

#include <iomanip>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::geometry;


// Convert from a SfM_Data format to another
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string
    sSfM_Data_Filename_In,
    sOutDir,
    sLocalFrameOrigin;

  Vec3 local_Frame_Origin;
  bool b_first_frame_origin;

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_option('o', sOutDir, "output_dir"));
  cmd.add(make_option('l', sLocalFrameOrigin, "local_frame_origin"));
  cmd.add(make_switch('f', "first_frame_origin"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
        << "[-i|--input_file] path to the input SfM_Data scene\n"
        << "[-o|--output_dir] path to the output SfM_Data scene (in local frame)\n"
        << "[-l|--local_frame_origin] \"x;y;z\" of local frame origin\n"
        << "[-f|--first_frame_origin] use position of first frame as origin\n"
        << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())
  {
    std::cerr << std::endl
      << "No output SfM_Data filename specified." << std::endl;
    return EXIT_FAILURE;
  }


  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Local origin
  b_first_frame_origin = cmd.used('f');

  if (sLocalFrameOrigin.empty() && !b_first_frame_origin)
  {
    std::cerr << std::endl
      << "No local frame origin specified." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    if (b_first_frame_origin)
    {
      if (sfm_data.poses.empty()) {
        std::cerr << "The provided scene does not contain any camera poses." << std::endl;
        return EXIT_FAILURE;
      }
      local_Frame_Origin =  (sfm_data.poses.cbegin()->second).center();
    }
    else
    {
      std::vector<std::string> vec_str;
      stl::split(sLocalFrameOrigin, ';', vec_str);
      if (vec_str.size() != 3)
      {
        std::cerr << "\n Missing ';' character in local frame origin" << std::endl;
        return EXIT_FAILURE;
      }
      // Check that all local frame origin values are valid numbers
      for (size_t i = 0; i < vec_str.size(); ++i)
      {
        double readvalue = 0.0;
        std::stringstream ss;
        ss.str(vec_str[i]);
        if (! (ss >> readvalue) )  {
          std::cerr << "\n Used an invalid not a number character in local frame origin" << std::endl;
          return EXIT_FAILURE;
        }
        local_Frame_Origin[i] = readvalue;
      }
    }
  }

  std::cout << "Using frame origin: " << local_Frame_Origin.transpose() << std::endl;

  // Define the transformation (Will substract the local_Frame_Origin):
  Similarity3 sim( Pose3(Mat3::Identity(), local_Frame_Origin), 1.0);
  // Apply the transformation to the sfm_data scene
  const bool b_transform_priors = true;
  ApplySimilarity(sim, sfm_data, b_transform_priors);

  // Save changed sfm data
    //-- Export to disk computed scene (data & visualizable results)
  std::cout << "...Export SfM_Data to disk." << std::endl;
  if (!Save(sfm_data,
            stlplus::create_filespec(sOutDir, "sfm_data_local", ".bin"),
            ESfM_Data(ALL))
    || !Save(sfm_data,
             stlplus::create_filespec(sOutDir, "cloud_and_poses_local", ".ply"),
             ESfM_Data(ALL)))
  {
    std::cerr << "Cannot save the resulting sfm_data scene." << std::endl;
  }

  std::ofstream file_LocalFrameOrigin(stlplus::create_filespec(sOutDir, "local_frame_origin", ".txt"));
  file_LocalFrameOrigin << std::setprecision(8) << std::fixed;
  file_LocalFrameOrigin << local_Frame_Origin << "\n";
  file_LocalFrameOrigin.close();

  return EXIT_SUCCESS;
}
