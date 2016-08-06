// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/utils/alignment.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>
#include <sstream>
#include <vector>

using namespace openMVG;
using namespace openMVG::sfm;

static bool parseAlignScale(const std::string& alignScale, double& S, Mat3& R, Vec3& t);

// Convert from a SfM_Data format to another
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename_In;
  std::string sSfM_Data_Filename_InRef;
  std::string sSfM_Data_Filename_Out;
  std::string sSfm_Data_YAlignScale;

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_option('r', sSfM_Data_Filename_InRef, "reference_file"));
  cmd.add(make_option('y', sSfm_Data_YAlignScale, "y_align_scale"));
  cmd.add(make_option('o', sSfM_Data_Filename_Out, "output_file"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
        << "[-i|--input_file] path to the input SfM_Data scene to align.\n"
        << "[-r|--reference_file] path to the scene used as the reference coordinate system\n"
        << "[-y|--y_align_scale] align [X,Y,Z] to +Y-axis, scale scene by S; syntax: X,Y,Z;S \n"
        << "[-o|--output_file] path to the output SfM_Data scene\n"
        << "\t .json, .bin, .xml, .ply, .baf"
#if HAVE_ALEMBIC
           ", .abc"
#endif
           "\n"
        << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }
  

  if (sSfM_Data_Filename_In.empty() ||
      sSfM_Data_Filename_Out.empty())
  {
    std::cerr << "Invalid input or output filename." << std::endl;
    return EXIT_FAILURE;
  }

  if (sSfM_Data_Filename_InRef.empty() &&
      sSfm_Data_YAlignScale.empty())
  {
    std::cerr << "At least one of -y and -r must be specified." << std::endl;
    return EXIT_FAILURE;
  }
  
  if (!sSfM_Data_Filename_InRef.empty() &&
      !sSfm_Data_YAlignScale.empty())
  {
    std::cerr << "Must specify exactly one of alignment and reference scene." << std::endl;
    return EXIT_FAILURE;
  }

  // Load input scene
  SfM_Data sfm_data_in;
  if (!Load(sfm_data_in, sSfM_Data_Filename_In, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  double S;
  Mat3 R;
  Vec3 t;

  if (!sSfM_Data_Filename_InRef.empty())
  {
    // Load reference scene
    SfM_Data sfm_data_inRef;
    if (!Load(sfm_data_inRef, sSfM_Data_Filename_InRef, ESfM_Data(ALL)))
    {
      std::cerr << std::endl
        << "The reference SfM_Data file \"" << sSfM_Data_Filename_InRef << "\" cannot be read." << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "Search similarity transformation." << std::endl;
    bool hasValidSimilarity = computeSimilarity(sfm_data_in, sfm_data_inRef, &S, &R, &t);
    if(!hasValidSimilarity)
    {
      std::cerr << std::endl
        << "Failed to find similarity between the 2 SfM scenes:"
        << "\"" << sSfM_Data_Filename_In << "\", "
        << "\"" << sSfM_Data_Filename_InRef << "\""
        << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (!parseAlignScale(sSfm_Data_YAlignScale, S, R, t))
  {
    std::cerr << std::endl << "Failed to parse align/scale argument.";
    return EXIT_FAILURE;
  }
  
  std::cout << "Apply transformation:" << std::endl;
  std::cout << " - Scale: " << S << std::endl;
  std::cout << " - Rotation:\n" << R << std::endl;
  std::cout << " - Translate: " << t.transpose() << std::endl;

  applyTransform(sfm_data_in, S, R, t);
  
  std::cout << "Save into \"" << sSfM_Data_Filename_Out << "\"" << std::endl;
  
  // Export the SfM_Data scene in the expected format
  if (!Save(sfm_data_in, sSfM_Data_Filename_Out, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
      << "An error occurred while trying to save \"" << sSfM_Data_Filename_Out << "\"." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

static bool parseAlignScale(const std::string& alignScale, double& S, Mat3& R, Vec3& t)
{
  double rx, ry, rz;
  
  {
    char delim[3];
    std::istringstream iss(alignScale);
    if (!(iss >> rx >> delim[0] >> ry >> delim[1] >> rz >> delim[2] >> S))
      return false;
    if (delim[0] != ',' || delim[1] != ',' || delim[2] != ';')
      return false;
  }
  
  auto q = ::openMVG::Quaternion::FromTwoVectors(Vec3(rx, ry, rz), Vec3(0, 1, 0));
  R = q.matrix();

  t = Vec3::Zero();
  
  return true;
}
