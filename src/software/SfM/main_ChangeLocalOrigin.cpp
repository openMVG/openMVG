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
    sLocalFrameOrigin,
    sTransformationMatrix;

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_option('o', sOutDir, "output_dir"));
  cmd.add(make_option('l', sLocalFrameOrigin, "local_frame_origin"));
  cmd.add(make_option('T', sTransformationMatrix, "transformation_matrix"));
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
        << "[-T|--transformation_matrix] \"t11;t12;t13;t14...t41;t42;t43;t44\" row-wise elements of transformation matrix cur_T_prev\n"
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

  // Transformation
  Similarity3 sim_transform;
  bool b_first_frame_origin = cmd.used('f');
  bool b_local_frame_origin = !sLocalFrameOrigin.empty();
  bool b_transformation_matrix = !sTransformationMatrix.empty();

  if (!(b_first_frame_origin || b_local_frame_origin || b_transformation_matrix))
  {
    std::cerr << std::endl
      << "No operation specified." << std::endl;
    return EXIT_FAILURE;
  }

  if (b_first_frame_origin)
  {
    // Use first pose as new origin
    if (sfm_data.poses.empty()) {
      std::cerr << "The provided scene does not contain any camera poses." << std::endl;
      return EXIT_FAILURE;
    }
    Vec3 local_frame_origin =  (sfm_data.poses.cbegin()->second).center();
    sim_transform = Similarity3(Pose3(Mat3::Identity(), local_frame_origin), 1.0);

    std::cout << "Using first frame as origin: " << local_frame_origin.transpose() << std::endl;
  }
  else if(b_local_frame_origin)
  {
    // Using provided position as origin
    Vec3 local_frame_origin;
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
      local_frame_origin[i] = readvalue;
    }
    // Create similarity transform
    sim_transform = Similarity3(Pose3(Mat3::Identity(), local_frame_origin), 1.0);
    std::cout << "Transform using position as origin: " << local_frame_origin.transpose() << std::endl;
  }
  else if(b_transformation_matrix)
  {
    // Using provided transformation
    Mat4 T_matrix;

    std::vector<std::string> vec_str;
    stl::split(sTransformationMatrix, ';', vec_str);
    if (vec_str.size() != 16)
    {
      std::cerr << "\n Missing ';' character in the transformation matrix" << std::endl;
      return EXIT_FAILURE;
    }
    // Check that all transformation matrix values are valid numbers
    for (size_t i = 0; i < vec_str.size(); ++i)
    {
      double readvalue = 0.0;
      std::stringstream ss;
      ss.str(vec_str[i]);
      if (! (ss >> readvalue) )  {
        std::cerr << "\n Used an invalid not a number character in the transformation matrix" << std::endl;
        return EXIT_FAILURE;
      }
      T_matrix((i/4),(i%4)) = readvalue;
    }

    // Check if x,y,z dimension scale values are equal (requirement!)
    Mat3 rot_matrix = T_matrix.block(0,0,3,3);
    Vec3 center = -rot_matrix.transpose() * T_matrix.block(0,3,3,1);
    double scale_factor = -1;
    for (int col_i = 0; col_i < 3; ++col_i)
    {
      double scale_col = rot_matrix.block(0,col_i,3,1).norm();
      if (scale_factor > 0 && fabs(scale_factor-scale_col) > 1e-4)
      {
        std::cout << "S: " << scale_factor << " S2: " << scale_col<<"\n";
        std::cerr << "\n Scale factor for X,Y,Z dimensions is not unified" << std::endl;
        return EXIT_FAILURE;
      }
      if (scale_factor <=0)
      {
        scale_factor = scale_col;
      }
    }

    // Get rotation matrix without scale
    for (int col_i = 0; col_i < 3; ++col_i)
    {
      rot_matrix.block(0,col_i,3,1) = rot_matrix.block(0,col_i,3,1) / rot_matrix.block(0,col_i,3,1).norm();
    }

    // Create similarity transform
    sim_transform = Similarity3(Pose3(rot_matrix, center), scale_factor);

    // Print info
    std::cout << "Transform using T as transformation matrix:\n" << T_matrix << std::endl;
  }
  
  // Apply the transformation to the sfm_data scene
  const bool b_transform_priors = true;
  ApplySimilarity(sim_transform, sfm_data, b_transform_priors);

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

  // Output detail of transform

  if (b_first_frame_origin || b_local_frame_origin)
  {
    std::ofstream file_LocalFrameOrigin(stlplus::create_filespec(sOutDir, "local_frame_origin", ".txt"));
    file_LocalFrameOrigin << std::setprecision(8) << std::fixed;
    file_LocalFrameOrigin << sim_transform.pose_.center() << "\n";
    file_LocalFrameOrigin.close();
  }
  else if(b_transformation_matrix)
  {
    std::ofstream file_LocalFrameOrigin(stlplus::create_filespec(sOutDir, "local_frame_origin", ".txt"));
    file_LocalFrameOrigin << std::setprecision(8) << std::fixed;
    file_LocalFrameOrigin << sim_transform.pose_.center() << "\n";
    file_LocalFrameOrigin.close();

    Mat4 T_matrix = Mat4::Identity();
    T_matrix.block(0,0,3,3) = sim_transform.pose_.rotation() * sim_transform.scale_;
    T_matrix.block(0,3,3,1) = sim_transform.pose_.translation();

    std::ofstream file_T(stlplus::create_filespec(sOutDir, "applied_transformation", ".txt"));
    file_T << std::setprecision(8) << std::fixed;
    file_T << T_matrix << "\n";
    file_T.close();

  }

  return EXIT_SUCCESS;
}
