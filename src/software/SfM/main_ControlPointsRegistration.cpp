// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#include "software/SfM/SfMPlyHelper.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/cmdLine/cmdLine.h"

#include <iostream>
#include <iomanip>
#include <string>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

using namespace std;

int main(int argc, char **argv)
{
  std::string
    sSfM_Data_Filename_In,
    sSfM_Data_Filename_Out;

  CmdLine cmd;
  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_option('o', sSfM_Data_Filename_Out, "output_file"));

  try
  {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string& s)
  {
    std::cerr
      << "Usage: " << argv[0] << '\n'
      << " Control Points registration of a SfM Data scene,\n"
      << "[-i|--input_file] path to the input SfM_Data scene\n"
      << "[-o|--output_file] path to the output SfM_Data scene\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sSfM_Data_Filename_In.empty() || sSfM_Data_Filename_Out.empty())
  {
    std::cerr << "Invalid input or output filename." << std::endl;
    return EXIT_FAILURE;
  }
  
  //
  // Load a SfM scene
  // - Compute the registration between the control point observations and reference position

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
  {
    std::cerr
      << "\nThe input SfM_Data file \"" << sSfM_Data_Filename_In
      << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---
  // registration (coarse):
  // - compute the 3D points corresponding to the control point observation for the SfM scene
  // - compute a coarse registration between the controls points & the triangulated point
  // - transform the scene according the found transformation
  //---
  std::map<IndexT, Vec3> vec_control_points, vec_triangulated;
  std::map<IndexT, double> vec_triangulation_errors;
  for (const auto & control_point_it : sfm_data.control_points)
  {
    const Landmark & landmark = control_point_it.second;
    //Triangulate the observations:
    const Observations & obs = landmark.obs;

    if (obs.empty()) // Error if this control point has no observation (never seen in any of the images)
    {
      std::cout << "Control point has no observations, ignoring control point" << std::endl;
      continue;
    }

    std::vector<Vec3> bearing;
    std::vector<Mat34> poses;
    bearing.reserve(obs.size());
    poses.reserve(obs.size());
    for (const auto & obs_it : obs)
    {
      const View * view = sfm_data.views.at(obs_it.first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const openMVG::geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
      const Vec2 pt = obs_it.second.x;
      bearing.emplace_back((*cam)(cam->get_ud_pixel(pt)));
      poses.emplace_back(pose.asMatrix());
    }
    const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
    Vec4 Xhomogeneous;
    if (!TriangulateNViewAlgebraic(bearing_matrix, poses, &Xhomogeneous))
    {
      std::cout << "Invalid triangulation, ignoring control point" << std::endl;
      continue;
    }
    const Vec3 X = Xhomogeneous.hnormalized();
    // Test validity of the hypothesis (front of the cameras):
    bool bChierality = true;
    int i(0);
    double reprojection_error_sum(0.0);
    for (const auto & obs_it : obs)
    {
      const View * view = sfm_data.views.at(obs_it.first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      bChierality &= CheiralityTest(bearing[i], pose, X);
      const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Vec2 pt = obs_it.second.x;
      const Vec2 residual = cam->residual(pose(X), pt);
      reprojection_error_sum += residual.norm();
      ++i;
    }
    if (bChierality) // Keep the point only if it has a positive depth
    {
      vec_triangulated[control_point_it.first] = X;
      vec_control_points[control_point_it.first] = landmark.X;
      vec_triangulation_errors[control_point_it.first] = reprojection_error_sum / (double)bearing.size();
    }
    else
    {
      std::cout << "Control Point cannot be triangulated (not in front of the cameras)" << std::endl;
    }
  }

  if (vec_control_points.size() < 3)
  {
    std::cerr << "Insufficient number of triangulated control points." << std::endl;
    return EXIT_FAILURE;
  }

  // compute the similarity
  {
    // data conversion to appropriate container
    Mat x1(3, vec_control_points.size()),
      x2(3, vec_control_points.size());

    size_t i = 0;
    for (const auto & control_point_it : sfm_data.control_points)
    {
      const IndexT CPIndex = control_point_it.first;
      if (vec_triangulated.find(CPIndex) != vec_triangulated.end())
      {
        x1.col(i) = vec_triangulated[CPIndex];
        x2.col(i) = vec_control_points[CPIndex];
        i++;
      }
    }

    std::cout
      << "Control points observation triangulations:\n"
      << x1 << std::endl << std::endl
      << "Control points coords:\n"
      << x2 << std::endl << std::endl;

    Vec3 t;
    Mat3 R;
    double S;
    if (openMVG::geometry::FindRTS(x1, x2, &S, &t, &R))
    {
      openMVG::geometry::Refine_RTS(x1, x2, &S, &t, &R);
      std::cout << "Found transform:\n"
        << " scale: " << S << "\n"
        << " rotation:\n" << R << "\n"
        << " translation: " << t.transpose() << std::endl;


      //--
      // Apply the found transformation as a 3D Similarity transformation matrix // S * R * X + t
      //--

      const openMVG::geometry::Similarity3 sim(geometry::Pose3(R, -R.transpose() * t / S), S);
      openMVG::sfm::ApplySimilarity(sim, sfm_data);

      // Display some statistics:
      std::stringstream os;
      for (Landmarks::const_iterator iterL = sfm_data.control_points.begin();
        iterL != sfm_data.control_points.end(); ++iterL)
      {
        const IndexT CPIndex = iterL->first;
        // If the control point has not been used, continue...
        if (vec_triangulation_errors.find(CPIndex) == vec_triangulation_errors.end())
          continue;

        os
          << "CP index: " << CPIndex << "\n"
          << "CP triangulation error: " << vec_triangulation_errors[CPIndex] << " pixel(s)\n"
          << "CP registration error: "
          << (sim(vec_triangulated[CPIndex]) - vec_control_points[CPIndex]).norm() << " user unit(s)" << "\n\n";
      }
      std::cout << os.str();
    }
    else
    {
      std::cerr << "Registration failed. Please check your Control Points coordinates." << std::endl;
      return EXIT_FAILURE;
    }
  }

  //---
  // Bundle adjustment with GCP
  //---
  {
    using namespace openMVG::sfm;
    Bundle_Adjustment_Ceres::BA_Ceres_options options;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    Control_Point_Parameter control_point_opt(20.0, true);
    if (!bundle_adjustment_obj.Adjust(sfm_data,
          Optimize_Options
          (
            cameras::Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
            Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
            Structure_Parameter_Type::ADJUST_ALL, // Adjust structure
            control_point_opt // Use GCP and weight more their observation residuals
          )
        )
      )
    {
      std::cerr << "BA with GCP failed." << std::endl;
      return EXIT_FAILURE;
    }
  }


  // Export the SfM_Data scene in the expected format
  if (Save(
        sfm_data,
        sSfM_Data_Filename_Out.c_str(),
        ESfM_Data(ALL)))
  {
    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr
      << std::endl
      << "An error occured while trying to save \""
      << sSfM_Data_Filename_Out << "\"." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_FAILURE;
}
