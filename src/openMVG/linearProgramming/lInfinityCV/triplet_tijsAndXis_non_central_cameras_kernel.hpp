// Copyright (c) 2016 Pierre MOULON, Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_FromXi_Ri_non_central_cameras.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"

namespace openMVG {
namespace non_central_camera {

/// A triplet of 3 poses (represent the pose of the 3 non central cameras)
struct PoseTripletErrorModel
{
  Mat3  R1, R2, R3;
  Vec3  t1, t2, t3;

  static double Error
  (
    const PoseTripletErrorModel & t,
    const std::vector < std::vector <double> > & pointInfo,
    const std::vector<Mat3> & vec_local_rotations, // rig subcamera rotations
    const std::vector<Vec3> & vec_local_centers // rig subcamera translations
  )
  {
    // List corresponding observations & projection matrices (for triangulation)
    std::vector < std::pair <Mat34, Vec2> >  views;

    for( size_t i = 0 ; i < pointInfo.size() ; ++ i)
    {
      // extract subcamera rotations and translation
      const size_t I = (size_t) pointInfo[i][2]; // intrinsic id
      const Mat3 RI = vec_local_rotations[I];
      const Vec3 tI = -RI * vec_local_centers[I];

      // compute projection matrix
      Mat34 P ;
      switch( (int) pointInfo[i][3] )  // pose Id
      {
        // if first pose
        case 0 :
          P = HStack(RI * t.R1, RI * t.t1 + tI);
          break;

        // if second pose
        case 1 :
          P = HStack(RI * t.R2, RI * t.t2 + tI);
          break;

        // if third pose
        case 2 :
          P = HStack(RI * t.R3, RI * t.t3 + tI);
          break;
        default:
          std::cerr << "Invalid index" << std::endl;
      };

      views.emplace_back(P, Vec2(pointInfo[i][0], pointInfo[i][1]));
    }

    // Setup the observation to triangulate
    Triangulation triangulationObj;
    for( size_t i = 0 ; i < views.size(); ++i )
      triangulationObj.add ( views[i].first, views[i].second );

    const Vec3 X = triangulationObj.compute();

    //- Return error
    double max_error = 0.0;

    // Compute the LInfinity squared residual & return it
    for( size_t i = 0 ; i < views.size(); ++i )
    {
      const Mat34 & P = views[i].first;
      const Vec2  & pt = views[i].second;
      max_error = std::max( (Project(P, X) - pt ).squaredNorm(), max_error );
    }
    return max_error;
  }
};

/// Solve the translations and the structure of a non-central-camera triplet which have known rotations
struct translations_Triplet_Solver
{
  enum { MINIMUM_SAMPLES = 4 };
  enum { MAX_MODELS = 1 };
  // Solve the computation of translations & structure of a triplet of non central cameras from image observations.
  static void Solve
  (
    const std::vector< std::vector < std::vector <double> > > pt,
    const std::vector<Mat3> & vec_local_rotations, // subcameras rotations (in the non central camera referential frame)
    const std::vector<Vec3> & vec_local_centers,   // subcameras center    (in the non central camera referential frame)
    const std::vector<Mat3> & vec_global_rotations,// global pose rotations of the rig (the non central camera)
    std::vector<non_central_camera::PoseTripletErrorModel> *P,
    const double ThresholdUpperBound
  )
  {
    //Build the megaMatMatrix
    size_t n_obs = 0;
    for (size_t i=0; i < pt.size() ; ++i)
      n_obs += pt[i].size();

    Mat megaMat(5, n_obs);
    {
      size_t cpt = 0;
      for (size_t i = 0; i  < pt.size(); ++i)
      {
        for(size_t j = 0; j < pt[i].size(); ++j)
        {
          // feature x,y, 3d point idx, intrinsic idx and pose idx
          megaMat.col(cpt++) << pt[i][j][0], pt[i][j][1], i, pt[i][j][2], pt[i][j][3];
        }
      }
    }

    //-- Solve the LInfinity translation and structure from Rotation and points data.
    std::vector<double> vec_solution((3 + MINIMUM_SAMPLES)*3);

    using namespace openMVG::lInfinityCV;

#ifdef OPENMVG_HAVE_MOSEK
    MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
#else
    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
#endif

    Rig_Translation_Structure_L1_ConstraintBuilder cstBuilder(vec_global_rotations, megaMat, vec_local_rotations, vec_local_centers);
    double gamma;
    if (BisectionLP<Rig_Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      LPsolver,
      cstBuilder,
      &vec_solution,
      ThresholdUpperBound, // admissible residual
      0.0, 1e-8, 2, &gamma, false))
    {
      non_central_camera::PoseTripletErrorModel PTemp;
      PTemp.R1 = vec_global_rotations[0]; PTemp.t1 = Vec3(vec_solution[0], vec_solution[1], vec_solution[2]);
      PTemp.R2 = vec_global_rotations[1]; PTemp.t2 = Vec3(vec_solution[3], vec_solution[4], vec_solution[5]);
      PTemp.R3 = vec_global_rotations[2]; PTemp.t3 = Vec3(vec_solution[6], vec_solution[7], vec_solution[8]);

      P->push_back(PTemp);
    }
  }

  // Compute the residual of reprojections
  static double Error
  (
    const non_central_camera::PoseTripletErrorModel & Tensor,
    const std::vector < std::vector <double> > & featInfo,
    const std::vector<Mat3> & vec_local_rotations,
    const std::vector<Vec3> & vec_local_centers
  )
  {
    return non_central_camera::PoseTripletErrorModel::Error(
      Tensor, featInfo, vec_local_rotations, vec_local_centers);
  }
};

}  // namespace non_central_camera
}  // namespace openMVG

