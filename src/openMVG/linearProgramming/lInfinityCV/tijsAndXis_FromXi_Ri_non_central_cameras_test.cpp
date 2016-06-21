// Copyright (c) 2016 Pierre MOULON, Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/test_data_sets.hpp"
#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"

#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_FromXi_Ri_non_central_cameras.hpp"

#include <iostream>
#include <vector>

using namespace openMVG;

using namespace linearProgramming;
using namespace lInfinityCV;

/*
TEST(Translation_Structure_L_Infinity, OSICLP_SOLVER) {

  const size_t nPoses = 3;
  const size_t rig_size = 2;
  const size_t nViews = rig_size * nPoses;
  const size_t nbPoints = 4;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity_translation.ply");
  //-- Test triangulation of all the point
  NViewDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  //Build the megaMatMatrix
  int n_obs = rig_size * nPoses * nbPoints;

  Mat megaMat(5, n_obs);
  {
    size_t cpt = 0;
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t k = 0 ; k < rig_size; ++k )
        for(size_t j = 0; j < nPoses; ++j)
        {
          megaMat(0, cpt) = d2._x[j * rig_size + k].col(i)(0); // feature x
          megaMat(1, cpt) = d2._x[j * rig_size + k].col(i)(1); // feature y
          megaMat(2, cpt) = i;    // point index
          megaMat(3, cpt) = k;    // sub pose index
          megaMat(4, cpt) = j;    // pose index
          ++cpt;
        }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nPoses + nbPoints)*3);
    using namespace openMVG::lInfinityCV;

    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));

    Rig_Translation_Structure_L1_ConstraintBuilder cstBuilder(d2._R, megaMat, d2._rotations, d2._offsets);
    EXPECT_TRUE ( (BisectionLP<Rig_Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      LPsolver,
      cstBuilder,
      &vec_solution,
      1.0,//admissibleResidual,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nPoses; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i].transpose() * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nPoses;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t j = 0; j < nPoses; ++j)
        for(size_t k = 0 ; k < rig_size; ++k )
        {
          xk = Project(d2.P(j,k), Vec3(d2._X.col(i)));
          xsum += Vec2(( xk - d2._x[j * rig_size + k].col(i)).array().pow(2));
        }

    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity_translation.ply");
}
*/
/*
TEST(Translation_Structure_L_Infinity, OSICLP_SOLVER_K) {

  const size_t nPoses = 3;
  const size_t rig_size = 2;
  const size_t nViews = rig_size * nPoses ;
  const size_t nbPoints = 4;
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
    nPoseDatasetConfigurator(1000,1000,500,500,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity_translation_K.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  // check that the scene is well constructed
  double  triangulation_error = 0.0;
  for (size_t i = 0; i  < nbPoints; ++i)
  {
    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    std::vector < std::pair <Mat34, Vec2> >  views;

    for(size_t j = 0; j < nPoses; ++j)
      for(size_t k = 0 ; k < rig_size; ++k )
      {
        //compute projection matrices
        Vec2 pt;
        pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1);
        views.push_back(std::make_pair ( d.P(j,k), pt ) );
      }

    // update triangulation object
    for( size_t i = 0 ; i < views.size(); ++i )
      triangulationObj.add ( views[i].first, views[i].second );

    const Vec3 X = triangulationObj.compute();
    triangulation_error += (X -d._X.col(i)).norm();
  }

  // Check that triangulation point are near to the inital data
  EXPECT_NEAR(0.0, triangulation_error, 1e-10);

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  //Build the megaMatMatrix
  int n_obs = rig_size * nPoses * nbPoints;

  Mat megaMat(5, n_obs);
  {
    size_t cpt = 0;
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t k = 0 ; k < rig_size; ++k )
        for(size_t j = 0; j < nPoses; ++j)
        {
          //compute projection matrices
          Vec3 pt;
          pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1), 1.0;

          pt = d2._K[k].inverse() * pt ;

          megaMat(0, cpt) = pt(0) / pt(2); // feature x
          megaMat(1, cpt) = pt(1) / pt(2); // feature y
          megaMat(2, cpt) = i;    // point index
          megaMat(3, cpt) = k;    // sub pose index
          megaMat(4, cpt) = j;    // pose index
          ++cpt;
        }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nPoses + nbPoints)*3);
    using namespace openMVG::lInfinityCV;

#ifdef OPENMVG_HAVE_MOSEK
    MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
#else
    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
#endif

    Rig_Translation_Structure_L1_ConstraintBuilder cstBuilder(d2._R, megaMat, d2._rotations, d2._offsets);
    EXPECT_TRUE ( (BisectionLP<Rig_Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      LPsolver,
      cstBuilder,
      &vec_solution,
      1.0,//admissibleResidual,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nPoses; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i].transpose() * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nPoses;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t j = 0; j < nPoses; ++j)
        for(size_t k = 0 ; k < rig_size; ++k )
        {
          xk = Project(d2.P(j,k), Vec3(d2._X.col(i)));
          xsum += Vec2(( xk - d2._x[j * rig_size + k].col(i)).array().pow(2));
        }

    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4 * 1000 ); // update tolerance to be in pixels
  }

  d2.ExportToPLY("test_After_Infinity_translation_K.ply");
}

TEST(Center_Structure_L_Infinity, OSICLP_SOLVER) {

  const size_t nPoses = 3;
  const size_t rig_size = 2;
  const size_t nViews = rig_size * nPoses;
  const size_t nbPoints = 4;
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
    nPoseDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity_center.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  // check that the scene is well constructed
  double  triangulation_error = 0.0;
  for (size_t i = 0; i  < nbPoints; ++i)
  {
    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    std::vector < std::pair <Mat34, Vec2> >  views;

    for(size_t j = 0; j < nPoses; ++j)
      for(size_t k = 0 ; k < rig_size; ++k )
      {
        //compute projection matrices
        Vec2 pt;
        pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1);
        views.push_back(std::make_pair ( d.P(j,k), pt ) );
      }

    // update triangulation object
    for( size_t i = 0 ; i < views.size(); ++i )
      triangulationObj.add ( views[i].first, views[i].second );

    const Vec3 X = triangulationObj.compute();
    triangulation_error += (X -d._X.col(i)).norm();
  }

  // Check that triangulation point are near to the inital data
  EXPECT_NEAR(0.0, triangulation_error, 1e-10);

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  //Build the megaMatMatrix
  int n_obs = rig_size * nPoses * nbPoints;

  Mat megaMat(5, n_obs);
  {
    size_t cpt = 0;
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t k = 0 ; k < rig_size; ++k )
        for(size_t j = 0; j < nPoses; ++j)
        {
          megaMat(0, cpt) = d2._x[j * rig_size + k].col(i)(0); // feature x
          megaMat(1, cpt) = d2._x[j * rig_size + k].col(i)(1); // feature y
          megaMat(2, cpt) = i;    // point index
          megaMat(3, cpt) = k;    // sub pose index
          megaMat(4, cpt) = j;    // pose index
          ++cpt;
        }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution( (nPoses*3 + n_obs) );
    using namespace openMVG::lInfinityCV;

#ifdef OPENMVG_HAVE_MOSEK
    MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
#else
    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
#endif

    Rig_Center_Structure_L1_ConstraintBuilder cstBuilder(d2._R, megaMat, d2._rotations, d2._offsets);
    EXPECT_TRUE ( (BisectionLP<Rig_Center_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      LPsolver,
      cstBuilder,
      &vec_solution,
      1.0,//admissibleResidual,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nPoses; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i].transpose() * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nPoses;
        Vec3  bearing_vector ;
        bearing_vector << d2._x[0].col(i)(0), d2._x[0].col(i)(1), 1.0;
        const Vec3  X_from_depth = ( vec_solution[index + i*nPoses*rig_size]
                                     * d2._R[0].transpose() * d2._rotations[0].transpose() * bearing_vector )
                                   + d2._R[0].transpose() * d2._offsets[0] + d2._C[0];
        d2._X.col(i) = X_from_depth;
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t j = 0; j < nPoses; ++j)
        for(size_t k = 0 ; k < rig_size; ++k )
        {
          xk = Project(d2.P(j,k), Vec3(d2._X.col(i)));
          xsum += Vec2(( xk - d2._x[j * rig_size + k].col(i)).array().pow(2));
        }

    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity_center.ply");
}

TEST(Center_Structure_L_Infinity, OSICLP_SOLVER_K) {

  const size_t nPoses = 3;
  const size_t rig_size = 2;
  const size_t nViews = rig_size * nPoses ;
  const size_t nbPoints = 4;
  const NPoseDataSet d = NRealisticPosesRing(nPoses, nViews, rig_size, nbPoints,
    nPoseDatasetConfigurator(1000,1000,500,500,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity_center_K.ply");
  //-- Test triangulation of all the point
  NPoseDataSet d2 = d;

  // check that the scene is well constructed
  double  triangulation_error = 0.0;
  for (size_t i = 0; i  < nbPoints; ++i)
  {
    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    std::vector < std::pair <Mat34, Vec2> >  views;

    for(size_t j = 0; j < nPoses; ++j)
      for(size_t k = 0 ; k < rig_size; ++k )
      {
        //compute projection matrices
        Vec2 pt;
        pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1);
        views.push_back(std::make_pair ( d.P(j,k), pt ) );
      }

    // update triangulation object
    for( size_t i = 0 ; i < views.size(); ++i )
      triangulationObj.add ( views[i].first, views[i].second );

    const Vec3 X = triangulationObj.compute();
    triangulation_error += (X -d._X.col(i)).norm();
  }

  // Check that triangulation point are near to the inital data
  EXPECT_NEAR(0.0, triangulation_error, 1e-10);

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  //Build the megaMatMatrix
  int n_obs = rig_size * nPoses * nbPoints;

  Mat megaMat(5, n_obs);
  {
    size_t cpt = 0;
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t k = 0 ; k < rig_size; ++k )
        for(size_t j = 0; j < nPoses; ++j)
        {
          //compute projection matrices
          Vec3 pt;
          pt << d2._x[j * rig_size + k].col(i)(0), d2._x[j * rig_size + k].col(i)(1), 1.0;

          pt = d2._K[k].inverse() * pt ;

          megaMat(0, cpt) = pt(0) / pt(2); // feature x
          megaMat(1, cpt) = pt(1) / pt(2); // feature y
          megaMat(2, cpt) = i;    // point index
          megaMat(3, cpt) = k;    // sub pose index
          megaMat(4, cpt) = j;    // pose index
          ++cpt;
        }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution( (nPoses*3 + n_obs) );
    using namespace openMVG::lInfinityCV;

#ifdef OPENMVG_HAVE_MOSEK
    MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
#else
    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
#endif

    Rig_Center_Structure_L1_ConstraintBuilder cstBuilder(d2._R, megaMat, d2._rotations, d2._offsets);
    EXPECT_TRUE ( (BisectionLP<Rig_Center_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      LPsolver,
      cstBuilder,
      &vec_solution,
      1.0,//admissibleResidual,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nPoses; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i].transpose() * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nPoses;
        Vec3  bearing_vector ;
        bearing_vector << d2._x[0].col(i)(0), d2._x[0].col(i)(1), 1.0;
        const Vec3  X_from_depth = ( vec_solution[index + i*nPoses*rig_size]
                                    * d2._R[0].transpose() * d2._rotations[0].transpose()
                                    * d2._K[0].inverse() * bearing_vector )
                                    + d2._R[0].transpose() * d2._offsets[0] + d2._C[0];
        d2._X.col(i) = X_from_depth;
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i  < nbPoints; ++i)
      for(size_t j = 0; j < nPoses; ++j)
        for(size_t k = 0 ; k < 1; ++k )
        {
          xk = Project(d2.P(j,k), Vec3(d2._X.col(i)));
          xsum += Vec2(( xk - d2._x[j * rig_size + k].col(i)).array().pow(2));
        }

    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4 * 1000 ); // update tolerance to be in pixels
  }

  d2.ExportToPLY("test_After_Infinity_center_K.ply");
}
*/

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
