// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"


#include <iostream>
#include <vector>

using namespace openMVG;

using namespace linearProgramming;
using namespace lInfinityCV;

TEST(Translation_Structure_L_Infinity, OSICLP_SOLVER) {

  const size_t nViews = 3;
  const size_t nbPoints = 6;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NViewDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  Mat megaMat(4, d._n*d._x[0].cols());
  {
    size_t cpt = 0;
    for (size_t i=0; i<d._n;++i)
    {
      const size_t camIndex = i;
      for (Mat2X::Index j=0; j < d._x[0].cols(); ++j)
      {
        megaMat(0,cpt) = d._x[camIndex].col(j)(0);
        megaMat(1,cpt) = d._x[camIndex].col(j)(1);
        megaMat(2,cpt) = j;
        megaMat(3,cpt) = camIndex;
        cpt++;
      }
    }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nViews + nbPoints)*3);

    OSI_CLP_SolverWrapper wrapperOSICLPSolver(vec_solution.size());
    Translation_Structure_L1_ConstraintBuilder cstBuilder( d._R, megaMat);
    EXPECT_TRUE(
      (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      wrapperOSICLPSolver,
      cstBuilder,
      &vec_solution,
      1.0,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nViews; ++i)
      {
        size_t index = i*3;
        d2._t[i] = Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i] * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nViews;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i < d2._n; ++i) {
      for (Mat2X::Index k = 0; k < d._x[0].cols(); ++k)
      {
        xk = Project(d2.P(i), Vec3(d2._X.col(k)));
        xsum += Vec2(( xk - d2._x[i].col(k)).array().pow(2));
      }
    }
    const double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}

TEST(Translation_Structure_L_Infinity, OSICLP_SOLVER_K) {

  const size_t nViews = 3;
  const size_t nbPoints = 6;
  const NViewDataSet d = NRealisticCamerasRing(nViews, nbPoints,
    nViewDatasetConfigurator(1000,1000,500,500,5,0)); // Suppose a camera with Unit matrix as K

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NViewDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  //Create the mega matrix
  Mat megaMat(4, d._n*d._x[0].cols());
  {
    size_t cpt = 0;
    for (size_t i=0; i < d._n;++i)
    {
      const size_t camIndex = i;
      for (size_t j=0; j < (size_t)d._x[0].cols(); ++j)
      {
        megaMat(0,cpt) = d._x[camIndex].col(j)(0);
        megaMat(1,cpt) = d._x[camIndex].col(j)(1);
        megaMat(2,cpt) = j;
        megaMat(3,cpt) = camIndex;
        cpt++;
      }
    }
  }

  // Solve the problem and check that fitted value are good enough
  {
    std::vector<double> vec_solution((nViews + nbPoints)*3);

    std::vector<Mat3> vec_KR(d._R);
    for (size_t i=0;i < nViews; ++i)
      vec_KR[i] = d._K[0] * d._R[i];

    OSI_CLP_SolverWrapper wrapperOSICLPSolver(vec_solution.size());
    Translation_Structure_L1_ConstraintBuilder cstBuilder( vec_KR, megaMat);
    EXPECT_TRUE(
      (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
      wrapperOSICLPSolver,
      cstBuilder,
      &vec_solution,
      1.0,
      0.0))
    );

    // Move computed value to dataset for residual estimation.
    {
      //-- Fill the ti
      for (size_t i=0; i < nViews; ++i)
      {
        size_t index = i*3;
        d2._t[i] = d._K[0].inverse() * Vec3(vec_solution[index], vec_solution[index+1], vec_solution[index+2]);
        // Change Ci to -Ri*Ci
        d2._C[i] = -d2._R[i] * d2._t[i];
      }

      //-- Now the Xi :
      for (size_t i=0; i < nbPoints; ++i) {
        size_t index = 3*nViews;
        d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
      }
    }

    // Compute residuals L2 from estimated parameter values :
    Vec2 xk, xsum(0.0,0.0);
    for (size_t i = 0; i < d2._n; ++i) {
      for (size_t k = 0; k < (size_t)d._x[0].cols(); ++k)
      {
        xk = Project(d2.P(i), Vec3(d2._X.col(k)));
        xsum += Vec2(( xk - d2._x[i].col(k)).array().pow(2));
      }
    }
    double dResidual2D = (xsum.array().sqrt().sum());

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-4);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
