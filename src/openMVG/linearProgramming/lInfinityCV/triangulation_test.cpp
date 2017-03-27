// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/lInfinityCV/triangulation.hpp"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/test_data_sets.hpp"

#include "testing/testing.h"

#include <iostream>
#include <vector>

using namespace openMVG;
using namespace linearProgramming;
using namespace lInfinityCV;

TEST(lInfinityCV, Triangulation_OSICLPSOLVER) {

  const NViewDataSet d = NRealisticCamerasRing(6, 10,
    nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  std::vector<Mat34> vec_Pi;

  d.ExportToPLY("test_Before_Infinity_Triangulation_OSICLP.ply");
  //-- Test triangulation of all the point
  NViewDataSet d2 = d;
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation

  for (size_t i = 0; i < d._n; ++i)
    vec_Pi.push_back(d.P(i));

  for (Mat2X::Index k = 0; k < d._x[0].cols(); ++k)
  {
    Mat2X x_ij;
    x_ij.resize(2,d._n);
    for (size_t i = 0; i < d._n; ++i)
      x_ij.col(i) = d._x[i].col(k);

    std::vector<double> vec_solution(3);

    OSI_CLP_SolverWrapper wrapperOSICLPSolver(3);
    Triangulation_L1_ConstraintBuilder cstBuilder(vec_Pi, x_ij);
    // Use bisection in order to find the global optimum and so find the
    //  best triangulated point under the L_infinity norm
    EXPECT_TRUE(
      (BisectionLP<Triangulation_L1_ConstraintBuilder,LP_Constraints>(
      wrapperOSICLPSolver,
      cstBuilder,
      &vec_solution,
      1.0,
      0.0))
    );

    Vec3 XSolution(vec_solution[0], vec_solution[1], vec_solution[2]);
    d2._X.col(k) = XSolution.transpose();

    // Compute residuals L2 from estimated parameter values :
    const Vec3 & X = XSolution;
    Vec2 x1, xsum(0.0,0.0);
    for (size_t i = 0; i < d2._n; ++i) {
      x1 = Project(d2.P(i), X);
      xsum += Vec2((x1-x_ij.col(i)).array().pow(2));
    }
    const double dResidual2D = (xsum.array().sqrt().sum());

    // Residual LInfinity between GT 3D point and found one
    const double dResidual3D = DistanceLInfinity(XSolution, Vec3(d._X.col(k)));

    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual2D, 1e-5);
    EXPECT_NEAR(0.0, dResidual3D, 1e-5);
  }
  d2.ExportToPLY("test_After_Infinity_Triangulation_OSICLP.ply");
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
