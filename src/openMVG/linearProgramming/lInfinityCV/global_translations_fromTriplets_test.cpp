// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTij.hpp"

#include "openMVG/multiview/translation_averaging_test.hpp"
#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::linearProgramming;
using namespace lInfinityCV;
using namespace std;

TEST(translation_averaging, globalTi_from_tijs_Triplets) {

  const int focal = 1000;
  const int principal_Point = 500;
  //-- Setup a circular camera rig or "cardiod".
  const int iNviews = 12;
  const int iNbPoints = 6;

  const bool bCardiod = true;
  const bool bRelative_Translation_PerTriplet = true;
  std::vector<openMVG::RelativeInfo_Vec > vec_relative_estimates;

  const NViewDataSet d =
    Setup_RelativeTranslations_AndNviewDataset
    (
      vec_relative_estimates,
      focal, principal_Point, iNviews, iNbPoints,
      bCardiod, bRelative_Translation_PerTriplet
    );

  d.ExportToPLY("global_translations_from_triplets_GT.ply");
  visibleCamPosToSVGSurface(d._C, "global_translations_from_triplets_GT.svg");

  //-- Compute the global translations from the triplets of heading directions
  //-   with the L_infinity optimization

  std::vector<double> vec_solution(iNviews*3 + vec_relative_estimates.size() + 1);
  double gamma = -1.0;

  //- a. Setup the LP solver,
  //- b. Setup the constraints generator (for the dedicated L_inf problem),
  //- c. Build constraints and solve the problem,
  //- d. Get back the estimated parameters.

  //- a. Setup the LP solver,
  OSI_CLP_SolverWrapper solverLP(vec_solution.size());

  //- b. Setup the constraints generator (for the dedicated L_inf problem),
  Tifromtij_ConstraintBuilder cstBuilder(vec_relative_estimates);

  //- c. Build constraints and solve the problem (Setup constraints and solver)
  LP_Constraints_Sparse constraint;
  cstBuilder.Build(constraint);
  solverLP.setup(constraint);
  //-- Solving
  EXPECT_TRUE(solverLP.solve()); // the linear program must have a solution

  //- d. Get back the estimated parameters.
  solverLP.getSolution(vec_solution);
  gamma = vec_solution[vec_solution.size()-1];

  //--
  //-- Unit test checking about the found solution
  //--
  EXPECT_NEAR(0.0, gamma, 1e-6); // Gamma must be 0, no noise, perfect data have been sent

  std::cout << "Found solution with gamma = " << gamma << std::endl;

  //-- Get back computed camera translations
  std::vector<double> vec_camTranslation(iNviews*3,0);
  std::copy(&vec_solution[0], &vec_solution[iNviews*3], &vec_camTranslation[0]);

  //-- Get back computed lambda factors
  std::vector<double> vec_camRelLambdas(&vec_solution[iNviews*3], &vec_solution[iNviews*3 + vec_relative_estimates.size()]);
  // lambda factors must be equal to 1.0 (no compression, no dilation);
  EXPECT_NEAR(vec_relative_estimates.size(), std::accumulate (vec_camRelLambdas.begin(), vec_camRelLambdas.end(), 0.0), 1e-6);

  // Get back the camera translations in the global frame:
  std::cout << std::endl << "Camera centers (Computed): " << std::endl;
  for (size_t i = 0; i < iNviews; ++i)
  {
    const Vec3 C_GT = d._C[i] - d._C[0]; //First camera supposed to be at Identity

    const Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
    const Mat3 & Ri = d._R[i];
    const Vec3 C_computed = - Ri.transpose() * t;

    //-- Check that found camera position is equal to GT value
    EXPECT_NEAR(0.0, DistanceLInfinity(C_computed, C_GT), 1e-6);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
