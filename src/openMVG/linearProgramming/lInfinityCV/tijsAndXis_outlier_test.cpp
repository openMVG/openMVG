// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"
#include "testing/testing.h"

#include "openMVG/multiview/projection.hpp"

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"

#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri_noise.hpp"

#include <iostream>
#include <vector>

using namespace openMVG;

using namespace linearProgramming;
using namespace lInfinityCV;

TEST(Translation_Structure_L_Infinity_Noisy, Outlier_OSICLP_SOLVER) {

  const int nViews = 5;
  const int nbPoints = 5;
  const double focalValue = 1000;
  const double cx = 500, cy = 500;

  const NViewDataSet d =
    //NRealisticCamerasRing(nViews, nbPoints,
    NRealisticCamerasCardioid(nViews, nbPoints,
    nViewDatasetConfigurator(focalValue,focalValue,cx,cy,5,0));

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NViewDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  {
    // Prepare Rotation matrix (premultiplied by K)
    // Apply the K matrix to the rotation
    Mat3 K;
    K << focalValue, 0, cx,
      0, focalValue, cy,
      0, 0, 1;
    std::vector<Mat3> vec_KRotation(nViews);
    for (size_t i = 0; i < nViews; ++i) {
      vec_KRotation[i] = K * d._R[i];
    }

    //set two point as outlier
    d2._x[0].col(0)(0) +=10; //Camera 0, measurement 0, noise on X coord
    d2._x[3].col(3)(1) -=8;  //Camera 3, measurement 3, noise on Y coord

    //Create the mega matrix
    Mat megaMat(4, d._n*d._x[0].cols());
    {
      int cpt = 0;
      for (int i=0; i<d._n;++i)
      {
        const int camIndex = i;
        for (int j=0; j<d._x[0].cols(); ++j)
        {
          megaMat(0,cpt) = d2._x[camIndex].col(j)(0);
          megaMat(1,cpt) = d2._x[camIndex].col(j)(1);
          megaMat(2,cpt) = j;
          megaMat(3,cpt) = camIndex;
          cpt++;
        }
      }
    }

    double admissibleResidual = 1.0/focalValue;
    std::vector<double> vec_solution((nViews + nbPoints + megaMat.cols())*3);
    OSI_CLP_SolverWrapper wrapperOSICLPSolver(vec_solution.size());
    TiXi_withNoise_L1_ConstraintBuilder cstBuilder( vec_KRotation, megaMat);
    std::cout << std::endl << "Bisection returns : " << std::endl;
    std::cout<< BisectionLP<TiXi_withNoise_L1_ConstraintBuilder, LP_Constraints_Sparse>(
            wrapperOSICLPSolver,
            cstBuilder,
            &vec_solution,
            admissibleResidual,
            0.0, 1e-8);

    std::cout << "Found solution:\n";
    std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

    //-- First the ti and after the Xi :

    //-- Fill the ti
    for (int i=0; i < nViews; ++i)  {
      d2._t[i] = K.inverse() * Vec3(vec_solution[3*i], vec_solution[3*i+1], vec_solution[3*i+2]);
      // Change Ci to -Ri*Ci
      d2._C[i] = -d2._R[i].inverse() * d2._t[i];
    }

    for (int i=0; i < nbPoints; ++i)  {
      size_t index = 3*nViews;
      d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
    }

    // Compute residuals L2 from estimated parameter values :
    std::cout << std::endl << "Residual : " << std::endl;
    Vec2 xk;
    double xsum = 0.0;
    for (int i = 0; i < d2._n; ++i) {
        std::cout << "\nCamera : " << i << " \t:";
        for(int k = 0; k < d._x[0].cols(); ++k)
        {
          xk = Project(d2.P(i),  Vec3(d2._X.col(k)));
          double residual = (xk - d2._x[i].col(k)).norm();
          std::cout << Vec2(( xk - d2._x[i].col(k)).array().pow(2)).array().sqrt().mean() <<"\t";
          //-- Check that were measurement are not noisy the residual is small
          //  check were the measurement are noisy, the residual is important
          //if ((i != 0 && k != 0) || (i!=3 && k !=3))
          if ((i != 0 && k != 0) && (i!=3 && k !=3)) {
            EXPECT_NEAR(0.0, residual, 1e-5);
            xsum += residual;
          }
        }
        std::cout << std::endl;
    }
    double dResidual = xsum / (d2._n*d._x[0].cols());
    std::cout << std::endl << "Residual mean in not noisy measurement: " << dResidual << std::endl;
    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual, 1e-1);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}

#ifdef OPENMVG_HAVE_MOSEK
TEST(Translation_Structure_L_Infinity_Noisy, Outlier_MOSEK) {

  const int nViews = 5;
  const int nbPoints = 5;
  const double focalValue = 1000;
  const double cx = 500, cy = 500;

  const NViewDataSet d =
    //NRealisticCamerasRing(nViews, nbPoints,
    NRealisticCamerasCardioid(nViews, nbPoints,
    nViewDatasetConfigurator(focalValue,focalValue,cx,cy,5,0));

  d.ExportToPLY("test_Before_Infinity.ply");
  //-- Test triangulation of all the point
  NViewDataSet d2 = d;

  //-- Set to 0 the future computed data to be sure of computation results :
  d2._X.fill(0); //Set _Xi of dataset 2 to 0 to be sure of new data computation
  fill(d2._t.begin(), d2._t.end(), Vec3(0.0,0.0,0.0));

  {
    // Prepare Rotation matrix (premultiplied by K)
    // Apply the K matrix to the rotation
    Mat3 K;
    K << focalValue, 0, cx,
      0, focalValue, cy,
      0, 0, 1;
    std::vector<Mat3> vec_KRotation(nViews);
    for (size_t i = 0; i < nViews; ++i) {
      vec_KRotation[i] = K * d._R[i];
    }

    //set two point as outlier
    d2._x[0].col(0)(0) +=10; //Camera 0, measurement 0, noise on X coord
    d2._x[3].col(3)(1) -=8;  //Camera 3, measurement 3, noise on Y coord

    //Create the mega matrix
    Mat megaMat(4, d._n*d._x[0].cols());
    {
      int cpt = 0;
      for (int i=0; i<d._n;++i)
      {
        const int camIndex = i;
        for (int j=0; j<d._x[0].cols(); ++j)
        {
          megaMat(0,cpt) = d2._x[camIndex].col(j)(0);
          megaMat(1,cpt) = d2._x[camIndex].col(j)(1);
          megaMat(2,cpt) = j;
          megaMat(3,cpt) = camIndex;
          cpt++;
        }
      }
    }

    double admissibleResidual = 1.0/focalValue;
    std::vector<double> vec_solution((nViews + nbPoints + megaMat.cols())*3);
    MOSEK_SolveWrapper wrapperMosek(vec_solution.size());
    TiXi_withNoise_L1_ConstraintBuilder cstBuilder( vec_KRotation, megaMat);
    std::cout << std::endl << "Bisection returns : " << std::endl;
    std::cout<< BisectionLP<TiXi_withNoise_L1_ConstraintBuilder, LP_Constraints_Sparse>(
            wrapperMosek,
            cstBuilder,
            &vec_solution,
            admissibleResidual,
            0.0, 1e-8);

    std::cout << "Found solution:\n";
    std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

    //-- First the ti and after the Xi :

    //-- Fill the ti
    for (int i=0; i < nViews; ++i)  {
      d2._t[i] = K.inverse() * Vec3(vec_solution[3*i], vec_solution[3*i+1], vec_solution[3*i+2]);
      // Change Ci to -Ri*Ci
      d2._C[i] = -d2._R[i].inverse() * d2._t[i];
    }

    for (int i=0; i < nbPoints; ++i)  {
      size_t index = 3*nViews;
      d2._X.col(i) = Vec3(vec_solution[index+i*3], vec_solution[index+i*3+1], vec_solution[index+i*3+2]);
    }

    // Compute residuals L2 from estimated parameter values :
    std::cout << std::endl << "Residual : " << std::endl;
    Vec2 xk;
    double xsum = 0.0;
    for (int i = 0; i < d2._n; ++i) {
        std::cout << "\nCamera : " << i << " \t:";
        for(int k = 0; k < d._x[0].cols(); ++k)
        {
          xk = Project(d2.P(i),  Vec3(d2._X.col(k)));
          double residual = (xk - d2._x[i].col(k)).norm();
          std::cout << Vec2(( xk - d2._x[i].col(k)).array().pow(2)).array().sqrt().mean() <<"\t";
          //-- Check that were measurement are not noisy the residual is small
          //  check were the measurement are noisy, the residual is important
          //if ((i != 0 && k != 0) || (i!=3 && k !=3))
          if ((i != 0 && k != 0) && (i!=3 && k !=3)) {
            EXPECT_NEAR(0.0, residual, 1e-6);
            xsum += residual;
          }
        }
        std::cout << std::endl;
    }
    double dResidual = xsum / (d2._n*d._x[0].cols());
    std::cout << std::endl << "Residual mean in not noisy measurement: " << dResidual << std::endl;
    // Check that 2D re-projection and 3D point are near to GT.
    EXPECT_NEAR(0.0, dResidual, 1e-1);
  }

  d2.ExportToPLY("test_After_Infinity.ply");
}
#endif // OPENMVG_HAVE_MOSEK

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
