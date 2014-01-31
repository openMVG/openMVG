
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTriplets.hpp"

#include "openMVG/multiview/essential.hpp"

#include "openMVG/graph/triplet_finder.hpp"
using namespace openMVG::graphUtils;

#include "third_party/vectorGraphics/svgDrawer.hpp"
using namespace svg;

#include "openMVG/multiview/test_data_sets.hpp"
#include "testing/testing.h"

#include <fstream>
#include <map>
#include <utility>
#include <vector>

using namespace openMVG;
using namespace openMVG::linearProgramming;
using namespace lInfinityCV;
using namespace std;

int modifiedMod(int number, int modulus)
{
   int result = number % modulus;
   if (result < 0) result += modulus;
   return result;
}

//-- Export a series of camera positions to a SVG surface of specified squared size
void visibleCamPosToSVGSurface(
  const std::vector<Vec3> & vec_Ci,
  const std::string & fileName)
{
  Mat points(3, vec_Ci.size());
  for(size_t i = 0; i  < vec_Ci.size(); ++i)
  {
    points.col(i) = vec_Ci[i];
  }

  Vec mean, variance;
  MeanAndVarianceAlongRows(points, &mean, &variance);

  double xfactor = sqrt(2.0 / variance(0));
  double yfactor = sqrt(2.0 / variance(2));

  std::vector<Vec3> out = vec_Ci;
  for(size_t i = 0; i  < vec_Ci.size(); ++i)
  {
    out[i](0) = ((out[i](0) * xfactor) + -xfactor * mean(0)) * 30 + 100;
    out[i](2) = ((out[i](2) * yfactor) + -yfactor * mean(2)) * 30 + 100;
  }

  if (!fileName.empty())
  {
    const double size = 200;
    svgDrawer svgSurface_GT(size,size);
    for(size_t i = 0; i  < vec_Ci.size(); ++i)
    {
      svgSurface_GT.drawCircle(out[i](0), out[i](2),
                               3,svgStyle().stroke("black",0.2).fill("red"));
    }
    std::ostringstream osSvgGT;
    osSvgGT << fileName;
    std::ofstream svgFileGT( osSvgGT.str().c_str());
    svgFileGT << svgSurface_GT.closeSvgFile().str();
  }
}

TEST(translation_averaging, globalTi_from_tijs_Triplets) {

  int focal = 1000;
  int principal_Point = 500;
  //-- Setup a circular camera rig or "cardiod".
  const int iNviews = 12;
  const int iNbPoints = 6;
  NViewDataSet d =
    //NRealisticCamerasRing(
      NRealisticCamerasCardioid(
      iNviews, iNbPoints,
      nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0));

  d.ExportToPLY("global_translations_from_triplets_GT.ply");

  visibleCamPosToSVGSurface(d._C, "global_translations_from_triplets_GT.svg");

  // List sucessives triplets of the large loop of camera
  std::vector< graphUtils::Triplet > vec_triplets;
  for (size_t i = 0; i < iNviews; ++i)
  {
    const size_t iPlus1 = modifiedMod(i+1,iNviews);
    const size_t iPlus2 = modifiedMod(i+2,iNviews);
    //-- sort the triplet index to have a monotic ascending series of value
    size_t triplet[3] = {i, iPlus1, iPlus2};
    std::sort(&triplet[0], &triplet[3]);
    vec_triplets.push_back(Triplet(triplet[0],triplet[1],triplet[2]));
  }

  //- For each triplet compute relative translations and rotations motions
  std::vector<openMVG::lInfinityCV::relativeInfo > vec_initialEstimates;

  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    size_t I = triplet.i, J = triplet.j , K = triplet.k;

    //-- Build camera alias over GT translations and rotations:
    const Mat3 & RI = d._R[I];
    const Mat3 & RJ = d._R[J];
    const Mat3 & RK = d._R[K];
    const Vec3 & ti = d._t[I];
    const Vec3 & tj = d._t[J];
    const Vec3 & tk = d._t[K];

    //-- Build relatives motions (that feeds the Linear program formulation)
    {
      Mat3 RijGt;
      Vec3 tij;
      RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
      vec_initialEstimates.push_back(
        std::make_pair(std::make_pair(I, J), std::make_pair(RijGt, tij)));

      Mat3 RjkGt;
      Vec3 tjk;
      RelativeCameraMotion(RJ, tj, RK, tk, &RjkGt, &tjk);
      vec_initialEstimates.push_back(
        std::make_pair(std::make_pair(J, K), std::make_pair(RjkGt, tjk)));

      Mat3 RikGt;
      Vec3 tik;
      RelativeCameraMotion(RI, ti, RK, tk, &RikGt, &tik);
      vec_initialEstimates.push_back(
        std::make_pair(std::make_pair(I, K), std::make_pair(RikGt, tik)));
    }
  }

  //-- Compute the global translations from the triplets of heading directions
  //-   with the L_infinity optimization

  std::vector<double> vec_solution(iNviews*3 + vec_initialEstimates.size()/3 + 1);
  double gamma = -1.0;

  //- a. Setup the LP solver,
  //- b. Setup the constraints generator (for the dedicated L_inf problem),
  //- c. Build constraints and solve the problem,
  //- d. Get back the estimated parameters.

  //- a. Setup the LP solver,
  OSI_CLP_SolverWrapper solverLP(vec_solution.size());

  //- b. Setup the constraints generator (for the dedicated L_inf problem),
  Tifromtij_ConstraintBuilder_OneLambdaPerTrif cstBuilder(vec_initialEstimates);

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
  std::vector<double> vec_camRelLambdas(&vec_solution[iNviews*3], &vec_solution[iNviews*3 + vec_initialEstimates.size()/3]);
  // lambda factors must be equal to 1.0 (no compression, no dilation);
  EXPECT_NEAR(vec_initialEstimates.size()/3, std::accumulate (vec_camRelLambdas.begin(), vec_camRelLambdas.end(), 0.0), 1e-6);

  // True camera C
  std::cout << std::endl << "Camera centers (Ground truth): " << std::endl;
  for (size_t i = 0; i < iNviews; ++i)
  {
    std::cout << i << ": " << d._C[i].transpose() - d._C[0].transpose() << std::endl;
  }

  // Get back the camera translations in the global frame:
  std::cout << std::endl << "Camera centers (Computed): " << std::endl;
  for (size_t i = 0; i < iNviews; ++i)
  {
    const Vec3 C_GT = d._C[i].transpose() - d._C[0].transpose(); //First camera supposed to be at Identity

    Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
    const Mat3 & Ri = d._R[i];
    const Vec3 C_computed = - Ri.transpose() * t;

    //-- Check that found camera position is equal to GT value
    EXPECT_NEAR(0.0, DistanceLInfinity(C_computed, C_GT), 1e-6);

    std::cout << i << ": " << C_computed.transpose() << std::endl;
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
