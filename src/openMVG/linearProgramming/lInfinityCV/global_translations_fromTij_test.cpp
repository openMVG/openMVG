
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTij.hpp"

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

TEST(translation_averaging, globalTi_from_tijs) {

  const int focal = 1000;
  const int principal_Point = 500;
  //-- Setup a circular camera rig or "cardiod".
  const int iNviews = 12;
  const int iNbPoints = 6;
  NViewDataSet d =
    //NRealisticCamerasRing(
    NRealisticCamerasCardioid(
        iNviews, iNbPoints,
        nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0));

  d.ExportToPLY("global_translations_from_Tij_GT.ply");

  visibleCamPosToSVGSurface(d._C, "global_translations_from_Tij_GT.svg");

  //- For each relative translations and rotations motions
  std::vector<openMVG::relativeInfo > vec_initialEstimates;

  //-- Setup initial pair that will be considered (link each camera to the two next)
  std::vector< std::pair<size_t,size_t> > map_pairs;
  for (size_t i = 0; i < iNviews; ++i)  {
    for (size_t j=i; j<=i+2; ++j)
    {
      const size_t jj = modifiedMod(j,iNviews);
      if (i != jj)
        map_pairs.push_back(make_pair(i,jj));
    }
  }
  
  for (std::vector< std::pair<size_t,size_t> >::const_iterator
    iter = map_pairs.begin();
    iter != map_pairs.end();
    ++iter)
  {
    const size_t I = (*iter).first;
    const size_t J = (*iter).second;
   
    //-- Build camera alias over GT translations and rotations:
    const Mat3 & RI = d._R[I];
    const Mat3 & RJ = d._R[J];
    const Vec3 & ti = d._t[I];
    const Vec3 & tj = d._t[J];

    //-- Build relative motions (that feeds the Linear program formulation)
    {
      Mat3 RijGt;
      Vec3 tij;
      RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
      //-- normalize tij (keep only direction)
      tij.normalize();
      vec_initialEstimates.push_back(
        std::make_pair(std::make_pair(I, J), std::make_pair(RijGt, tij)));
    }
  }
  
  //-- Compute the global translations from the translation heading directions
  //-   with the L_infinity optimization
  // 3*NCam*[X,Y,Z] ; Ncam*[Lambda], [gamma]
  std::vector<double> vec_solution(iNviews*3 + vec_initialEstimates.size() + 1);

  //- a. Setup the LP solver,
  //- b. Setup the constraints generator (for the dedicated L_inf problem),
  //- c. Build constraints and solve the problem,
  //- d. Get back the estimated parameters.

  //- a. Setup the LP solver,
  OSI_CLP_SolverWrapper solverLP(vec_solution.size());

  //- b. Setup the constraints generator (for the dedicated L_inf problem),
  Tifromtij_ConstraintBuilder cstBuilder(vec_initialEstimates);

  //- c. Build constraints and solve the problem (Setup constraints and solver)
  LP_Constraints_Sparse constraint;
  cstBuilder.Build(constraint);
  solverLP.setup(constraint);
  //-- Solving
  EXPECT_TRUE(solverLP.solve()); // the linear program must have a solution

  //- d. Get back the estimated parameters.
  solverLP.getSolution(vec_solution);
  const double gamma = vec_solution[vec_solution.size()-1];

  //--
  //-- Unit test checking about the found solution
  //--
  EXPECT_NEAR(0.0, gamma, 1e-6); // Gamma must be 0, no noise, perfect data have been sent

  std::cout << "Found solution with gamma = " << gamma << std::endl;

  //-- Get back computed camera translations
  std::vector<double> vec_camTranslation(iNviews*3,0);
  std::copy(&vec_solution[0], &vec_solution[iNviews*3], &vec_camTranslation[0]);

  //-- Get back computed lambda factors
  std::vector<double> vec_camRelLambdas(&vec_solution[iNviews*3], &vec_solution[iNviews*3 + vec_initialEstimates.size()]);

  // Check validity of the camera centers:
  // Check the direction since solution if found up to a scale
  for (size_t i = 0; i < iNviews; ++i)
  {
    const Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
    const Mat3 & Ri = d._R[i];
    const Vec3 C_computed = - Ri.transpose() * t;

    const Vec3 C_GT = d._C[i] - d._C[0];

    //-- Check that found camera position is equal to GT value
    if (i==0)  {
      EXPECT_MATRIX_NEAR(C_computed, C_GT, 1e-6);
    }
    else  {
     EXPECT_NEAR(0.0, DistanceLInfinity(C_computed.normalized(), C_GT.normalized()), 1e-6);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

