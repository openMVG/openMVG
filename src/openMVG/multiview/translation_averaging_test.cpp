// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"

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

TEST(translation_averaging, globalTi_from_tijs_Triplets_ECCV14) {

  int focal = 1000;
  int principal_Point = 500;
  //-- Setup a circular camera rig or "cardioid".
  const int iNviews = 12;
  const int iNbPoints = 6;
  NViewDataSet d =
    //NRealisticCamerasRing(
    NRealisticCamerasCardioid(
    iNviews, iNbPoints,
    nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0));

  d.ExportToPLY("global_translations_from_triplets_GT.ply");

  visibleCamPosToSVGSurface(d._C, "global_translations_from_triplets_GT.svg");

  // List successive triplets of the large loop of camera
  std::vector< graphUtils::Triplet > vec_triplets;
  for (size_t i = 0; i < iNviews; ++i)
  {
    const size_t iPlus1 = modifiedMod(i+1,iNviews);
    const size_t iPlus2 = modifiedMod(i+2,iNviews);
    //-- sort the triplet index to have a monotonic ascending series of value
    size_t triplet[3] = {i, iPlus1, iPlus2};
    std::sort(&triplet[0], &triplet[3]);
    vec_triplets.push_back(Triplet(triplet[0],triplet[1],triplet[2]));
  }

  //- For each triplet compute relative translations and rotations motions
  std::vector<openMVG::relativeInfo > vec_initialEstimates;

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
  //-   with the Kyle method
  std::vector<int> vec_edges;
  vec_edges.reserve(vec_initialEstimates.size() * 2);
  std::vector<double> vec_poses;
  vec_poses.reserve(vec_initialEstimates.size() * 3);
  std::vector<double> vec_weights;
  vec_weights.reserve(vec_initialEstimates.size());

  for(int i=0; i < vec_initialEstimates.size(); ++i)
  {
    const openMVG::relativeInfo & rel = vec_initialEstimates[i];
    vec_edges.push_back(rel.first.first);
    vec_edges.push_back(rel.first.second);

    const Vec3 EdgeDirection = -(d._R[rel.first.second].transpose() * rel.second.second.normalized());

    vec_poses.push_back(EdgeDirection(0));
    vec_poses.push_back(EdgeDirection(1));
    vec_poses.push_back(EdgeDirection(2));

    vec_weights.push_back(1.0);
  }

  const double function_tolerance=1e-7, parameter_tolerance=1e-8;
  const int max_iterations = 500;

  const double loss_width = 0.0;

  std::vector<double> X(iNviews*3);

  EXPECT_TRUE(
    solve_translations_problem(
      &vec_edges[0],
      &vec_poses[0],
      &vec_weights[0],
      vec_initialEstimates.size(),
      loss_width,
      &X[0],
      function_tolerance,
      parameter_tolerance,
      max_iterations));

  // Get back the camera translations in the global frame:
  for (size_t i = 0; i < iNviews; ++i)
  {
    if (i==0) {  //First camera supposed to be at Identity
      const Vec3 C0(X[0], X[1], X[2]);
      EXPECT_NEAR(0.0, DistanceLInfinity(C0, Vec3(0,0,0)), 1e-6);
    }
    else  {
      const Vec3 t_GT = (d._C[i] - d._C[0]);

      const Vec3 CI(X[i*3], X[i*3+1], X[i*3+2]);
      const Vec3 C0(X[0], X[1], X[2]);
      const Vec3 t_computed = CI - C0;

      //-- Check that vector are colinear
      EXPECT_NEAR(0.0, DistanceLInfinity(t_computed.normalized(), t_GT.normalized()), 1e-6);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

