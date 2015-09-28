// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"

#include "openMVG/graph/triplet_finder.hpp"
using namespace openMVG::graph;

#include "third_party/vectorGraphics/svgDrawer.hpp"
using namespace svg;

#include "openMVG/multiview/test_data_sets.hpp"
#include "testing/testing.h"

#include <fstream>
#include <map>
#include <utility>
#include <vector>

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include <numeric>

using namespace openMVG;
using namespace std;

int modifiedMod
(
  int number,
  int modulus
)
{
  int result = number % modulus;
  if (result < 0) result += modulus;
  return result;
}

//-- Export a series of camera positions to a SVG surface of specified squared size
void visibleCamPosToSVGSurface
(
  const std::vector<Vec3> & vec_Ci,
  const std::string & fileName
)
{
  Mat points(3, vec_Ci.size());
  for(size_t i = 0; i  < vec_Ci.size(); ++i)
  {
    points.col(i) = vec_Ci[i];
  }

  Vec mean, variance;
  MeanAndVarianceAlongRows(points, &mean, &variance);

  const double xfactor = sqrt(2.0 / variance(0));
  const double yfactor = sqrt(2.0 / variance(2));

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

NViewDataSet Setup_RelativeTranslations_AndNviewDataset
(
  std::vector<openMVG::relativeInfo > & vec_relative_estimates,
  const int focal = 1000,
  const int principal_Point = 500,
  const int iNviews = 12,
  const int iNbPoints = 6,
  const bool bCardiod = true,
  const bool bRelative_Translation_PerTriplet = false
)
{
  const NViewDataSet d =
    bCardiod ?
      NRealisticCamerasCardioid(
        iNviews, iNbPoints,
        nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0))
      :NRealisticCamerasRing(
        iNviews, iNbPoints,
        nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0));

  // Init the relative pair of motion depending of the asked configuration:
  // -2-view: bearing direction,
  // -3-view: triplet of bearing direction.
  std::vector< std::pair<size_t,size_t> > map_pairs;
  if (bRelative_Translation_PerTriplet)
  {
    std::vector< graph::Triplet > vec_triplets;
    for (IndexT i = 0; i < iNviews; ++i)
    {
      const IndexT iPlus1 = modifiedMod(i+1,iNviews);
      const IndexT iPlus2 = modifiedMod(i+2,iNviews);
      //-- sort the triplet index to have a monotonic ascending series of value
      IndexT triplet[3] = {i, iPlus1, iPlus2};
      std::sort(&triplet[0], &triplet[3]);
      vec_triplets.push_back(Triplet(triplet[0],triplet[1],triplet[2]));
    }
    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graph::Triplet & triplet = vec_triplets[i];
      const size_t I = triplet.i, J = triplet.j , K = triplet.k;

      map_pairs.emplace_back(I,J);
      map_pairs.emplace_back(J,K);
      map_pairs.emplace_back(I,K);
    }
  }
  else
  {
    //-- Setup initial pairs that will be considered (link each camera to the two next)
    for (size_t i = 0; i < iNviews; ++i)
    {
      for (size_t j=i; j<=i+2; ++j)
      {
        const size_t jj = modifiedMod(j,iNviews);
        if (i != jj)
          map_pairs.emplace_back(i,jj);
      }
    }
  }
  // Compute all the required relative motions
  for (const std::pair<size_t,size_t> & iter : map_pairs)
  {
    const size_t I = iter.first;
    const size_t J = iter.second;

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
      vec_relative_estimates.push_back(
        std::make_pair(std::make_pair(I, J), std::make_pair(RijGt, tij)));
    }
  }
  return d;
}

