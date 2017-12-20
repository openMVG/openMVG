// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_TEST_HPP
#define OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_TEST_HPP

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "openMVG/graph/triplet_finder.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/numeric/numeric.h"

#include "testing/testing.h"
#include "third_party/vectorGraphics/svgDrawer.hpp"

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
  const std::vector<openMVG::Vec3> & vec_Ci,
  const std::string & fileName
)
{
  openMVG::Mat points(3, vec_Ci.size());
  for (size_t i = 0; i  < vec_Ci.size(); ++i)
  {
    points.col(i) = vec_Ci[i];
  }

  openMVG::Vec mean, variance;
  openMVG::MeanAndVarianceAlongRows(points, &mean, &variance);

  const double xfactor = sqrt(2.0 / variance(0));
  const double yfactor = sqrt(2.0 / variance(2));

  std::vector<openMVG::Vec3> out = vec_Ci;
  for (size_t i = 0; i  < vec_Ci.size(); ++i)
  {
    out[i](0) = ((out[i](0) * xfactor) + -xfactor * mean(0)) * 30 + 100;
    out[i](2) = ((out[i](2) * yfactor) + -yfactor * mean(2)) * 30 + 100;
  }

  if (!fileName.empty())
  {
    const double size = 200;
    svg::svgDrawer svgSurface_GT(size,size);
    for (size_t i = 0; i  < vec_Ci.size(); ++i)
    {
      svgSurface_GT.drawCircle(out[i](0), out[i](2),
                               3,svg::svgStyle().stroke("black",0.2).fill("red"));
    }
    std::ostringstream osSvgGT;
    osSvgGT << fileName;
    std::ofstream svgFileGT( osSvgGT.str().c_str());
    svgFileGT << svgSurface_GT.closeSvgFile().str();
  }
}

openMVG::NViewDataSet Setup_RelativeTranslations_AndNviewDataset
(
  std::vector<openMVG::RelativeInfo_Vec > & vec_relative_estimates,
  const int focal = 1000,
  const int principal_Point = 500,
  const int iNviews = 12,
  const int iNbPoints = 6,
  const bool bCardiod = true,
  const bool bRelative_Translation_PerTriplet = false
)
{
  const openMVG::NViewDataSet d =
    bCardiod ?
      NRealisticCamerasCardioid(
        iNviews, iNbPoints,
        openMVG::nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0))
      :NRealisticCamerasRing(
        iNviews, iNbPoints,
        openMVG::nViewDatasetConfigurator(focal,focal,principal_Point,principal_Point,5,0));

  // Init the relative pair of motion depending of the asked configuration:
  // -2-view: bearing direction,
  // -3-view: triplet of bearing direction.
  std::vector<openMVG::Pair_Vec > motion_group;
  if (bRelative_Translation_PerTriplet)
  {
    std::vector<openMVG::graph::Triplet> vec_triplets;
    for (int i = 0; i < iNviews; ++i)
    {
      const openMVG::IndexT iPlus1 = modifiedMod(i+1,iNviews);
      const openMVG::IndexT iPlus2 = modifiedMod(i+2,iNviews);
      //-- sort the triplet index to have a monotonic ascending series of value
      openMVG::IndexT triplet[3] = {openMVG::IndexT(i), iPlus1, iPlus2};
      std::sort(&triplet[0], &triplet[3]);
      vec_triplets.emplace_back(triplet[0],triplet[1],triplet[2]);
    }
    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const openMVG::graph::Triplet & triplet = vec_triplets[i];
      const size_t I = triplet.i, J = triplet.j , K = triplet.k;

      openMVG::Pair_Vec pairs;
      pairs.emplace_back(I,J);
      pairs.emplace_back(J,K);
      pairs.emplace_back(I,K);
      motion_group.push_back(pairs);
    }
  }
  else
  {
    //-- Setup initial pairs that will be considered (link each camera to the two next)
    for (int i = 0; i < iNviews; ++i)
    {
      for (int j=i; j<=i+2; ++j)
      {
        const int jj = modifiedMod(j,iNviews);
        if (i != jj)
        {
          const openMVG::Pair_Vec pairs = {openMVG::Pair(i,jj)};
          motion_group.push_back(pairs);
        }
      }
    }
  }
  // Compute all the required relative motions
  for (const auto & iterGroup : motion_group)
  {
    openMVG::RelativeInfo_Vec relative_motion;
    for (const openMVG::Pair & iter : iterGroup)
    {
      const size_t I = iter.first;
      const size_t J = iter.second;

      //-- Build camera alias over GT translations and rotations:
      const openMVG::Mat3 & RI = d._R[I];
      const openMVG::Mat3 & RJ = d._R[J];
      const openMVG::Vec3 & ti = d._t[I];
      const openMVG::Vec3 & tj = d._t[J];

      //-- Build relative motions (that feeds the Linear program formulation)
      {
        openMVG::Mat3 RijGt;
        openMVG::Vec3 tij;
        openMVG::RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
        relative_motion.emplace_back(std::make_pair(I, J), std::make_pair(RijGt, tij));
      }
    }
    vec_relative_estimates.push_back(relative_motion);
  }
  return d;
}

#endif // OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_TEST_HPP
