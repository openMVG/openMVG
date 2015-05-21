
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/features/feature.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/matching/indMatch.hpp"

using namespace openMVG;

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <vector>
#include <map>

using namespace openMVG::matching;

class ImageCollectionGeometricFilter
{
  public:
    ImageCollectionGeometricFilter() : _feat_provider(NULL)
  {
  }

  ImageCollectionGeometricFilter(Features_Provider * feat_provider)
    : _feat_provider(feat_provider)
  { }

  /// Filter all putative correspondences according the templated geometric filter
  template <typename GeometricFilterT>
  void Filter(
    const GeometricFilterT & geometricFilter,  // geometric filter functor
    PairWiseMatches & map_PutativesMatchesPair, // putative correspondences to filter
    PairWiseMatches & map_GeometricMatches,
    const std::vector<std::pair<size_t, size_t> > & vec_imagesSize) const
  {
    C_Progress_display my_progress_bar( map_PutativesMatchesPair.size() );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < (int)map_PutativesMatchesPair.size(); ++i)
    {
      PairWiseMatches::const_iterator iter = map_PutativesMatchesPair.begin();
      advance(iter,i);

      const size_t iIndex = iter->first.first;
      const size_t jIndex = iter->first.second;
      const std::vector<IndMatch> & vec_PutativeMatches = iter->second;

      // Load features of Inth and Jnth images
      const features::PointFeatures & kpSetI = _feat_provider->getFeatures(iIndex);
      const features::PointFeatures & kpSetJ = _feat_provider->getFeatures(jIndex);

      //-- Copy point to array in order to estimate fundamental matrix :
      const size_t n = vec_PutativeMatches.size();
      Mat xI(2,n), xJ(2,n);

      for (size_t i=0; i < vec_PutativeMatches.size(); ++i)  {
        const features::PointFeature & imaA = kpSetI[vec_PutativeMatches[i]._i];
        const features::PointFeature & imaB = kpSetJ[vec_PutativeMatches[i]._j];
        xI.col(i) = Vec2f(imaA.coords()).cast<double>();
        xJ.col(i) = Vec2f(imaB.coords()).cast<double>();
      }

      //-- Apply the geometric filter
      {
        std::vector<size_t> vec_inliers;
        geometricFilter.Fit(
          iter->first,
          xI, vec_imagesSize[iIndex],
          xJ, vec_imagesSize[jIndex],
          vec_inliers);

        if(!vec_inliers.empty())
        {
          std::vector<IndMatch> vec_filteredMatches;
          vec_filteredMatches.reserve(vec_inliers.size());
          for (size_t i=0; i < vec_inliers.size(); ++i)  {
            vec_filteredMatches.push_back( vec_PutativeMatches[vec_inliers[i]] );
          }
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
          {
            map_GeometricMatches[std::make_pair(iIndex,jIndex)] = vec_filteredMatches;
          }
        }
      }
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
      {
        ++my_progress_bar;
      }
    }
  }

  private:
  // Features per image
  Features_Provider * _feat_provider;
};

