
// Copyright (c) 2012, 2013 Lionel MOISAN.
// Copyright (c) 2012, 2013 Pascal MONASSE.
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_H_
#define OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_H_

//-------------------
// Generic implementation of ACRANSAC
//-------------------
// The A contrario parametrization have been first explained in [1] and
//  later extended to generic model estimation in [2] (with a demonstration for
//  the homography) and extended and use at large scale for Structure from
//  Motion in [3].
//
//--
//  [1] Lionel Moisan, Berenger Stival,
//  A probalistic criterion to detect rigid point matches between
//  two images and estimate the fundamental matrix.
//  IJCV 04.
//--
//  [2] Lionel Moisan, Pierre Moulon, Pascal Monasse.
//  Automatic Homographic Registration of a Pair of Images,
//    with A Contrario Elimination of Outliers
//  Image Processing On Line (IPOL), 2012.
//  http://dx.doi.org/10.5201/ipol.2012.mmm-oh
//--
//  [3] Pierre Moulon, Pascal Monasse and Renaud Marlet.
//  Adaptive Structure from Motion with a contrario mode estimation.
//  In 11th Asian Conference on Computer Vision (ACCV 2012)
//--


#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <vector>

#include "openMVG/robust_estimation/rand_sampling.hpp"

namespace openMVG {
namespace robust{

/// logarithm (base 10) of binomial coefficient
template <typename T>
static T logcombi
(
  size_t k,
  size_t n,
  const std::vector<T> & vec_log10 // lookuptable in [0,n+1]
)
{
  if (k>=n || k<=0) return(0.0f);
  if (n-k<k) k=n-k;
  T r = 0.0f;
  for (size_t i = 1; i <= k; ++i)
    r += vec_log10[n-i+1] - vec_log10[i];
  return r;
}

/// tabulate logcombi(.,n)
template<typename Type>
static void makelogcombi_n
(
  size_t n,
  std::vector<Type> & l,
  std::vector<Type> & vec_log10 // lookuptable [0,n+1]
)
{
  l.resize(n+1);
  for (size_t k = 0; k <= n; ++k)
    l[k] = logcombi<Type>(k, n, vec_log10);
}

/// tabulate logcombi(k,.)
template<typename Type>
static void makelogcombi_k
(
  size_t k,
  size_t nmax,
  std::vector<Type> & l,
  std::vector<Type> & vec_log10 // lookuptable [0,n+1]
)
{
  l.resize(nmax+1);
  for (size_t n = 0; n <= nmax; ++n)
    l[n] = logcombi<Type>(k, n, vec_log10);
}

template <typename Type>
static void makelogcombi
(
  size_t k,
  size_t n,
  std::vector<Type> & vec_logc_k,
  std::vector<Type> & vec_logc_n)
{
  // compute a lookuptable of log10 value for the range [0,n+1]
  std::vector<Type> vec_log10(n + 1);
  for (size_t k = 0; k <= n; ++k)
    vec_log10[k] = log10((Type)k);

  makelogcombi_n(n, vec_logc_n, vec_log10);
  makelogcombi_k(k, n, vec_logc_k, vec_log10);
}


/// Distance and associated index
typedef std::pair<double,size_t> ErrorIndex;

/// Find best NFA and its index wrt square error threshold in e.
static ErrorIndex bestNFA(
  int startIndex, //number of point required for estimation
  double logalpha0,
  const std::vector<ErrorIndex>& e,
  double loge0,
  double maxThreshold,
  const std::vector<float> &logc_n,
  const std::vector<float> &logc_k,
  double multError = 1.0)
{
  ErrorIndex bestIndex(std::numeric_limits<double>::infinity(), startIndex);
  const size_t n = e.size();
  for(size_t k=startIndex+1; k<=n && e[k-1].first<=maxThreshold; ++k) {
    const double logalpha = logalpha0 + multError * log10(e[k-1].first + std::numeric_limits<float>::min());
    ErrorIndex index(loge0 +
      logalpha * (double)(k-startIndex) +
      logc_n[k] +
      logc_k[k], k);

    if(index.first < bestIndex.first)
      bestIndex = index;
  }
  return bestIndex;
}

/// Pick a random sample
/// \param sizeSample The size of the sample.
/// \param vec_index  The possible data indices.
/// \param sample The random sample of sizeSample indices (output).
static void UniformSample(int sizeSample,
  const std::vector<size_t> &vec_index,
  std::vector<size_t> *sample)
{
  sample->resize(sizeSample);
  random_sample(sizeSample, vec_index.size(), sample);
  for(int i = 0; i < sizeSample; ++i)
    (*sample)[i] = vec_index[ (*sample)[i] ];
}

/**
 * @brief ACRANSAC routine (ErrorThreshold, NFA)
 *
 * @param[in] kernel model and metric object
 * @param[out] vec_inliers points that fit the estimated model
 * @param[in] nIter maximum number of consecutive iterations
 * @param[out] model returned model if found
 * @param[in] precision upper bound of the precision (squared error)
 * @param[in] bVerbose display console log
 *
 * @return (errorMax, minNFA)
 */
template<typename Kernel>
std::pair<double, double> ACRANSAC(const Kernel &kernel,
  std::vector<size_t> & vec_inliers,
  size_t nIter = 1024,
  typename Kernel::Model * model = NULL,
  double precision = std::numeric_limits<double>::infinity(),
  bool bVerbose = false)
{
  vec_inliers.clear();

  const size_t sizeSample = Kernel::MINIMUM_SAMPLES;
  const size_t nData = kernel.NumSamples();
  if(nData <= (size_t)sizeSample)
    return std::make_pair(0.0,0.0);

  const double maxThreshold = (precision==std::numeric_limits<double>::infinity()) ?
    std::numeric_limits<double>::infinity() :
    precision * kernel.normalizer2()(0,0) * kernel.normalizer2()(0,0);

  std::vector<ErrorIndex> vec_residuals(nData); // [residual,index]
  std::vector<double> vec_residuals_(nData);
  std::vector<size_t> vec_sample(sizeSample); // Sample indices

  // Possible sampling indices [0,..,nData] (will change in the optimization phase)
  std::vector<size_t> vec_index(nData);
  std::iota(vec_index.begin(), vec_index.end(), 0);

  // Precompute log combi
  const double loge0 = log10((double)Kernel::MAX_MODELS * (nData-sizeSample));
  std::vector<float> vec_logc_n, vec_logc_k;
  makelogcombi(sizeSample, nData, vec_logc_k, vec_logc_n);

  // Output parameters
  double minNFA = std::numeric_limits<double>::infinity();
  double errorMax = std::numeric_limits<double>::infinity();

  // Reserve 10% of iterations for focused sampling
  size_t nIterReserve = nIter/10;
  nIter -= nIterReserve;

  bool bACRansacMode = (precision == std::numeric_limits<double>::infinity());

  // Main estimation loop.
  for (size_t iter=0; iter < nIter; ++iter) {
    if (bACRansacMode)
      UniformSample(sizeSample, vec_index, &vec_sample); // Get random sample
    else
      UniformSample(sizeSample, nData, &vec_sample); // Get random sample

    std::vector<typename Kernel::Model> vec_models; // Up to max_models solutions
    kernel.Fit(vec_sample, &vec_models);

    // Evaluate models
    bool better = false;
    for (size_t k = 0; k < vec_models.size(); ++k)  {
      // Residuals computation and ordering
      kernel.Errors(vec_models[k], vec_residuals_);

      if (!bACRansacMode)
      {
        unsigned int nInlier = 0;
        for (size_t i = 0; i < nData; ++i)
        {
         if (vec_residuals_[i] <= maxThreshold)
          ++nInlier;
        }
        if (nInlier > 2.5 * sizeSample) // does the model is meaningful
          bACRansacMode = true;
        if (!bACRansacMode && nIter > nIterReserve)
          nIter = 0;
      }
      if (bACRansacMode)
      {
        for (size_t i = 0; i < nData; ++i)  {
          const double error = vec_residuals_[i];
          vec_residuals[i] = ErrorIndex(error, i);
        }
        std::sort(vec_residuals.begin(), vec_residuals.end());

        // Most meaningful discrimination inliers/outliers
        const ErrorIndex best = bestNFA(
          sizeSample,
          kernel.logalpha0(),
          vec_residuals,
          loge0,
          maxThreshold,
          vec_logc_n,
          vec_logc_k,
          kernel.multError());

        if (best.first < minNFA /*&& vec_residuals[best.second-1].first < errorMax*/)  {
          // A better model was found
          better = true;
          minNFA = best.first;
          vec_inliers.resize(best.second);
          for (size_t i=0; i<best.second; ++i)
            vec_inliers[i] = vec_residuals[i].second;
          errorMax = vec_residuals[best.second-1].first; // Error threshold
          if(model) *model = vec_models[k];

          if(bVerbose)  {
            std::cout << "  nfa=" << minNFA
              << " inliers=" << best.second << "/" << nData
              << " precisionNormalized=" << errorMax
              << " precision=" << kernel.unormalizeError(errorMax)
              << " (iter=" << iter;
            std::cout << ",sample=";
            std::copy(vec_sample.begin(), vec_sample.end(),
              std::ostream_iterator<size_t>(std::cout, ","));
            std::cout << ")" <<std::endl;
          }
        }
      }
    }

    // ACRANSAC optimization: draw samples among best set of inliers so far
    if (bACRansacMode && ((better && minNFA<0) || (iter+1==nIter && nIterReserve))) {
      if (vec_inliers.empty()) { // No model found at all so far
        ++nIter; // Continue to look for any model, even not meaningful
        --nIterReserve;
      } else {
        // ACRANSAC optimization: draw samples among best set of inliers so far
        vec_index = vec_inliers;
        if(nIterReserve) {
            nIter = iter + 1 + nIterReserve;
            nIterReserve = 0;
        }
      }
    }
  }

  if(minNFA >= 0)
    vec_inliers.clear();

  if (!vec_inliers.empty())
  {
    if (model)
      kernel.Unnormalize(model);
    errorMax = kernel.unormalizeError(errorMax);
  }

  return std::make_pair(errorMax, minNFA);
}

} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_H_
