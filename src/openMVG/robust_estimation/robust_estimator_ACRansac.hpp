
// Copyright (c) 2012, 2013 Lionel MOISAN.
// Copyright (c) 2012, 2013 Pascal MONASSE.
// Copyright (c) 2012, 2016 Pierre MOULON.

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
//  the homography) and extended to be used in large scale Structure from
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
#include "third_party/histogram/histogram.hpp"

namespace openMVG {
namespace robust{

namespace acransac_nfa_internal {

/// logarithm (base 10) of binomial coefficient
template <typename T>
static T logcombi
(
  size_t k,
  size_t n,
  const std::vector<T> & vec_log10 // lookuptable in [0,n+1]
)
{
  if (k>=n || k<=0) return T(0);
  if (n-k<k) k=n-k;
  T r(0);
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
  std::vector<Type> & vec_logc_n
)
{
  // compute a lookuptable of log10 value for the range [0,n+1]
  std::vector<Type> vec_log10(n + 1);
  for (size_t k = 0; k <= n; ++k)
    vec_log10[k] = log10((Type)k);

  makelogcombi_n(n, vec_logc_n, vec_log10);
  makelogcombi_k(k, n, vec_logc_k, vec_log10);
}

template <typename Kernel>
class NFA_Interface
{
public:
  /**
   * @brief NFA_Interface constructor
   * @param[in] kernel Template kernel model estimator & residual error evaluation interface
   * @param[in] dmaxThreshold Upper bound of the residual error (default infinity)
   * @param[in] bquantified_nfa_evaluation Tell if NFA evaluation is using the quantified or exhaustive evaluation method.
   *  An upper bound different from infinity must be provided to be set to true.
   */
  NFA_Interface
  (
    const Kernel & kernel,
    const double dmaxThreshold = std::numeric_limits<double>::infinity(),
    const bool bquantified_nfa_evaluation = false
  ):
    m_residuals(kernel.NumSamples()),
    m_kernel(kernel),
    m_bquantified_nfa_evaluation(bquantified_nfa_evaluation),
    m_max_threshold(dmaxThreshold)
  {
    // Precompute log combi
    m_loge0 = log10((double)Kernel::MAX_MODELS * (kernel.NumSamples()-Kernel::MINIMUM_SAMPLES));
    makelogcombi(Kernel::MINIMUM_SAMPLES, kernel.NumSamples(), m_logc_k, m_logc_n);
  };

  std::vector<double> & residuals()
  { return m_residuals;}

  /**
   * @brief Evaluation of the NFA (Number of False Alarm)
   *  for the given residual distribution.
   * The NFA can be evaluated in two way:
   * 1. If an upper bound of the threshold is provided:
   *  - The NFA is estimated by using quantified residual values.
   * 2. No upper bound => m_max_threshold == infinity:
   *  - The NFA is estimated by using all the residual errors.
   *
   * @param[out] inliers inlier indices list (updated if a better NFA is found)
   * @param[in, out] nfa_threshold Found NFA and corresponding error Threshold
   *  (updated if the current estimated NFA is lower than the existing one)
   *  For the first run it can be set to {std::numeric_limits<double>::infinity(), 0.0}
   *
   * @return true if a better NFA is found.
   */
  bool ComputeNFA_and_inliers
  (
    std::vector<size_t> & inliers,
    std::pair<double,double> & nfa_threshold
  );

private:

  /// residual array
  std::vector<double> m_residuals;
  /// [residual,index] array -> used in the exhaustive nfa computation mode
  std::vector<std::pair<double,int> > m_sorted_residuals;

  /// Combinatorial log
  std::vector<float> m_logc_n, m_logc_k;
  /// A-Contrario Epsilon 0 value
  double m_loge0;

  /// Kernel (model estimation interface)
  const Kernel & m_kernel;
  /// Tell if the NFA is computed in the quantified or "exhaustive" mode
  const bool m_bquantified_nfa_evaluation;
  /// upper bound of the maximum authorized residual value
  const double m_max_threshold;
};

template <typename Kernel>
bool
NFA_Interface<Kernel>::ComputeNFA_and_inliers
(
    std::vector<size_t> & inliers,
    /// NFA and residual threshold
    std::pair<double,double> & nfa_threshold
)
{
  // A-Contrario computation of the most meaningful discrimination inliers/outliers.
  // Two computation mode are implemented:
  // - A quantified computation
  //    (valuable is an upper bound of the maximal tolerated residual is provided)
  // - An exhaustive computation that evaluate all the possible NFA values
  //    i.e. (for every k from the n value of the datum).
  if (m_bquantified_nfa_evaluation)
  {
    // Find the best NFA (Number of False Alarm) score
    //  by using residual errors sorted in a histogram.
    // This version avoid:
    //   - to sort explicitly the residual error array,
    //   - to compute the NFA for every sample of the datum.
    const int nBins = 20;
    Histogram<double> histo(0.0f, m_max_threshold, nBins);
    histo.Add(m_residuals.begin(), m_residuals.end());

    // Compute NFA scoring from the cumulative histogram

    typedef std::pair<double,double> nfa_thresholdT; // NFA and residual threshold
    nfa_thresholdT current_best_nfa(std::numeric_limits<double>::infinity(), 0.0);
    unsigned int cumulative_count = 0;
    const std::vector<size_t> & frequencies = histo.GetHist();
    const std::vector<double> residual_val = histo.GetXbinsValue();
    for (int bin = 0; bin < nBins; ++bin)
    {
      cumulative_count += frequencies[bin];
      if (cumulative_count > Kernel::MINIMUM_SAMPLES
          && residual_val[bin] > std::numeric_limits<float>::epsilon())
      {
        const double logalpha = m_kernel.logalpha0()
          + m_kernel.multError() * log10(residual_val[bin]
          + std::numeric_limits<float>::epsilon());
        const nfa_thresholdT current_nfa( m_loge0
          + logalpha * (double)(cumulative_count - Kernel::MINIMUM_SAMPLES)
          + m_logc_n[cumulative_count]
          + m_logc_k[cumulative_count], residual_val[bin]);
        // Keep the best NFA iff it is meaningful ( NFA < 0 ) and better than the existing one
        if(current_nfa.first < current_best_nfa.first && current_nfa.first < 0)
          current_best_nfa = current_nfa;
      }
    }
    // If the current NFA is better than the previous
    // - update the sample inlier index list.
    if (current_best_nfa.first < nfa_threshold.first)
    {
      nfa_threshold.first = current_best_nfa.first; // NFA score
      nfa_threshold.second = current_best_nfa.second; // Corresponding threshold

      inliers.clear();
      for (size_t index = 0; index < m_kernel.NumSamples(); ++index)
      {
        if (m_residuals[index] <= nfa_threshold.second)
          inliers.push_back(index);
      }
      return inliers.size() > Kernel::MINIMUM_SAMPLES;
    }
  }
  else // exhaustive computation
  {
    // Residuals sorting (ascending order while keeping original point indexes)
    {
      m_sorted_residuals.clear();
      m_sorted_residuals.reserve(m_kernel.NumSamples());
      for (size_t i = 0; i < m_kernel.NumSamples(); ++i)
      {
        m_sorted_residuals.emplace_back(m_residuals[i], i);
      }
      std::sort(m_sorted_residuals.begin(), m_sorted_residuals.end());
    }

    // Find best NFA and its index wrt square error threshold in m_sorted_residuals.
    typedef std::pair<double,int> nfa_indexT;
    nfa_indexT current_best_nfa(std::numeric_limits<double>::infinity(), Kernel::MINIMUM_SAMPLES);
    const size_t n = m_kernel.NumSamples();
    for(size_t k=Kernel::MINIMUM_SAMPLES+1;
        k<=n && m_sorted_residuals[k-1].first<=m_max_threshold;
        ++k) // Compute the NFA for all k in [minimal_sample+1,n]
    {
      const double logalpha = m_kernel.logalpha0()
        + m_kernel.multError() * log10(m_sorted_residuals[k-1].first
        + std::numeric_limits<float>::epsilon());
      const nfa_indexT current_nfa( m_loge0
        + logalpha * (double)(k - Kernel::MINIMUM_SAMPLES)
        + m_logc_n[k]
        + m_logc_k[k], k);

      if(current_nfa.first < current_best_nfa.first)
        current_best_nfa = current_nfa;
    }

    // If the current NFA is better than the previous
    // - update the sample inlier index list.
    if (current_best_nfa.first < nfa_threshold.first)
    {
      nfa_threshold.first = current_best_nfa.first;
      nfa_threshold.second = m_sorted_residuals[current_best_nfa.second-1].first;

      inliers.resize(current_best_nfa.second);
      for (size_t i=0; i<current_best_nfa.second; ++i)
      {
        inliers[i] = m_sorted_residuals[i].second;
      }
      return true;
    }
  }
  return false;
}
}  // namespace acransac_nfa_internal

/**
 * @brief ACRANSAC routine (ErrorThreshold, NFA)
 * If an upper bound of the threshold is provided:
 *  - The NFA is estimated by using quantified residual values.
 * Else (no upper bound => precision == infinity):
 *  - The NFA is estimated by using all the residual errors.
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
  const unsigned int num_max_iteration = 1024,
  typename Kernel::Model * model = NULL,
  double precision = std::numeric_limits<double>::infinity(),
  bool bVerbose = false)
{
  vec_inliers.clear();

  const unsigned int sizeSample = Kernel::MINIMUM_SAMPLES;
  const unsigned int nData = kernel.NumSamples();
  if (nData <= sizeSample)
    return std::make_pair(0.0,0.0);

  //--
  // Sampling:
  // Possible sampling indices [0,..,nData] (will change in the optimization phase)
  std::vector<size_t> vec_index(nData);
  std::iota(vec_index.begin(), vec_index.end(), 0);
  // Sample indices (used for model evaluation)
  std::vector<size_t> vec_sample(sizeSample);

  const double maxThreshold = (precision==std::numeric_limits<double>::infinity()) ?
    std::numeric_limits<double>::infinity() :
    precision * kernel.normalizer2()(0,0) * kernel.normalizer2()(0,0);

  // Initialize the NFA computation interface
  // (quantified NFA computation is used if a valid upper bound is provided)
  acransac_nfa_internal::NFA_Interface<Kernel> nfa_interface
    (kernel, maxThreshold, (precision!=std::numeric_limits<double>::infinity()));

  // Output parameters
  double minNFA = std::numeric_limits<double>::infinity();
  double errorMax = std::numeric_limits<double>::infinity();

  //--
  // Local optimization:
  // Reserve 10% of iterations for focused sampling
  int nIterReserve = num_max_iteration/10;
  unsigned int nIter = num_max_iteration - nIterReserve;

  //--
  // An early exit is used when an upper bound threshold is provided:
  // - If a model can be found with a valid point support by using MAX-CONSENSUS
  //    on a short number of iteration, then we enable AC-RANSAC,
  //    else we do an early exit since there is a very few chance to find a valid model.
  bool bACRansacMode = (precision == std::numeric_limits<double>::infinity());

  //--
  // Main estimation loop.
  for (unsigned int iter=0; iter < nIter && iter < num_max_iteration; ++iter)
  {
    // Get random samples
    if (bACRansacMode)
      UniformSample(sizeSample, &vec_index, &vec_sample);
    else
      UniformSample(sizeSample, nData, &vec_sample);

    // Fit model(s). Can find up to Kernel::MAX_MODELS solution(s)
    std::vector<typename Kernel::Model> vec_models;
    kernel.Fit(vec_sample, &vec_models);

    // Evaluate model(s)
    bool better = false;
    for (unsigned int k = 0; k < static_cast<unsigned int>(vec_models.size()); ++k)
    {
      // Compute residual values
      kernel.Errors(vec_models[k], nfa_interface.residuals());

      if (!bACRansacMode)
      {
        // MAX-CONSENSUS checking (does a model with some support is existing)
        unsigned int nInlier = 0;
        for (size_t i = 0; i < nData; ++i)
        {
         if (nfa_interface.residuals()[i] <= maxThreshold)
          ++nInlier;
        }
        if (nInlier > 2.5 * sizeSample) // does the model is meaningful
          bACRansacMode = true;
      }

      if (bACRansacMode)
      {
        // NFA evaluation; If better than the previous: update scoring & inliers indices
        std::pair<double,double> nfa_threshold(minNFA, 0.0);
        const bool b_better_model_found =
          nfa_interface.ComputeNFA_and_inliers(vec_inliers, nfa_threshold);

        if (b_better_model_found)
        {
          better = true;
          minNFA = nfa_threshold.first;
          errorMax = nfa_threshold.second;
          if(model) *model = vec_models[k];

          if(bVerbose)
          {
            std::cout << "  nfa=" << minNFA
              << " inliers=" << vec_inliers.size() << "/" << nData
              << " precisionNormalized=" << errorMax
              << " precision=" << kernel.unormalizeError(errorMax)
              << " (iter=" << iter
              << " ,sample=";
            std::copy(vec_sample.begin(), vec_sample.end(),
              std::ostream_iterator<size_t>(std::cout, ","));
            std::cout << ")" << std::endl;
          }
        }
      }
    }

    // Early exit test -> no meaningful model found so far
    //  see explanation above
    if (!bACRansacMode && iter > nIterReserve*2)
    {
      nIter = 0; // No more round will be performed
      continue;
    }

    // ACRANSAC optimization: draw samples among best set of inliers so far
    if (bACRansacMode && ((better && minNFA<0) || (iter+1==nIter && nIterReserve > 0)))
    {
      if (vec_inliers.empty())
      {
        // No model found at all so far
        ++nIter; // Continue to look for any model, even not meaningful
        --nIterReserve;
      }
      else
      {
        // ACRANSAC optimization: draw samples among best set of inliers so far
        vec_index = vec_inliers;
        if(nIterReserve) {
            // reduce the number of iteration
            // next iterations will be dedicated to local optimization
            nIter = iter + 1 + nIterReserve;
            nIterReserve = 0;
        }
      }
    }
  }

  if(minNFA >= 0) // no meaningful model found so far
    vec_inliers.clear();

  if (!vec_inliers.empty())
  {
    // Un-normalize the model and the associated NFA threshold
    if (model)
      kernel.Unnormalize(model);
    errorMax = kernel.unormalizeError(errorMax);
  }

  return std::make_pair(errorMax, minNFA);
}

} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_H_
