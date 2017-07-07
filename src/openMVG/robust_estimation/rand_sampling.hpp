// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_HPP
#define OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_HPP

#include <algorithm>
#include <cstdlib>
#include <random>
#include <type_traits>
#include <vector>

namespace openMVG {
namespace robust{

/**
* Pick a random subset of the integers [0, total), in random order.
* Note that this can behave badly if num_samples is close to total; runtime
* could be unlimited!
*
* This uses a quadratic rejection strategy and should only be used for small
* num_samples.
*
* \param[in] num_samples      The number of samples to produce.
* \param[in] total_samples    The number of available samples.
* \param[in] random_generator The random number generator.
* \param[out] samples         num_samples of numbers in [0, total_samples) is placed
*                             here on return.
*/
template <class RandomGeneratorT, typename SamplingType>
inline void UniformSample
(
  const uint32_t num_samples,
  const uint32_t total_samples,
  RandomGeneratorT &random_generator,
  std::vector<SamplingType> *samples
)
{
  static_assert(std::is_integral<SamplingType>::value, "SamplingType must be an integral type");

  std::uniform_int_distribution<SamplingType> distribution(0, total_samples-1);
  samples->resize(0);
  while (samples->size() < num_samples) {
    const auto sample = distribution(random_generator);
    bool bFound = false;
    for (size_t j = 0; j < samples->size() && !bFound; ++j) {
      bFound = (*samples)[j] == sample;
    }
    if (!bFound) {
      samples->push_back(sample);
    }
  }
}

/**
* Pick n random sample from an array of indices.
* Use a Fisher Yates sampling (shuffling) to avoid picking the same index many time.
*
*
* \param[in] num_samples The number of randomly picked value in the vec_index array.
* \param[in] random_generator The random number generator.
* \param[out] vec_index An array of unique index value. The function shuffle this vector.
* \param[out] samples Output randomly picked value.
* \return true if the sampling can be performed
*/
template<typename T, class RandomGeneratorT, typename SamplingType = uint32_t>
bool UniformSample
(
  const size_t num_samples,
  RandomGeneratorT &random_generator,
  std::vector<T> * vec_index, // the array that provide the index (will be shuffled)
  std::vector<T> * samples // output found indices
)
{
  static_assert(std::is_integral<SamplingType>::value, "SamplingType must be an integral type");

  if (num_samples > vec_index->size() ||
      vec_index->size() > static_cast<size_t>(std::numeric_limits<SamplingType>::max()))
    return false;

  const SamplingType last_idx (vec_index->size() - 1);
  for (SamplingType i = 0; i < num_samples; ++i)
  {
    std::uniform_int_distribution<SamplingType> distribution(i, last_idx);
    const SamplingType sample = distribution(random_generator);

    std::swap((*vec_index)[i], (*vec_index)[sample]);
  }
  samples->resize(num_samples);
  for (size_t i=0; i < num_samples; ++i)
  {
    (*samples)[i] = (*vec_index)[i];
  }
  return true;
}


} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_HPP
