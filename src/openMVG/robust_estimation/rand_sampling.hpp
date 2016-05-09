
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_H_
#define OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_H_

#include <random>
#include <stdlib.h>
#include <vector>

namespace openMVG {
namespace robust{

namespace
{
  std::default_random_engine random_generator;
}

/**
* Pick a random subset of the integers [0, total), in random order.
* Note that this can behave badly if num_samples is close to total; runtime
* could be unlimited!
*
* This uses a quadratic rejection strategy and should only be used for small
* num_samples.
*
* \param num_samples   The number of samples to produce.
* \param total_samples The number of available samples.
* \param samples       num_samples of numbers in [0, total_samples) is placed
*                      here on return.
*/
inline void UniformSample
(
  size_t num_samples,
  size_t total_samples,
  std::vector<size_t> *samples
)
{
  std::uniform_int_distribution<unsigned int> distribution(0, total_samples-1);
  samples->resize(0);
  while (samples->size() < num_samples) {
    const unsigned int sample = distribution(random_generator);
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
* \param num_samples The number of randomly picked value in the vec_index array.
* \param vec_index An array of unique index value. The function shuffle this vector.
* \param samples Output randomly picked value.
* \return true if the sampling can be performed
*/
//
template<typename T>
inline bool UniformSample
(
  const uint32_t num_samples,
  std::vector<T> * vec_index, // the array that provide the index (will be shuffled)
  std::vector<T> * samples // output found indices
)
{
  if (vec_index->size() < num_samples)
    return false;

  const uint32_t last_idx (vec_index->size() - 1);
  for (uint32_t i = 0; i < num_samples; ++i)
  {
    std::uniform_int_distribution<uint32_t> distribution(i, last_idx);
    const uint32_t sample = distribution(random_generator);

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
#endif // OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_H_
