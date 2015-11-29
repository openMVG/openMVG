
// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_H_
#define OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_H_

#include <set>
#include <cstdlib>
#include <random>
#include <cassert>

namespace openMVG {
namespace robust{

/**
* Pick a random subset of the integers [0, total), in random order.
*
* This uses a quadratic rejection strategy and should only be used for small
* num_samples.
*
* \param num_samples   The number of samples to produce.
* \param total_samples The number of samples available.
* \param samples       num_samples of numbers in [0, total_samples) is placed
*                      here on return.
* \warning Argument values should respect: num_samples <= total_samples
*/
template<typename IntT>
static void UniformSample(
  size_t num_samples,
  size_t total_samples,
  std::set<IntT> *samples)
{
  assert(num_samples <= total_samples);
  std::random_device rd;
  std::default_random_engine e1(rd());
  std::uniform_int_distribution<IntT> uniform_dist(0, total_samples-1);
  while (samples->size() < num_samples)
  {
    IntT sample = uniform_dist(e1);
    samples->insert(sample);
  }
}

template<typename IntT>
static void UniformSample(
  size_t num_samples,
  size_t total_samples,
  std::vector<IntT> *samples)
{
  assert(num_samples <= total_samples);
  samples->resize(0);
  samples->reserve(num_samples);
  std::random_device rd;
  std::default_random_engine e1(rd());
  std::uniform_int_distribution<IntT> uniform_dist(0, total_samples-1);
  std::set<IntT> set_samples;
  while (set_samples.size() < num_samples)
  {
    IntT sample = uniform_dist(e1);
    if(set_samples.count(sample) == 0)
    {
      set_samples.insert(sample);
      samples->push_back(sample);
    }
  }
}

/// Get a (sorted) random sample of size X in [0:n-1]
static void random_sample(size_t X, size_t n, std::vector<size_t> *samples)
{
  samples->resize(X);
  for(size_t i=0; i < X; ++i) {
    size_t r = (std::rand()>>3)%(n-i), j;
    for(j=0; j<i && r>=(*samples)[j]; ++j)
      ++r;
    size_t j0 = j;
    for(j=i; j > j0; --j)
      (*samples)[j] = (*samples)[j-1];
    (*samples)[j0] = r;
  }
}

} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATION_RAND_SAMPLING_H_
