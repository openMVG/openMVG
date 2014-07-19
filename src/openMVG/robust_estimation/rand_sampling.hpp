
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

#include <vector>
#include <stdlib.h>

namespace openMVG {
namespace robust{

using namespace std;

/**
* Pick a random subset of the integers [0, total), in random order.
* Note that this can behave badly if num_samples is close to total; runtime
* could be unlimited!
*
* This uses a quadratic rejection strategy and should only be used for small
* num_samples.
*
* \param num_samples   The number of samples to produce.
* \param total_samples The number of samples available.
* \param samples       num_samples of numbers in [0, total_samples) is placed
*                      here on return.
*/
static void UniformSample(
  size_t num_samples,
  size_t total_samples,
  std::vector<size_t> *samples)
{
  samples->resize(0);
  while (samples->size() < num_samples) {
    size_t sample = size_t(rand() % total_samples);
    bool bFound = false;
    for (size_t j = 0; j < samples->size(); ++j) {
      bFound = (*samples)[j] == sample;
      if (bFound) { //the picked index already exist
        break;
      }
    }
    if (!bFound) {
      samples->push_back(sample);
    }
  }
}

/// Get a (sorted) random sample of size X in [0:n-1]
/// samples array must be pre-allocated
static void random_sample(size_t X, size_t n, std::vector<size_t> *samples)
{
  samples->resize(X);
  for(size_t i=0; i < X; ++i) {
    size_t r = (rand()>>3)%(n-i), j;
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
