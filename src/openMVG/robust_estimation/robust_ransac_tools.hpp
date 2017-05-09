
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

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_RANSAC_TOOLS_HPP
#define OPENMVG_ROBUST_RANSAC_TOOLS_HPP

#include <cmath>

namespace openMVG {
namespace robust{

/// Number of samplings to have at least \a minProba probability of absence of
/// outlier in a sample of \a SampleSize elements.
inline size_t getNumSamples(
  double minProba,
  double outlierRatio,
  std::size_t SampleSize)
{
  return static_cast<std::size_t>( std::log(1.-minProba) /
    std::log(1.-std::pow(1.-outlierRatio, static_cast<int>(SampleSize))));
}

inline size_t IterationsRequired(
  std::size_t min_samples,
  double outliers_probability,
  double inlier_ratio)
{
  return static_cast<std::size_t>(
    std::log(outliers_probability) /
    std::log(1.0 - std::pow(inlier_ratio, static_cast<int>(min_samples))));
}

} // namespace robust
} // namespace openMVG

#endif // OPENMVG_ROBUST_RANSAC_TOOLS_HPP
