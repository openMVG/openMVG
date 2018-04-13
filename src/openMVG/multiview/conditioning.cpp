// Copyright (c) 2010 libmv authors.
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

#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {

// HZ 4.4.4 pag.109
void PreconditionerFromPoints(const Mat &points, Mat3 *T) {
  Vec mean, variance;
  MeanAndVarianceAlongRows(points, &mean, &variance);

  double xfactor = sqrt(2.0 / variance(0));
  double yfactor = sqrt(2.0 / variance(1));

  // If variance is equal to 0.0 set scaling factor to identity.
  // -> Else it will provide nan value (because division by 0).
  if (variance(0) < 1e-8)
    xfactor = mean(0) = 1.0;
  if (variance(1) < 1e-8)
    yfactor = mean(1) = 1.0;

  (*T) << xfactor, 0,       -xfactor * mean(0),
          0,       yfactor, -yfactor * mean(1),
          0,       0,        1;
}

void PreconditionerFromPoints(int width, int height, Mat3 *T) {
  // Build the normalization matrix
  const double dNorm = 1.0 / sqrt( static_cast<double>(width*height) );

  (*T) = Mat3::Identity();
  (*T)(0,0) = (*T)(1,1) = dNorm;
  (*T)(0,2) = -.5f*width*dNorm;
  (*T)(1,2) = -.5*height*dNorm;
}

void ApplyTransformationToPoints(const Mat &points,
                                 const Mat3 &T,
                                 Mat *transformed_points) {
  (*transformed_points) = (T * points.colwise().homogeneous()).colwise().hnormalized();
}

void NormalizePoints(const Mat &points,
                      Mat *normalized_points,
                      Mat3 *T,
                      int width,
                      int height) {
  PreconditionerFromPoints(width, height, T);
  ApplyTransformationToPoints(points, *T, normalized_points);
}


void NormalizePoints(const Mat &points,
                     Mat *normalized_points,
                     Mat3 *T) {
  PreconditionerFromPoints(points, T);
  ApplyTransformationToPoints(points, *T, normalized_points);
}

// Denormalize the results. See HZ page 282.
void UnnormalizerT::Unnormalize(const Mat3 &T1, const Mat3 &T2, Mat3 *H)  {
  *H = T2.transpose() * (*H) * T1;
}

// Denormalize the results. See HZ page 109.
void UnnormalizerI::Unnormalize(const Mat3 &T1, const Mat3 &T2, Mat3 *H)  {
  *H = T2.inverse() * (*H) * T1;
}

} // namespace openMVG
