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

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_CONDITIONNING_H_
#define OPENMVG_MULTIVIEW_CONDITIONNING_H_

#include "openMVG/numeric/numeric.h"

//-- Implementation of normalized coordinates.
// Normalization improve accuracy of results and provide benefits
//  that make scale and coordinate origin invariant.
// The implementation follows Algorithm 4.2 from HZ page 109.

namespace openMVG {

// Point conditioning :
void PreconditionerFromPoints(const Mat &points, Mat3 *T);

/// Normalize input point for a given T transform matrix
void ApplyTransformationToPoints(const Mat &points,
                                 const Mat3 &T,
                                 Mat *transformed_points);

// Normalize point in [-.5, .5] and return transformation matrix
void NormalizePoints(const Mat &points,
                     Mat *normalized_points,
                     Mat3 *T);

/// Point conditioning (compute Transformation matrix)
void PreconditionerFromPoints(int width, int height, Mat3 *T);

///  Normalize point rom image coordinates to [-.5, .5]
void NormalizePoints(const Mat &points,
                     Mat *normalized_points,
                     Mat3 *T, int width, int height);


/// Unnormalize using Inverse
struct UnnormalizerI {
  // Denormalize the results. See HZ page 109.
  static void Unnormalize(const Mat3 &T1, const Mat3 &T2, Mat3 *H);
};

/// Unnormalize using Transpose
struct UnnormalizerT {
  // Denormalize the results. See HZ page 109.
  static void Unnormalize(const Mat3 &T1, const Mat3 &T2, Mat3 *H);
};

} //namespace openMVG


#endif // OPENMVG_MULTIVIEW_CONDITIONNING_H_
