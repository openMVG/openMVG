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

#ifndef OPENMVG_MULTIVIEW_CONDITIONNING_HPP
#define OPENMVG_MULTIVIEW_CONDITIONNING_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

//-- Implementation of normalized coordinates.
// Normalization improve accuracy of results and provide benefits
//  that make scale and coordinate origin invariant.
// The implementation follows Algorithm 4.2 from HZ page 109.

namespace openMVG
{

/**
* @brief Compute conditioning matrix from input points
* @param points Input points
* @param[out] T Output conditioning matrix
*/
void PreconditionerFromPoints( const Mat &points, Mat3 *T );

/**
* @brief Normalize input point for a given T transform matrix
* @param points Input points to normalize
* @param T Input conditioning matrix
* @param[out] transformed_points transformed (i.e.: conditioned ) points
*/
void ApplyTransformationToPoints( const Mat &points,
                                  const Mat3 &T,
                                  Mat *transformed_points );

/**
* @brief Normalize point in [-.5, .5] and return transformation matrix
* @param points Input points
* @param[out] normalized_points Points after conditioning
* @param[out] T Conditioning matrix used to normalize input points
*/
void NormalizePoints( const Mat &points,
                      Mat *normalized_points,
                      Mat3 *T );


/**
* @brief Compute conditioning from a 2d range
* @param width First range upper bound
* @param height Second range upper bound
* @param[out] Transformation matrix
* @note Transformation compress input range to [ -sqrt(2); sqrt(2)]
* @note Range is [0;width]x[0;height]
*/
void PreconditionerFromPoints( int width, int height, Mat3 *T );

/**
* @brief Normalize point row image coordinates to [- sqrt(2); sqrt(2) ]
* @param points Input points
* @param[out] normalized_points Normalized points
* @param[out] T Normalization matrix used to normalize points
* @param width Range of points
* @param height Range of points
*/
void NormalizePoints( const Mat &points,
                      Mat *normalized_points,
                      Mat3 *T, int width, int height );


/**
* @brief Unnormalize using Inverse
*/
struct UnnormalizerI
{
  /**
  * @brief Denormalize the results.
  * @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
  * @see HZ page 109, H = T2^{-1} H T1
  * @param T1 Input transformation of first dataset
  * @param T2 Input transformation of second dataset
  * @param H Denormalization transformation
  */
  static void Unnormalize( const Mat3 &T1, const Mat3 &T2, Mat3 *H );
};

/**
* Unnormalize using Transpose
*/
struct UnnormalizerT
{
  /**
  * @brief Denormalize the results.
  * @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman
  * @see HZ page 282, H = T2^T H T1
  * @param T1 Input transformation of first dataset
  * @param T2 Input transformation of second dataset
  * @param H Denormalization transformation
  */
  static void Unnormalize( const Mat3 &T1, const Mat3 &T2, Mat3 *H );
};

} //namespace openMVG


#endif // OPENMVG_MULTIVIEW_CONDITIONNING_HPP
