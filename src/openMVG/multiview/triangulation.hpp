
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

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_H_
#define OPENMVG_MULTIVIEW_TRIANGULATION_H_

#include "openMVG/numeric/numeric.h"

namespace openMVG
{

/**
* @brief Linear DLT triangulation
* @brief P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 Point in first camera
* @param x2 Point in second camera
* @param[out] X_homogeneous Homogeneous triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT( const Mat34 &P1, const Vec2 &x1,
                     const Mat34 &P2, const Vec2 &x2,
                     Vec4 *X_homogeneous );

/**
* @brief Linear DLT triangulation
* @brief P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 Point in first camera
* @param x2 Point in second camera
* @param[out] X_euclidean Euclidean triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT( const Mat34 &P1, const Vec2 &x1,
                     const Mat34 &P2, const Vec2 &x2,
                     Vec3 *X_euclidean );

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_H_
