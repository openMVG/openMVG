
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

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_H
#define OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_H

#include "openMVG/numeric/numeric.h"

namespace openMVG {

  /// Compute a 3D position of a point from several images of it. In particular,
  ///  compute the projective point X in R^4 such that x = PX.
  /// Algorithm is the standard DLT; for derivation see appendix of Keir's thesis.
  void TriangulateNView(
    const Mat2X &x, // x's are 2D coordinates (x,y,1) in each image
    const std::vector< Mat34 > &Ps, // Ps are projective cameras
    Vec4 *X);

  // This method uses the algebraic distance approximation.
  // Note that this method works better when the 2D points are normalized
  // with an isotopic normalization.
  void TriangulateNViewAlgebraic(
    const Mat2X &x, // x's are 2D coordinates (x,y,1) in each image
    const std::vector< Mat34 > &Ps, // Ps are projective cameras.
    Vec4 *X);

  //Iterated linear method
  class Triangulation
	{
	public:

		size_t size() const {	return views.size();}

		void clear()  { views.clear();}

		void add(const Mat34& projMatrix, const Vec2 & p) {
			views.push_back( std::pair<Mat34, Vec2>(projMatrix,p) );
		}

    // Return squared L2 sum of error
		double error(const Vec3 &X) const;

		Vec3 compute(int iter = 3) const;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Accessors

		// These values are defined after a successful call to compute
		double minDepth() const { return zmin; }
		double maxDepth() const { return zmax; }
		double error()    const { return err; }

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Data members

	protected:
		mutable double zmin; // min depth, mutable since modified in compute(...) const;
    mutable double zmax; // max depth, mutable since modified in compute(...) const;
    mutable double err;  // re-projection error, mutable since modified in compute(...) const;
		std::vector< std::pair<Mat34, Vec2> > views; // Proj matrix and associated image point
	};

}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_H
