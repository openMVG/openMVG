
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

#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/multiview/projection.hpp"

namespace openMVG {

  void TriangulateNView(const Mat2X &x,
    const std::vector< Mat34 > &Ps,
    Vec4 *X) {
      Mat2X::Index nviews = x.cols();
      assert(nviews == Ps.size());

      Mat design = Mat::Zero(3*nviews, 4 + nviews);
      for (int i = 0; i < nviews; i++) {
        design.block<3, 4>(3*i, 0) = -Ps[i];
        design(3*i + 0, 4 + i) = x(0, i);
        design(3*i + 1, 4 + i) = x(1, i);
        design(3*i + 2, 4 + i) = 1.0;
      }
      Vec X_and_alphas;
      Nullspace(&design, &X_and_alphas);
      *X = X_and_alphas.head(4);
  }

  typedef Eigen::Matrix<double, 2, 3> Mat23;
  inline Mat23 SkewMatMinimal(const Vec2 &x) {
    Mat23 skew;
    skew <<
      0, -1,  x(1),
      1, 0, -x(0);
    return skew;
  }

  void TriangulateNViewAlgebraic(const Mat2X &x,
    const std::vector< Mat34 > &Ps,
    Vec4 *X) {
      Mat2X::Index nviews = x.cols();
      assert(nviews == Ps.size());

      Mat design(2*nviews, 4);
      for (int i = 0; i < nviews; i++) {
        design.block<2, 4>(2*i, 0) = SkewMatMinimal(x.col(i)) * Ps[i];
      }
      Nullspace(&design, X);
  }

	double Triangulation::error(const Vec3 &X) const
	{
		double squared_reproj_error = 0.0;
		for (std::size_t i=0;i<views.size();i++)
		{
			const Mat34& PMat = views[i].first;
			const Vec2 & xy = views[i].second;
			const Vec2 p = Project(PMat, X);
			squared_reproj_error += (xy-p).norm();
		}
		return squared_reproj_error;
	}

	// Camera triangulation using the iterated linear method
	Vec3 Triangulation::compute(int iter) const
	{
		const int nviews = int(views.size());
		assert(nviews>=2);

		// Iterative weighted linear least squares
		Mat3 AtA;
		Vec3 Atb, X;
		std::vector<double> weights(nviews,double(1.0));
		for (int it=0;it<iter;++it)
		{
			AtA.fill(0.0);
			Atb.fill(0.0);
			for (int i=0;i<nviews;++i)
			{
				const Mat34& PMat = views[i].first;
				const Vec2 & p = views[i].second;
				const double w = weights[i];

				Vec3 v1,v2;
				for (int j=0;j<3;j++)
				{
					v1[j] = w * ( PMat(0,j) - p(0) * PMat(2,j) );
					v2[j] = w * ( PMat(1,j) - p(1) * PMat(2,j) );
					Atb[j] += w * ( v1[j] * ( p(0) * PMat(2,3) - PMat(0,3) )
            + v2[j] * ( p(1) * PMat(2,3) - PMat(1,3) ) );
				}

				for (int k=0;k<3;k++)
				{
					for (int j=0;j<=k;j++)
					{
						const double v = v1[j] * v1[k] + v2[j] * v2[k];
						AtA(j,k) += v;
						if (j<k) AtA(k,j) += v;
					}
				}
			}

			X = AtA.inverse() * Atb;

			// Compute reprojection error, min and max depth, and update weights
			zmin = std::numeric_limits<double>::max();
			zmax = - std::numeric_limits<double>::max();
			err = 0;
			for (int i=0;i<nviews;++i)
			{
				const Mat34& PMat = views[i].first;
				const Vec2 & p = views[i].second;
				const Vec3 xProj = PMat * Vec4(X(0), X(1), X(2), 1.0);
				const double z = xProj(2);
        const Vec2 x = xProj.head<2>() / z;
				if (z<zmin) zmin = z;
				if (z>zmax) zmax = z;
				err += (p-x).norm();
				weights[i] = 1.0 / z;
			}
		}
		return X;
	}


}  // namespace openMVG

