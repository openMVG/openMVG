// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SPHERICAL_CAMERA_HPP
#define OPENMVG_SPHERICAL_CAMERA_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"

// [1] "Robust and accurate calibration of camera networks". PhD.
// Authors: Pierre MOULON
// Chapter 4.4.2 Application to estimation of relative orientation between two
//  spherical camera.

// Implementation of the AContrario estimation of the relative pose estimation
//  between two spherical panorama.

namespace openMVG {
namespace spherical_cam {

static void planarToSpherical(
  const Mat & planarCoords, // Input (x,y)' coords
  size_t width,             // Width of the 2D planar surface
  size_t height,            // Height of the 2D planar surface
  Mat & sphericalCoords)     // Output spherical coordinates (on the unit sphere)
{
  sphericalCoords.resize(3, planarCoords.cols());
  for (size_t iCol = 0; iCol < planarCoords.cols(); ++iCol)
  {
    const Vec2 & xy = planarCoords.col(iCol);
    double uval = xy(0) / (double)width;
    double vval = xy(1) / (double)height;

    sphericalCoords.col(iCol) =
      (Vec3 (sin(vval*M_PI)*cos(M_PI*(2.0*uval+0.5)),
            cos(vval*M_PI),
            sin(vval*M_PI)*sin(M_PI*(2.0*uval+0.5)))).normalized();
  }
}

using namespace std;

/**
 * Eight-point algorithm for solving for the essential matrix from normalized
 * image coordinates of point correspondences.
 * See page 294 in HZ Result 11.1.
 *
 */
struct EightPointRelativePoseSolver {
  enum { MINIMUM_SAMPLES = 8 };
  enum { MAX_MODELS = 1 };
  static void Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *pvec_E)
  {
    assert(3 == x1.rows());
    assert(8 <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());

    MatX9 A(x1.cols(), 9);
    EncodeEpipolarEquation(x1, x2, &A);

    Vec9 e;
    Nullspace(&A, &e);
    Mat3 E = Map<RMat3>(e.data());

    // Find the closest essential matrix to E in frobenius norm
    // E = UD'VT
    if (x1.cols() > 8) {
      Eigen::JacobiSVD<Mat3> USV(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Vec3 d = USV.singularValues();
      double a = d[0];
      double b = d[1];
      d << (a+b)/2., (a+b)/2., 0.0;
      E = USV.matrixU() * d.asDiagonal() * USV.matrixV().transpose();
    }
    pvec_E->push_back(E);

  }

  template<typename TMatX, typename TMatA>
  static inline void EncodeEpipolarEquation(const TMatX &x1, const TMatX &x2, TMatA *A) {
    for (int i = 0; i < x1.cols(); ++i) {
      (*A)(i, 0) = x2(0, i) * x1(0, i);  // 0 represents x coords,
      (*A)(i, 1) = x2(0, i) * x1(1, i);  // 1 represents y coords,
      (*A)(i, 2) = x2(0, i) * x1(2, i);  // 2 represents z coords.
      (*A)(i, 3) = x2(1, i) * x1(0, i);
      (*A)(i, 4) = x2(1, i) * x1(1, i);
      (*A)(i, 5) = x2(1, i) * x1(2, i);
      (*A)(i, 6) = x2(2, i) * x1(0, i);
      (*A)(i, 7) = x2(2, i) * x1(1, i);
      (*A)(i, 8) = x2(2, i) * x1(2, i);
    }
  }
};

// Return the angular error between [0; PI/2]
struct AngularError {
  static double Error(const Mat3 &model, const Vec3 &x1, const Vec3 &x2) {
    const Vec3 Em1 = (model * x1).normalized();
    double angleVal = (x2.transpose() * Em1);
    angleVal /= (x2.norm() * Em1.norm());
    return abs(asin(angleVal));
  }
};

class EssentialKernel_spherical
{
public:
  typedef Mat3 Model;
  enum { MINIMUM_SAMPLES = EightPointRelativePoseSolver::MINIMUM_SAMPLES };

  EssentialKernel_spherical(const Mat &x1, const Mat &x2) : x1_(x1), x2_(x2) {}

  void Fit(const vector<size_t> &samples, std::vector<Model> *models) const {
    const Mat x1 = ExtractColumns(x1_, samples);
    const Mat x2 = ExtractColumns(x2_, samples);

    assert(3 == x1.rows());
    assert(MINIMUM_SAMPLES <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());

    EightPointRelativePoseSolver::Solve(x1, x2, models);
  }

  size_t NumSamples() const {return x1_.cols();}

  /// Return the angular error (between 0 and PI/2)
  double Error(size_t sample, const Model &model) const {
    return AngularError::Error(model, x1_.col(sample), x2_.col(sample));
  }

  protected:
  const Mat & x1_, & x2_;
};

// Solve:
// [cross(x0,P0) X = 0]
// [cross(x1,P1) X = 0]
void TriangulateDLT(const Mat34 &P1, const Vec3 &x1,
                    const Mat34 &P2, const Vec3 &x2,
                    Vec4 *X_homogeneous) {
  Mat design(6,4);
  for (int i = 0; i < 4; ++i) {
    design(0,i) = -x1[2] * P1(1,i) + x1[1] * P1(2,i);
    design(1,i) =  x1[2] * P1(0,i) - x1[0] * P1(2,i);
    design(2,i) = -x1[1] * P1(0,i) + x1[0] * P1(1,i);

    design(3,i) = -x2[2] * P2(1,i) + x2[1] * P2(2,i);
    design(4,i) =  x2[2] * P2(0,i) - x2[0] * P2(2,i);
    design(5,i) = -x2[1] * P2(0,i) + x2[0] * P2(1,i);
  }
  Nullspace(&design, X_homogeneous);
}

void TriangulateDLT(const Mat34 &P1, const Vec3 &x1,
                    const Mat34 &P2, const Vec3 &x2,
                    Vec3 *X_euclidean) {
  Vec4 X_homogeneous;
  TriangulateDLT(P1, x1, P2, x2, &X_homogeneous);
  HomogeneousToEuclidean(X_homogeneous, X_euclidean);
}

} // namespace spherical_cam
} // namespace openMVG

#endif /* OPENMVG_SPHERICAL_CAMERA_HPP */
