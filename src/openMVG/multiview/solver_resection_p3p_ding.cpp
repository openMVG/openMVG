// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2024 Yaqing Ding
// Some of the scripts are based on Mark Shachkov (mark.shachkov@gmail.com) and
// the Lambda-twist P3P implementation Adapted to openMVG by Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p3p_ding.hpp"
#include "openMVG/multiview/projection.hpp"

#include <array>

namespace openMVG {
namespace euclidean_resection {

/**
 * @brief Refine a valid solution with a Gauss-Newton Solver.
 * @param lambda Valid lambdas
 * @param a12 is the squared distance between X1 and X2
 * @param a13 is the squared distance between X1 and X3
 * @param a23 is the squared distance between X2 and X3
 * @param b12 is the cosine of the angle between bearing vector 1 and bearing
 * vector 2
 * @param b13 is the cosine of the angle between bearing vector 1 and bearing
 * vector 3
 * @param b23 is the cosine of the angle between bearing vector 2 and bearing
 * vector 3 The paper note it rarely improve after two iterations. The original
 * implementation use 5 iterations.
 */
static void gauss_newton_refineL(double &lambda1, double &lambda2,
                                 double &lambda3, const double &a12,
                                 const double &a13, const double &a23,
                                 const double &b12, const double &b13,
                                 const double &b23) {
  for (int iter = 0; iter < 5; ++iter) {
    double r1 = (lambda1 * lambda1 - 2.0 * lambda1 * lambda2 * b12 +
                 lambda2 * lambda2 - a12);
    double r2 = (lambda1 * lambda1 - 2.0 * lambda1 * lambda3 * b13 +
                 lambda3 * lambda3 - a13);
    double r3 = (lambda2 * lambda2 - 2.0 * lambda2 * lambda3 * b23 +
                 lambda3 * lambda3 - a23);
    if (std::abs(r1) + std::abs(r2) + std::abs(r3) < 1e-10)
      return;
    double x11 = lambda1 - lambda2 * b12;
    double x12 = lambda2 - lambda1 * b12;
    double x21 = lambda1 - lambda3 * b13;
    double x23 = lambda3 - lambda1 * b13;
    double x32 = lambda2 - lambda3 * b23;
    double x33 = lambda3 - lambda2 * b23;
    double detJ = 0.5 / (x11 * x23 * x32 +
                         x12 * x21 * x33); // half minus inverse determinant
    // This uses the closed form of the inverse for the jacobian.
    // Due to the zero elements this actually becomes quite nice.
    lambda1 += (-x23 * x32 * r1 - x12 * x33 * r2 + x12 * x23 * r3) * detJ;
    lambda2 += (-x21 * x33 * r1 + x11 * x33 * r2 - x11 * x23 * r3) * detJ;
    lambda3 += (x21 * x32 * r1 - x11 * x32 * r2 - x12 * x21 * r3) * detJ;
  }
};

/**
 * @brief This function finds a single root of the cubic polynomial equation
 * @param c2 Coefficient of quadratic parameter
 * @param c1 Coefficient of linear parameter
 * @param c0 Coefficient of scalar parameter
 * @param root the single root output
 * @return true: if the Discriminant is positive
 *
 * h(r) = r^3 + c2*r^2 + c1*r + c0 = 0
 *
 */
bool solve_cubic_single_real(double c2, double c1, double c0, double &root) {
  double a = c1 - c2 * c2 / 3.0;
  double b = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1) / 27.0 + c0;
  double c = b * b / 4.0 + a * a * a / 27.0;
  if (c != 0) {
    if (c > 0) {
      c = std::sqrt(c);
      b *= -0.5;
      root = std::cbrt(b + c) + std::cbrt(b - c) - c2 / 3.0;
      return true;
    } else {
      c = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
      root =
          2.0 * std::sqrt(-a / 3.0) * std::cos(std::acos(c) / 3.0) - c2 / 3.0;
    }
  } else {
    root = -c2 / 3.0 + (a != 0 ? (3.0 * b / a) : 0);
  }
  return false;
}
/**
 * @brief This function decomposes a degenerate conic into a pair of lines
 * @param C The degenerate conic, which can be written as C = p*q^T + q*p^T
 * @return A pair of lines
 */
inline std::array<Vec3, 2> compute_pq(Mat3 C) {
  std::array<Vec3, 2> pq;
  Mat3 C_adj;

  C_adj(0, 0) = C(1, 2) * C(2, 1) - C(1, 1) * C(2, 2);
  C_adj(1, 1) = C(0, 2) * C(2, 0) - C(0, 0) * C(2, 2);
  C_adj(2, 2) = C(0, 1) * C(1, 0) - C(0, 0) * C(1, 1);
  C_adj(0, 1) = C(0, 1) * C(2, 2) - C(0, 2) * C(2, 1);
  C_adj(0, 2) = C(0, 2) * C(1, 1) - C(0, 1) * C(1, 2);
  C_adj(1, 0) = C_adj(0, 1);
  C_adj(1, 2) = C(0, 0) * C(1, 2) - C(0, 2) * C(1, 0);
  C_adj(2, 0) = C_adj(0, 2);
  C_adj(2, 1) = C_adj(1, 2);

  Vec3 v;
  if (C_adj(0, 0) > C_adj(1, 1)) {
    if (C_adj(0, 0) > C_adj(2, 2)) {
      v = C_adj.col(0) / std::sqrt(C_adj(0, 0));
    } else {
      v = C_adj.col(2) / std::sqrt(C_adj(2, 2));
    }
  } else if (C_adj(1, 1) > C_adj(2, 2)) {
    v = C_adj.col(1) / std::sqrt(C_adj(1, 1));
  } else {
    v = C_adj.col(2) / std::sqrt(C_adj(2, 2));
  }

  C(0, 1) -= v(2);
  C(0, 2) += v(1);
  C(1, 2) -= v(0);
  C(1, 0) += v(2);
  C(2, 0) -= v(1);
  C(2, 1) += v(0);

  pq[0] = C.col(0);
  pq[1] = C.row(0);

  return pq;
}

/**
 * @brief Compute the absolute pose of a camera using three 3D-to-2D
 * correspondences. Implementation of the paper "Revisiting the P3P Problem". Y
 * Ding, J Yang, V Larsson, C Olsson, K Åström. CVPR 2023
 *
 * @authors Y Ding, J Yang, V Larsson, C Olsson, K Åström.
 *
 * @param[in] bearing_vectors 3x3 matrix with UNITARY feature vectors (each
 * column is a vector)
 * @param[in] X  3x3 matrix with corresponding 3D world points (each column is a
 * point)
 * @param[out] rotation_translation_solutions vector that will contain the
 * solutions (up to 4 solutions)
 *
 * @return true if at least one solution is found, false if no solution was
 * found
 *
 */
bool computePosesDing(
    const Mat &bearing_vectors, const Mat &X,
    std::vector<std::tuple<Mat3, Vec3>> &rotation_translation_solutions) {
  // Extraction of 3D points vectors

  Vec3 X0 = X.col(0);
  Vec3 X1 = X.col(1);
  Vec3 X2 = X.col(2);

  Vec3 X01 = X0 - X1;
  Vec3 X02 = X0 - X2;
  Vec3 X12 = X1 - X2;

  double a01 = X01.squaredNorm();
  double a02 = X02.squaredNorm();
  double a12 = X12.squaredNorm();

  // Extraction of feature vectors
  Vec3 x0 = bearing_vectors.col(0);
  Vec3 x1 = bearing_vectors.col(1);
  Vec3 x2 = bearing_vectors.col(2);

  // Switch the columns of X and x so that BC = |X1 - X2| is the largest
  // one among {X01, X02, X12}
  if (a01 > a02) {
    if (a01 > a12) {
      std::swap(X0, X2);
      std::swap(x0, x2);
      std::swap(a01, a12);
      X01 = -X12;
      X02 = -X02;
    }
  } else if (a02 > a12) {
    std::swap(X0, X1);
    std::swap(x0, x1);
    std::swap(a02, a12);
    X01 = -X01;
    X02 = X12;
  }

  const double a12d = 1.0 / a12;
  const double a = a01 * a12d;
  const double b = a02 * a12d;

  const double m01 = x0.dot(x1);
  const double m02 = x0.dot(x2);
  const double m12 = x1.dot(x2);

  // Ugly parameters to simplify the calculation
  const double m12sq = -m12 * m12 + 1.0;
  const double m02sq = -1.0 + m02 * m02;
  const double m01sq = -1.0 + m01 * m01;
  const double ab = a * b;
  const double bsq = b * b;
  const double asq = a * a;
  const double m013 = -2.0 + 2.0 * m01 * m02 * m12;
  const double bsqm12sq = bsq * m12sq;
  const double asqm12sq = asq * m12sq;
  const double abm12sq = 2.0 * ab * m12sq;
  // coefficients of the cubic equation
  const double k3_inv = 1.0 / (bsqm12sq + b * m02sq);
  const double k2 =
      k3_inv * ((-1.0 + a) * m02sq + abm12sq + bsqm12sq + b * m013);
  const double k1 =
      k3_inv * (asqm12sq + abm12sq + a * m013 + (-1.0 + b) * m01sq);
  const double k0 = k3_inv * (asqm12sq + a * m01sq);

  double s;
  bool G = solve_cubic_single_real(k2, k1, k0, s);
  // The degenerate conic C = C1+s*C2
  Mat3 C;
  C(0, 0) = -a + s * (1 - b);
  C(0, 1) = -m02 * s;
  C(0, 2) = a * m12 + b * m12 * s;
  C(1, 0) = C(0, 1);
  C(1, 1) = s + 1;
  C(1, 2) = -m01;
  C(2, 0) = C(0, 2);
  C(2, 1) = C(1, 2);
  C(2, 2) = -a - b * s + 1;

  std::array<Vec3, 2> pq = compute_pq(C);

  Mat3 XX;

  XX << X01, X02, X01.cross(X02);
  XX = XX.inverse().eval();

  Vec3 v1, v2;
  Mat3 YY;
  int n_sols = 0;
  double d0, d1, d2;

  for (int i = 0; i < 2; ++i) {
    // [p0 p1 p2] * [1; x; y] = 0, or [p0 p1 p2] * [d2; d0; d1] = 0
    double p0 = pq[i](0);
    double p1 = pq[i](1);
    double p2 = pq[i](2);
    bool switch_12 = std::abs(p0) <= std::abs(p1);
    if (switch_12) {
      // eliminate d0
      double w0 = -p0 / p1;
      double w1 = -p2 / p1;
      double ca = 1.0 / (w1 * w1 - b);
      double cb = 2.0 * (b * m12 - m02 * w1 + w0 * w1) * ca;
      double cc = (w0 * w0 - 2 * m02 * w0 - b + 1.0) * ca;
      double taus[2];
      if (!P3PSolver_Nordberg::root2real(cb, cc, taus[0], taus[1]))
        continue;
      for (double tau : taus) {
        if (tau <= 0)
          continue;
        // positive only
        d2 = std::sqrt(a12 / (tau * (tau - 2.0 * m12) + 1.0));
        d1 = tau * d2;
        d0 = (w0 * d2 + w1 * d1);
        if (d0 < 0)
          continue;

        // refine the depths
        gauss_newton_refineL(d0, d1, d2, a01, a02, a12, m01, m02, m12);

        v1 = d0 * x0 - d1 * x1;
        v2 = d0 * x0 - d2 * x2;
        YY << v1, v2, v1.cross(v2);
        Mat3 R = YY * XX;

        rotation_translation_solutions.emplace_back(R, d0 * x0 - R * X0);

        ++n_sols;
      }
    } else {
      double w0 = -p1 / p0;
      double w1 = -p2 / p0;
      double ca = 1.0 / (-a * w1 * w1 + 2 * a * m12 * w1 - a + 1);
      double cb = 2 * (a * m12 * w0 - m01 - a * w0 * w1) * ca;
      double cc = (1 - a * w0 * w0) * ca;

      double taus[2];
      if (!P3PSolver_Nordberg::root2real(cb, cc, taus[0], taus[1]))
        continue;
      for (double tau : taus) {
        if (tau <= 0)
          continue;
        d0 = std::sqrt(a01 / (tau * (tau - 2.0 * m01) + 1.0));
        d1 = tau * d0;
        d2 = w0 * d0 + w1 * d1;

        if (d2 < 0)
          continue;

        gauss_newton_refineL(d0, d1, d2, a01, a02, a12, m01, m02, m12);
        v1 = d0 * x0 - d1 * x1;
        v2 = d0 * x0 - d2 * x2;
        YY << v1, v2, v1.cross(v2);
        Mat3 R = YY * XX;
        rotation_translation_solutions.emplace_back(R, d0 * x0 - R * X0);
        ++n_sols;
      }
    }
    // if n_sols > 0, we have found at least one solution. If G > 0, It means
    // that only one line from the pair of lines has intersections with the
    // conic. Hence, we can skip the second line.
    if (n_sols > 0 && G)
      break;
  }

  return n_sols > 0;
}

void P3PSolver_Ding::Solve(const Mat &bearing_vectors,
                           const Mat &X, // 3D points
                           std::vector<Mat34> *models) {
  assert(3 == bearing_vectors.rows());
  assert(3 == X.rows());
  assert(bearing_vectors.cols() == X.cols());
  std::vector<std::tuple<Mat3, Vec3>> rotation_translation_solutions;
  if (computePosesDing(bearing_vectors, X, rotation_translation_solutions)) {
    for (const auto &rot_trans_it : rotation_translation_solutions) {
      Mat34 P;
      P_From_KRt(Mat3::Identity(),          // intrinsics
                 std::get<0>(rot_trans_it), // rotation
                 std::get<1>(rot_trans_it), // translation
                 &P);
      models->push_back(P);
    }
  }
};

} // namespace euclidean_resection
} // namespace openMVG
