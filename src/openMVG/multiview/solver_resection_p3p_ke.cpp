// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p3p_ke.hpp"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/numeric/poly.h"

#include <array>
#include <cmath>
#include <complex>
#include <tuple>

using namespace openMVG;

namespace openMVG {
namespace euclidean_resection {

/**
* @brief Compute the absolute pose of a camera using three 3D-to-2D correspondences.
*  Implementation of the paper "An Efficient Algebraic Solution to the Perspective-Three-Point Problem". Tong Ke, Stergios Roumeliotis. CVPR 2017
*
* @authors Tong Ke, Stergios Roumeliotis (original code). Pierre Moulon (adaptation to openMVG).
*
* @param[in] bearing_vectors 3x3 matrix with UNITARY feature vectors (each column is a vector)
* @param[in] X_observations  3x3 matrix with corresponding 3D world points (each column is a point)
* @param[in] rotation_translation_solutions vector that will contain the solutions (up to 4 solutions)
*
* @return true if at least one solution is found, false if no solution was found
*
*/
bool computePoses
(
  const Mat & bearing_vectors,
  const Mat & X_observations,
  std::vector<std::tuple<Mat3, Vec3>> & rotation_translation_solutions
)
{
  //world point vectors
  const Vec3 w1 = X_observations.col(0);
  const Vec3 w2 = X_observations.col(1);
  const Vec3 w3 = X_observations.col(2);

  // k1
  const Vec3 u0 = w1 - w2;
  const double nu0 = u0.norm();
  const Vec3 k1 = u0.normalized();

  // bi
  const Vec3 b1 = bearing_vectors.col(0);
  const Vec3 b2 = bearing_vectors.col(1);
  const Vec3 b3 = bearing_vectors.col(2);
  // k3, tz
  Vec3 k3 = b1.cross(b2);
  const double nk3 = k3.norm();
  k3 = k3.normalized();

  const Vec3 tz = b1.cross(k3);

  // ui, vi
  const Vec3 v1 = b1.cross(b3);
  const Vec3 v2 = b2.cross(b3);

  const Vec3 u1 = w1 - w3;
  // coefficients related terms
  const double u1k1 = u1.dot(k1);
  const double k3b3 = k3.dot(b3);
  // f1i
  double f11 = k3b3;
  double f13 = k3.dot(v1);
  const double f15 = -u1k1 * f11;
  //delta
  const Vec3 nl = u1.cross(k1).normalized();
  const double delta = u1.cross(k1).norm();
  f11 *= delta;
  f13 *= delta;
  // f2i
  const double u2k1 = u1k1 - nu0;
  double f21 = tz.dot(v2);
  double f22 = nk3 * k3b3;
  double f23 = k3.dot(v2);
  const double f24 = u2k1 * f22;
  const double f25 = -u2k1 * f21;
  f21 *= delta;
  f22 *= delta;
  f23 *= delta;
  const double g1 = f13 * f22;
  const double g2 = f13 * f25 - f15 * f23;
  const double g3 = f11 * f23 - f13 * f21;
  const double g4 = -f13 * f24;
  const double g5 = f11 * f22;
  const double g6 = f11 * f25 - f15 * f21;
  const double g7 = -f15 * f24;
  const std::array<double, 5> coeffs = {
    {g5 * g5 + g1 * g1 + g3 * g3,
      2 * (g5 * g6 + g1 * g2 + g3 * g4),
      g6 * g6 + 2 * g5 * g7 + g2 * g2 + g4 * g4 - g1 * g1 - g3 * g3,
      2 * (g6 * g7 - g1 * g2 - g3 * g4),
      g7 * g7 - g2 * g2 - g4 * g4}
  };
  std::array<double, 4> s;
  solveQuarticPolynomial(coeffs, s);
  polishQuarticPolynomialRoots(coeffs, s);

  const Vec3 temp = k1.cross(nl);

  Mat3 Ck1nl;
  Ck1nl << k1, nl, temp;

  Mat3 Cb1k3tzT;
  Cb1k3tzT << b1.transpose(), k3.transpose(), tz.transpose();

  const Vec3 b3p = b3 * (delta / k3b3);

  for (const auto ctheta1p : s) {
    if (std::abs(ctheta1p) > 1)
      continue;
    const double stheta1p = ((k3b3 > 0) ? 1 : -1) * sqrt(1 - ctheta1p * ctheta1p);
    const double ntheta3 = stheta1p / ((g5 * ctheta1p + g6) * ctheta1p + g7);
    const double ctheta3 = (g1 * ctheta1p + g2) * ntheta3;
    const double stheta3 = (g3 * ctheta1p + g4) * ntheta3;

    Mat3 C13;
    C13 <<
      ctheta3,            0,         -stheta3,
      stheta1p * stheta3, ctheta1p,  stheta1p * ctheta3,
      ctheta1p * stheta3, -stheta1p, ctheta1p * ctheta3;

    const Mat3 R = (Ck1nl * C13) * Cb1k3tzT;
    const Vec3 rp3 = R.transpose() * w3; // R' * p3
    rotation_translation_solutions.emplace_back(R.transpose(), (b3p * stheta1p) - rp3);
  }

  return !rotation_translation_solutions.empty();
}

void P3PSolver_Ke::Solve
(
  const Mat & bearing_vectors,
  const Mat & X,
  std::vector<Mat34> * models
)
{
  std::vector<std::tuple<Mat3, Vec3>> rotation_translation_solutions;
  if (computePoses(bearing_vectors, X, rotation_translation_solutions))
  {
    for (const auto & rot_trans_it : rotation_translation_solutions) {
      Mat34 P;
      P_From_KRt( Mat3::Identity(),          // intrinsics
                  std::get<0>(rot_trans_it), // rotation
                  std::get<1>(rot_trans_it), // translation
                  &P);
      models->push_back(P);
    }
  }
}

double P3PSolver_Ke::Error
(
  const Mat34 & P,
  const Vec3 & bearing_vector,
  const Vec3 & pt3D
)
{
  const auto new_bearing = (P * pt3D.homogeneous()).normalized();
  return 1.0 - (bearing_vector.dot(new_bearing));
}

} // namespace euclidean_resection
} // namespace openMVG
