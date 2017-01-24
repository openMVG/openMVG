
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

#include "openMVG/multiview/solver_resection_kernel.hpp"

#include <cassert>

namespace openMVG {
namespace resection {
namespace kernel {

using namespace std;

void translate(
  const Mat3X & X,
  const Vec3 & vecTranslation,
  Mat3X * XPoints)
{
  XPoints->resize(X.rows(), X.cols());
  for (int i=0; i<X.cols(); ++i)
    XPoints->col(i) = X.col(i) + vecTranslation;
}

template <typename TMat, typename TVec>
double NullspaceRatio(TMat *A, TVec *nullspace) {
  if (A->rows() >= A->cols()) {
    Eigen::JacobiSVD<TMat> svd(*A, Eigen::ComputeFullV);
    (*nullspace) = svd.matrixV().col(A->cols()-1);
    return svd.singularValues()(A->cols()-2) / svd.singularValues()(0);
  }
  // Extend A with rows of zeros to make it square. It's a hack, but is
  // necessary until Eigen supports SVD with more columns than rows.
  TMat A_extended(A->cols(), A->cols());
  A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
  A_extended.block(0,0, A->rows(), A->cols()) = (*A);
  return Nullspace(&A_extended, nullspace);
}

/// Setup the Direct Linear Transform.
///  Use template in order to support fixed or dynamic sized matrix.
template<typename Matrix>
void BuildActionMatrix(Matrix & A, const Mat &pt2D, const Mat &XPoints)  {

  const size_t n = pt2D.cols();
  for (size_t i = 0; i < n; ++i) {
    size_t row_index = i * 2;
    const Vec3 & X = XPoints.col(i);
    const Vec2 & x = pt2D.col(i);
    A(row_index,  0) =  X(0);
    A(row_index,  1) =  X(1);
    A(row_index,  2) =  X(2);
    A(row_index,  3) =  1.0;
    A(row_index,  8) = -X(0) * x(0);
    A(row_index,  9) = -X(1) * x(0);
    A(row_index, 10) = -X(2) * x(0);
    A(row_index, 11) = -1.0 * x(0);

    row_index = i * 2 + 1;
    A(row_index,  4) =  X(0);
    A(row_index,  5) =  X(1);
    A(row_index,  6) =  X(2);
    A(row_index,  7) =  1.0;
    A(row_index,  8) = -X(0) * x(1);
    A(row_index,  9) = -X(1) * x(1);
    A(row_index, 10) = -X(2) * x(1);
    A(row_index, 11) = -1.0 * x(1);
  }
  // Normalize each row
  for (size_t i = 0; i < static_cast<size_t>(A.rows()); ++i)
    A.row(i).normalize();
}

void SixPointResectionSolver::Solve(
  const Mat &pt2D,
  const Mat &pt3d,
  vector<Mat34> *Ps,
  bool bcheck)
{
  assert(2 == pt2D.rows());
  assert(3 == pt3d.rows());
  assert(6 <= pt2D.cols());
  assert(pt2D.cols() == pt3d.cols());

  //-- Translate 3D points in order to have X0 = (0,0,0,1).
  Vec3 vecTranslation = - pt3d.col(0);
  Mat4 translationMatrix = Mat4::Identity();
  translationMatrix << 1, 0, 0, vecTranslation(0),
                       0, 1, 0, vecTranslation(1),
                       0, 0, 1, vecTranslation(2),
                       0, 0, 0, 1;
  Mat3X XPoints;
  translate(pt3d, vecTranslation, &XPoints);

  const size_t n = pt2D.cols();

  using Vec12 = Eigen::Matrix<double, 12, 1>;
  Vec12 p;
  double ratio = -1.0;
  if (n==6) {
    // In the case of minimal configuration we use fixed sized matrix to let
    //  Eigen and the compiler doing the maximum of optimization.
    using Mat12 = Eigen::Matrix<double, 12, 12>;
    Mat12 A = Mat12::Zero(12, 12);
    BuildActionMatrix(A, pt2D, XPoints);
    ratio = NullspaceRatio(&A, &p);
  }
  else  {
    Mat A = Mat::Zero(n*2, 12);
    BuildActionMatrix(A, pt2D, XPoints);
    ratio = NullspaceRatio(&A, &p);
  }
  if (bcheck) {
    if (ratio > 1e-5) //Assert that at least only one solution if found by SVD
    {
      Mat34 P = Map<Mat>(p.data(),4,3).transpose();
      P = P * translationMatrix;
      P /= P(2,3);

      Mat3 K, R;
      Vec3 t;
      KRt_From_P(P,&K,&R,&t);

      //Assert point in front of the cam
      size_t cpt = 0;
      for (size_t i = 0; i < n; ++i) {
        cpt += (Depth(R, t, pt3d.col(i))>0) ? 1 : 0;
      }
      if (cpt == n) {
        Ps->push_back(P);
      }
    }
  }
  else  {
    Mat34 P = Map<Mat>(p.data(),4,3).transpose();
    P = P * translationMatrix;
    P /= P(2,3);
    Ps->push_back(P);
  }
}

}  // namespace kernel
}  // namespace resection
}  // namespace openMVG

namespace openMVG {
namespace euclidean_resection {
namespace kernel {

// Selects 4 virtual control points using mean and PCA.
void SelectControlPoints(const Mat3X &X_world,
                         Mat *X_centered,
                         Mat34 *X_control_points) {
  size_t num_points = X_world.cols();

  // The first virtual control point, C0, is the centroid.
  Vec mean, variance;
  MeanAndVarianceAlongRows(X_world, &mean, &variance);
  X_control_points->col(0) = mean;

  // Computes PCA
  X_centered->resize (3, num_points);
  for (size_t c = 0; c < num_points; c++) {
    X_centered->col(c) = X_world.col (c) - mean;
  }
  Mat3 X_centered_sq = (*X_centered) * X_centered->transpose();
  Eigen::JacobiSVD<Mat3> X_centered_sq_svd(X_centered_sq, Eigen::ComputeFullU);
  Vec3 w = X_centered_sq_svd.singularValues();
  Mat3 u = X_centered_sq_svd.matrixU();
  for (size_t c = 0; c < 3; c++) {
    double k = sqrt (w (c) / num_points);
    X_control_points->col (c + 1) = mean + k * u.col (c);
  }
}

// Computes the barycentric coordinates for all real points
void ComputeBarycentricCoordinates(const Mat3X &X_world_centered,
                                   const Mat34 &X_control_points,
                                   Mat4X *alphas) {
  size_t num_points = X_world_centered.cols();
  Mat3 C2 ;
  for (size_t c = 1; c < 4; c++) {
    C2.col(c-1) = X_control_points.col(c) - X_control_points.col(0);
  }

  Mat3 C2inv = C2.inverse();
  Mat3X a = C2inv * X_world_centered;

  alphas->resize(4, num_points);
  alphas->setZero();
  alphas->block(1, 0, 3, num_points) = a;
  for (size_t c = 0; c < num_points; c++) {
    (*alphas)(0, c) = 1.0 - alphas->col(c).sum();
  }
}

// Estimates the coordinates of all real points in the camera coordinate frame
void ComputePointsCoordinatesInCameraFrame(
    const Mat4X &alphas,
    const Vec4 &betas,
    const Eigen::Matrix<double, 12, 12> &U,
    Mat3X *X_camera) {
  size_t num_points = alphas.cols();

  // Estimates the control points in the camera reference frame.
  Mat34 C2b; C2b.setZero();
  for (size_t cu = 0; cu < 4; cu++) {
    for (size_t c = 0; c < 4; c++) {
      C2b.col(c) += betas(cu) * U.block(11 - cu, c * 3, 1, 3).transpose();
    }
  }

  // Estimates the 3D points in the camera reference frame
  X_camera->resize(3, num_points);
  for (size_t c = 0; c < num_points; c++) {
    X_camera->col(c) = C2b * alphas.col(c);
  }

  // Check the sign of the z coordinate of the points (should be positive)
  size_t num_z_neg = 0;
  for (Mat3X::Index i = 0; i < X_camera->cols(); ++i) {
    if ((*X_camera)(2,i) < 0) {
      num_z_neg++;
    }
  }

  // If more than 50% of z are negative, we change the signs
  if (num_z_neg > 0.5 * X_camera->cols()) {
    C2b = -C2b;
    *X_camera = -(*X_camera);
  }
}

void AbsoluteOrientation(const Mat3X &X,
                         const Mat3X &Xp,
                         Mat3 *R,
                         Vec3 *t) {
  int num_points = static_cast<int>(X.cols());
  Vec3 C  = X.rowwise().sum() / num_points;   // Centroid of X.
  Vec3 Cp = Xp.rowwise().sum() / num_points;  // Centroid of Xp.

  // Normalize the two point sets.
  Mat3X Xn(3, num_points), Xpn(3, num_points);
  for( int i = 0; i < num_points; ++i ){
    Xn.col(i)  = X.col(i) - C;
    Xpn.col(i) = Xp.col(i) - Cp;
  }

  // Construct the N matrix (pg. 635).
  double Sxx = Xn.row(0).dot(Xpn.row(0));
  double Syy = Xn.row(1).dot(Xpn.row(1));
  double Szz = Xn.row(2).dot(Xpn.row(2));
  double Sxy = Xn.row(0).dot(Xpn.row(1));
  double Syx = Xn.row(1).dot(Xpn.row(0));
  double Sxz = Xn.row(0).dot(Xpn.row(2));
  double Szx = Xn.row(2).dot(Xpn.row(0));
  double Syz = Xn.row(1).dot(Xpn.row(2));
  double Szy = Xn.row(2).dot(Xpn.row(1));

  Mat4 N;
  N << Sxx + Syy + Szz, Syz - Szy,        Szx - Sxz,        Sxy - Syx,
       Syz - Szy,       Sxx - Syy - Szz,  Sxy + Syx,        Szx + Sxz,
       Szx - Sxz,       Sxy + Syx,       -Sxx + Syy - Szz,  Syz + Szy,
       Sxy - Syx,       Szx + Sxz,        Syz + Szy,       -Sxx - Syy + Szz;

  // Find the unit quaternion q that maximizes qNq. It is the eigenvector
  // corresponding to the largest eigenvalue.
  Vec4 q = N.jacobiSvd(Eigen::ComputeFullU).matrixU().col(0);

  // Retrieve the 3x3 rotation matrix.
  Vec4 qq = q.array() * q.array();
  double q0q1 = q(0) * q(1);
  double q0q2 = q(0) * q(2);
  double q0q3 = q(0) * q(3);
  double q1q2 = q(1) * q(2);
  double q1q3 = q(1) * q(3);
  double q2q3 = q(2) * q(3);

  (*R) << qq(0) + qq(1) - qq(2) - qq(3),
          2 * (q1q2 - q0q3),
          2 * (q1q3 + q0q2),
          2 * (q1q2+ q0q3),
          qq(0) - qq(1) + qq(2) - qq(3),
          2 * (q2q3 - q0q1),
          2 * (q1q3 - q0q2),
          2 * (q2q3 + q0q1),
          qq(0) - qq(1) - qq(2) + qq(3);

  // Fix the handedness of the R matrix.
  if (R->determinant() < 0) {
    R->row(2) = -R->row(2);
  }
  // Compute the final translation.
  *t = Cp - *R * C;
}


bool EuclideanResectionEPnP(const Mat2X &x_camera,
                            const Mat3X &X_world,
                            Mat3 *R, Vec3 *t) {
  assert(x_camera.cols() == X_world.cols());
  assert(x_camera.cols() > 3);
  size_t num_points = X_world.cols();

  // Select the control points.
  Mat34 X_control_points;
  Mat X_centered;
  SelectControlPoints(X_world, &X_centered, &X_control_points);

  // Compute the barycentric coordinates.
  Mat4X alphas(4, num_points);
  ComputeBarycentricCoordinates(X_centered, X_control_points, &alphas);

  // Estimates the M matrix with the barycentric coordinates
  Mat M(2 * num_points, 12);
  Eigen::Matrix<double, 2, 12> sub_M;
  for (size_t c = 0; c < num_points; c++) {
    double a0 = alphas(0, c);
    double a1 = alphas(1, c);
    double a2 = alphas(2, c);
    double a3 = alphas(3, c);
    double ui = x_camera(0, c);
    double vi = x_camera(1, c);
    M.block(2*c, 0, 2, 12) << a0, 0,
                              a0*(-ui), a1, 0,
                              a1*(-ui), a2, 0,
                              a2*(-ui), a3, 0,
                              a3*(-ui), 0,
                              a0, a0*(-vi), 0,
                              a1, a1*(-vi), 0,
                              a2, a2*(-vi), 0,
                              a3, a3*(-vi);
  }

  // TODO(julien): Avoid the transpose by rewriting the u2.block() calls.
  Eigen::JacobiSVD<Mat> MtMsvd(M.transpose()*M, Eigen::ComputeFullU);
  Eigen::Matrix<double, 12, 12> u2 = MtMsvd.matrixU().transpose();

  // Estimate the L matrix.
  Eigen::Matrix<double, 6, 3> dv1;
  Eigen::Matrix<double, 6, 3> dv2;
  Eigen::Matrix<double, 6, 3> dv3;
  Eigen::Matrix<double, 6, 3> dv4;

  dv1.row(0) = u2.block(11, 0, 1, 3) - u2.block(11, 3, 1, 3);
  dv1.row(1) = u2.block(11, 0, 1, 3) - u2.block(11, 6, 1, 3);
  dv1.row(2) = u2.block(11, 0, 1, 3) - u2.block(11, 9, 1, 3);
  dv1.row(3) = u2.block(11, 3, 1, 3) - u2.block(11, 6, 1, 3);
  dv1.row(4) = u2.block(11, 3, 1, 3) - u2.block(11, 9, 1, 3);
  dv1.row(5) = u2.block(11, 6, 1, 3) - u2.block(11, 9, 1, 3);
  dv2.row(0) = u2.block(10, 0, 1, 3) - u2.block(10, 3, 1, 3);
  dv2.row(1) = u2.block(10, 0, 1, 3) - u2.block(10, 6, 1, 3);
  dv2.row(2) = u2.block(10, 0, 1, 3) - u2.block(10, 9, 1, 3);
  dv2.row(3) = u2.block(10, 3, 1, 3) - u2.block(10, 6, 1, 3);
  dv2.row(4) = u2.block(10, 3, 1, 3) - u2.block(10, 9, 1, 3);
  dv2.row(5) = u2.block(10, 6, 1, 3) - u2.block(10, 9, 1, 3);
  dv3.row(0) = u2.block( 9, 0, 1, 3) - u2.block( 9, 3, 1, 3);
  dv3.row(1) = u2.block( 9, 0, 1, 3) - u2.block( 9, 6, 1, 3);
  dv3.row(2) = u2.block( 9, 0, 1, 3) - u2.block( 9, 9, 1, 3);
  dv3.row(3) = u2.block( 9, 3, 1, 3) - u2.block( 9, 6, 1, 3);
  dv3.row(4) = u2.block( 9, 3, 1, 3) - u2.block( 9, 9, 1, 3);
  dv3.row(5) = u2.block( 9, 6, 1, 3) - u2.block( 9, 9, 1, 3);
  dv4.row(0) = u2.block( 8, 0, 1, 3) - u2.block( 8, 3, 1, 3);
  dv4.row(1) = u2.block( 8, 0, 1, 3) - u2.block( 8, 6, 1, 3);
  dv4.row(2) = u2.block( 8, 0, 1, 3) - u2.block( 8, 9, 1, 3);
  dv4.row(3) = u2.block( 8, 3, 1, 3) - u2.block( 8, 6, 1, 3);
  dv4.row(4) = u2.block( 8, 3, 1, 3) - u2.block( 8, 9, 1, 3);
  dv4.row(5) = u2.block( 8, 6, 1, 3) - u2.block( 8, 9, 1, 3);

  Eigen::Matrix<double, 6, 10> L;
  for (size_t r = 0; r < 6; r++) {
    L.row(r) << dv1.row(r).dot(dv1.row(r)),
          2.0 * dv1.row(r).dot(dv2.row(r)),
                dv2.row(r).dot(dv2.row(r)),
          2.0 * dv1.row(r).dot(dv3.row(r)),
          2.0 * dv2.row(r).dot(dv3.row(r)),
                dv3.row(r).dot(dv3.row(r)),
          2.0 * dv1.row(r).dot(dv4.row(r)),
          2.0 * dv2.row(r).dot(dv4.row(r)),
          2.0 * dv3.row(r).dot(dv4.row(r)),
                dv4.row(r).dot(dv4.row(r));
  }
  Vec rho;
  rho.resize(6);
  rho << (X_control_points.col(0) - X_control_points.col(1)).squaredNorm(),
         (X_control_points.col(0) - X_control_points.col(2)).squaredNorm(),
         (X_control_points.col(0) - X_control_points.col(3)).squaredNorm(),
         (X_control_points.col(1) - X_control_points.col(2)).squaredNorm(),
         (X_control_points.col(1) - X_control_points.col(3)).squaredNorm(),
         (X_control_points.col(2) - X_control_points.col(3)).squaredNorm();

  // There are three possible solutions based on the three approximations of L
  // (betas). Below, each one is solved for then the best one is chosen.
  Mat3X X_camera;
  Mat3 K; K.setIdentity();
  std::vector<Mat3> Rs(3);
  std::vector<Vec3> ts(3);
  Vec rmse(3);

  bool bSol = false;

  // Find the first possible solution for R, t corresponding to:
  // Betas          = [b00 b01 b11 b02 b12 b22 b03 b13 b23 b33]
  // Betas_approx_1 = [b00 b01     b02         b03]
  Vec4 betas = Vec4::Zero();
  Eigen::Matrix<double, 6, 4> l_6x4;
  for (size_t r = 0; r < 6; r++) {
    l_6x4.row(r) << L(r, 0), L(r, 1), L(r, 3), L(r, 6);
  }
  Eigen::JacobiSVD<Mat> svd_of_l4(l_6x4,
                                  Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vec4 b4 = svd_of_l4.solve(rho);
  if ((l_6x4 * b4).isApprox(rho, 1e-3)) {
    if (b4(0) < 0) {
      b4 = -b4;
    }
    b4(0) =  std::sqrt(b4(0));
    betas << b4(0), b4(1) / b4(0), b4(2) / b4(0), b4(3) / b4(0);
    ComputePointsCoordinatesInCameraFrame(alphas, betas, u2, &X_camera);
    AbsoluteOrientation(X_world, X_camera, &Rs[0], &ts[0]);
    rmse(0) = RootMeanSquareError(x_camera, X_world, K, Rs[0], ts[0]);
    bSol = true;
  } else {
    //std::cerr << "First approximation of beta not good enough.";
    ts[0].setZero();
    rmse(0) = std::numeric_limits<double>::max();
  }

  // Find the second possible solution for R, t corresponding to:
  // Betas          = [b00 b01 b11 b02 b12 b22 b03 b13 b23 b33]
  // Betas_approx_2 = [b00 b01 b11]
  betas.setZero();
  Eigen::Matrix<double, 6, 3> l_6x3;
  l_6x3 = L.block(0, 0, 6, 3);
  Eigen::JacobiSVD<Mat> svdOfL3(l_6x3,
                                Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vec3 b3 = svdOfL3.solve(rho);
  if ((l_6x3 * b3).isApprox(rho, 1e-3)) {
    if (b3(0) < 0) {
      betas(0) = std::sqrt(-b3(0));
      betas(1) = (b3(2) < 0) ? std::sqrt(-b3(2)) : 0;
    } else {
      betas(0) = std::sqrt(b3(0));
      betas(1) = (b3(2) > 0) ? std::sqrt(b3(2)) : 0;
    }
    if (b3(1) < 0) {
      betas(0) = -betas(0);
    }
    betas(2) = 0;
    betas(3) = 0;
    ComputePointsCoordinatesInCameraFrame(alphas, betas, u2, &X_camera);
    AbsoluteOrientation(X_world, X_camera, &Rs[1], &ts[1]);
    rmse(1) = RootMeanSquareError(x_camera, X_world, K, Rs[1], ts[1]);
    bSol = true;
  } else {
    //std::cerr << "Second approximation of beta not good enough.";
    ts[1].setZero();
    rmse(1) = std::numeric_limits<double>::max();
  }

  // Find the third possible solution for R, t corresponding to:
  // Betas          = [b00 b01 b11 b02 b12 b22 b03 b13 b23 b33]
  // Betas_approx_3 = [b00 b01 b11 b02 b12]
  betas.setZero();
  Eigen::Matrix<double, 6, 5> l_6x5;
  l_6x5 = L.block(0, 0, 6, 5);
  Eigen::JacobiSVD<Mat> svdOfL5(l_6x5,
                                Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vec b5 = svdOfL5.solve(rho);
  if ((l_6x5 * b5).isApprox(rho, 1e-3)) {
    if (b5(0) < 0) {
      betas(0) = std::sqrt(-b5(0));
      if (b5(2) < 0) {
        betas(1) = std::sqrt(-b5(2));
      } else {
        b5(2) = 0;
      }
    } else {
      betas(0) = std::sqrt(b5(0));
      if (b5(2) > 0) {
        betas(1) = std::sqrt(b5(2));
      } else {
        b5(2) = 0;
      }
    }
    if (b5(1) < 0) {
      betas(0) = -betas(0);
    }
    betas(2) = b5(3) / betas(0);
    betas(3) = 0;
    ComputePointsCoordinatesInCameraFrame(alphas, betas, u2, &X_camera);
    AbsoluteOrientation(X_world, X_camera, &Rs[2], &ts[2]);
    rmse(2) = RootMeanSquareError(x_camera, X_world, K, Rs[2], ts[2]);
    bSol = true;
  } else {
    //std::cerr << "Third approximation of beta not good enough.";
    ts[2].setZero();
    rmse(2) = std::numeric_limits<double>::max();
  }

  // Finally, with all three solutions, select the (R, t) with the best RMSE.
  size_t n = 0;
  if (rmse(1) < rmse(0)) {
    n = 1;
  }
  if (rmse(2) < rmse(n)) {
    n = 2;
  }
  if (bSol)  { //If at least one solution have been found
    *R = Rs[n];
    *t = ts[n];
    return true;
  }

  return false;
}


}  // namespace kernel
}  // namespace euclidean_resection
}  // namespace openMVG

