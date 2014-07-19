// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_3D_REGISTRATION_7DOF_H_
#define OPENMVG_GEOMETRY_3D_REGISTRATION_7DOF_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/lm.hpp"

namespace openMVG {
namespace geometry {

/** 3D rigid transformation estimation (7 dof)
 * Compute a Scale Rotation and Translation rigid transformation.
 * This transformation provide a distortion-free transformation
 * using the following formulation Xb = S * R * Xa + t.
 * Haralick, Robert, Shapiro, Linda. Computer and Robot Vision book, 1992.
 *
 * \param[in] x1 The first 3xN matrix of euclidean points
 * \param[in] x2 The second 3xN matrix of euclidean points
 * \param[out] S The scale factor
 * \param[out] t The 3x1 translation
 * \param[out] R The 3x3 rotation
 *
 * \return true if the transformation estimation has succeeded
 *
 * \note Need at least 3 points
 */

static bool FindRTS(const Mat &x1,
  const Mat &x2,
  double * S,
  Vec3 * t,
  Mat3 * R)
{
  assert(3 == x1.rows());
  assert(3 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  const int n = static_cast<int>(x1.cols());

  // Compute scale factor
  double S1 =0.0, S2=0.0;
  S1 = (x1.array()*x1.array()).sum();
  S2 = (x2.array()*x2.array()).sum();

  *S = sqrt(S2/S1);

  // Compute centroid and variance
  Vec meanx1, meanx2;
  Vec varx1, varx2;
  MeanAndVarianceAlongRows(x1, &meanx1, &varx1);
  MeanAndVarianceAlongRows(x2, &meanx2, &varx2);

  // Compute the n*n correlation matrix
  Mat p1c(3,n) , p2c(3,n);

  for(int i = 0; i < n; ++i){
    p1c.col(i) = x1.col(i) - meanx1;
    p2c.col(i) = x2.col(i) - meanx2;
  }

  Mat H = p1c * p2c.transpose();

  Eigen::JacobiSVD<Mat> svdObj(H,Eigen::ComputeFullV | Eigen::ComputeFullU);
  Mat V = svdObj.matrixV();
  Mat U = svdObj.matrixU();
  *R = V * U.transpose();

  if( R->determinant() < 0){
    return false;
  }

  //Compute the optimal translation
  *t = meanx2 - (*S) * (*R) * meanx1;

  return true;
}

// Eigen LM functor to refine translation, Rotation and Scale parameter.
struct lm_SRTRefine_functor : Functor<double>
{
  lm_SRTRefine_functor(int inputs, int values,
    const Mat &x1, const Mat &x2,
    const double &S, const Mat3 & R, const Vec &t): Functor<double>(inputs,values),
    _x1(x1), _x2(x2), _t(t), _R(R), _S(S) { }

  int operator()(const Vec &x, Vec &fvec) const
  {
    // convert x to rotation matrix and a translation vector and a Scale factor
    // x = {tx,ty,tz,anglex,angley,anglez,S}
    Vec3 transAdd = x.block<3,1>(0,0);
    Vec3 rot = x.block<3,1>(3,0);
    double Sadd = x(6);

    //Build the rotation matrix
    Mat3 Rcor =
      (Eigen::AngleAxis<double>(rot(0), Vec3::UnitX())
      * Eigen::AngleAxis<double>(rot(1), Vec3::UnitY())
      * Eigen::AngleAxis<double>(rot(2), Vec3::UnitZ())).toRotationMatrix();

    const Mat3 nR  = _R*Rcor;
    const Vec3 nt = _t+transAdd;
    const double nS = _S+Sadd;

    // Evaluate re-projection errors
    Vec3 proj;
    for (Mat::Index i = 0; i < _x1.cols(); ++i)
    {
      proj = _x2.col(i) -  (nS *  nR * (_x1.col(i)) + nt);
      fvec[i*3]   = proj(0);
      fvec[i*3+1] = proj(1);
      fvec[i*3+2] = proj(2);
    }
    return 0;
  }

  Mat _x1, _x2;
  Vec3 _t;
  Mat3 _R;
  double _S;
};

// Eigen LM functor to refine Rotation.
struct lm_RRefine_functor : Functor<double>
{
  lm_RRefine_functor(int inputs, int values,
    const Mat &x1, const Mat &x2,
    const double &S, const Mat3 & R, const Vec &t): Functor<double>(inputs,values),
    _x1(x1), _x2(x2), _t(t), _R(R), _S(S) { }

  int operator()(const Vec &x, Vec &fvec) const
  {
    // convert x to rotation matrix
    // x = {anglex,angley,anglez}
    Vec3 rot = x.block<3,1>(0,0);

    //Build the rotation matrix
    Mat3 Rcor =
      (Eigen::AngleAxis<double>(rot(0), Vec3::UnitX())
      * Eigen::AngleAxis<double>(rot(1), Vec3::UnitY())
      * Eigen::AngleAxis<double>(rot(2), Vec3::UnitZ())).toRotationMatrix();

    const Mat3 nR  = _R*Rcor;
    const Vec3 nt = _t;
    const double nS = _S;

    // Evaluate re-projection errors
    Vec3 proj;
    for (Mat::Index i = 0; i < _x1.cols(); ++i)
    {
      proj = _x2.col(i) -  (nS *  nR * (_x1.col(i)) + nt);
      fvec[i*3]   = proj(0);
      fvec[i*3+1] = proj(1);
      fvec[i*3+2] = proj(2);
    }
    return 0;
  }

  Mat _x1, _x2;
  Vec3 _t;
  Mat3 _R;
  double _S;
};

/** 3D rigid transformation refinement using LM
 * Refine the Scale, Rotation and translation rigid transformation
 * using a Levenberg-Marquardt opimization.
 *
 * \param[in] x1 The first 3xN matrix of euclidean points
 * \param[in] x2 The second 3xN matrix of euclidean points
 * \param[out] S The initial scale factor
 * \param[out] t The initial 3x1 translation
 * \param[out] R The initial 3x3 rotation
 *
 * \return none
 */
static void Refine_RTS(const Mat &x1,
  const Mat &x2,
  double * S,
  Vec3 * t,
  Mat3 * R)
{
  {
    lm_SRTRefine_functor functor(7, 3*x1.cols(), x1, x2, *S, *R, *t);

    Eigen::NumericalDiff<lm_SRTRefine_functor> numDiff(functor);

    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lm_SRTRefine_functor> > lm(numDiff);
    lm.parameters.maxfev = 1000;
    Vec xlm = Vec::Zero(7); // The deviation vector {tx,ty,tz,rotX,rotY,rotZ,S}

    lm.minimize(xlm);

    Vec3 transAdd = xlm.block<3,1>(0,0);
    Vec3 rot = xlm.block<3,1>(3,0);
    double SAdd = xlm(6);

    //Build the rotation matrix
    Mat3 Rcor =
      (Eigen::AngleAxis<double>(rot(0), Vec3::UnitX())
      * Eigen::AngleAxis<double>(rot(1), Vec3::UnitY())
      * Eigen::AngleAxis<double>(rot(2), Vec3::UnitZ())).toRotationMatrix();

    *R = (*R)*Rcor;
    *t = (*t)+transAdd;
    *S = (*S)+SAdd;
  }

  // Refine rotation
  {
    lm_RRefine_functor functor(3, 3*x1.cols(), x1, x2, *S, *R, *t);

    Eigen::NumericalDiff<lm_RRefine_functor> numDiff(functor);

    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lm_RRefine_functor> > lm(numDiff);
    lm.parameters.maxfev = 1000;
    Vec xlm = Vec::Zero(3); // The deviation vector {rotX,rotY,rotZ}

    lm.minimize(xlm);

    Vec3 rot = xlm.block<3,1>(0,0);

    //Build the rotation matrix
    Mat3 Rcor =
      (Eigen::AngleAxis<double>(rot(0), Vec3::UnitX())
      * Eigen::AngleAxis<double>(rot(1), Vec3::UnitY())
      * Eigen::AngleAxis<double>(rot(2), Vec3::UnitZ())).toRotationMatrix();

    *R = (*R)*Rcor;
  }
}

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_REGISTRATION_7DOF_H_
