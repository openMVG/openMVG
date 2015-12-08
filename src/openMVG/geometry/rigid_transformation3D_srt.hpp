// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_3D_REGISTRATION_7DOF_H_
#define OPENMVG_GEOMETRY_3D_REGISTRATION_7DOF_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/lm.hpp"

#include <openMVG/robust_estimation/robust_estimator_ACRansac.hpp>

namespace openMVG {
namespace geometry {

/**
 * @brief Compose a similarity matrix given a scale factor, a rotation matrix and 
 * a translation vector.
 * 
 * @param[in] S The scale factor.
 * @param[in] t The translation vector.
 * @param[in] R The rotation matrix.
 * @param[out] RTS The 4x4 similarity matrix [S*R | t; 0 0 0 1].
 */
static void composeRTS(double &S, const Vec3 &t, const Mat3 &R, Mat4 &RTS)
{
  RTS.topLeftCorner(3, 3) = S*R;
  RTS.topRightCorner(3, 1) = t;
}

/**
 * @brief Decompose a similarity matrix into its scale, rotation and translation parts
 * 
 * @param[in] RTS The similarity matrix to decompose.
 * @param[out] S The scale factor.
 * @param[out] t The translation part.
 * @param[out] R The rotation part.
 * @return true if the input matrix is a similarity matrix.
 */
static bool decomposeRTS(const Mat4 &RTS, double &S, Vec3 &t, Mat3 &R)
{
  // Check critical cases
  R = RTS.topLeftCorner<3, 3>();
  if (R.determinant() < 0)
    return false;
  S = pow(R.determinant(), 1.0 / 3.0);
  // Check for degenerate case (if all points have the same value...)
  if (S < std::numeric_limits<double>::epsilon())
    return false;

  // Extract transformation parameters
  R /= S;
  t = RTS.topRightCorner<3, 1>();
  return true;
}

/** 3D rigid transformation estimation (7 dof)
 * Compute a Scale Rotation and Translation rigid transformation.
 * This transformation provide a distortion-free transformation
 * using the following formulation Xb = S * R * Xa + t.
 * "Least-squares estimation of transformation parameters between two point patterns",
 * Shinji Umeyama, PAMI 1991, DOI: 10.1109/34.88573
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
  double &S,
  Vec3 &t,
  Mat3 &R)
{
  if (x1.cols() < 3 || x2.cols() < 3)
    return false;

  assert(3 == x1.rows());
  assert(3 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  // Get the transformation via Umeyama's least squares algorithm. This returns
  // a matrix of the form:
  // [ s * R t]
  // [ 0 1]
  // from which we can extract the scale, rotation, and translation.
  const Eigen::Matrix4d transform = Eigen::umeyama(x1, x2, true);

  return decomposeRTS(transform, S, t, R);
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
  double &S,
  Vec3 &t,
  Mat3 &R)
{
  {
    lm_SRTRefine_functor functor(7, 3*x1.cols(), x1, x2, S, R, t);

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

    R = R * Rcor;
    t = t + transAdd;
    S = S + SAdd;
  }

  // Refine rotation
  {
    lm_RRefine_functor functor(3, 3*x1.cols(), x1, x2, S, R, t);

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

    R = R * Rcor;
  }
}

/**
 * @brief the Solver to use for ACRansac
 */
struct RTSSolver
{

  enum
  {
    MINIMUM_SAMPLES = 3
  };

  enum
  {
    MAX_MODELS = 1
  };
  
  // Solve the RTS problem
  static void Solve(const Mat &pts1, const Mat &pts2, std::vector<Mat4> *models)
  {
    models->push_back(Eigen::umeyama(pts1, pts2, true));
  }

  // Compute the residual of the transformation 
  static double Error(const Mat4 &RTS, const Vec3 & pt1, const Vec3 & pt2)
  {
    const Mat3 &RS = RTS.topLeftCorner<3, 3>();
    const Vec3 &t = RTS.topRightCorner<3, 1>();
    return (pt2 - (RS*pt1 + t)).norm();
  }
};

/**
 * @brief A functor that computes the squared error between points transformed by
 * a similarity
 */
struct RTSSquaredResidualError 
{
  // Return the squared error
  static double Error(const Mat4 &RTS, const Vec3 & pt1, const Vec3 & pt2)
  {
    const Mat3 &RS = RTS.topLeftCorner<3, 3>();
    const Vec3 &t = RTS.topRightCorner<3, 1>();
    return (pt2 - (RS*pt1 + t)).squaredNorm();
  }
};

/**
 * @brief The kernel to use for ACRansac
 */
template <typename SolverArg,
typename ErrorArg,
typename ModelArg = Mat4>
class ACKernelAdaptor_PointsRegistrationSRT
{
public:
  typedef SolverArg Solver;
  typedef ModelArg Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptor_PointsRegistrationSRT(const Mat & xA,
                                        const Mat & xB) :
             x1_(xA), x2_(xB), logalpha0_(log10(M_PI)) //@todo  WTF?
  {
    assert(3 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());
    
    // @todo normalize points?
  }

  enum
  {
    MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES
  };

  enum
  {
    MAX_MODELS = Solver::MAX_MODELS
  };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const
  {
    const Mat x1 = ExtractColumns(x1_, samples);
    const Mat x2 = ExtractColumns(x2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(size_t sample, const Model &model) const
  {
    return Square(ErrorT::Error(model, x1_.col(sample), x2_.col(sample)));
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const
  {
    vec_errors.resize(x1_.cols());
    for(size_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = Square(ErrorT::Error(model, x1_.col(sample), x2_.col(sample)));
  }

  size_t NumSamples() const { return static_cast<size_t> (x1_.cols()); }

  void Unnormalize(Model * model) const { } //-- Do nothing, no normalization 

  double logalpha0() const { return logalpha0_; }

  double multError() const { return 1.;}

  Mat3 normalizer1() const { return Mat3::Identity(); }

  Mat3 normalizer2() const { return Mat3::Identity(); }

  double unormalizeError(double val) const { return sqrt(val); }

private:
  Mat x1_, x2_;       // Normalized input data
  double logalpha0_;  // Alpha0 is used to make the error scale invariant
};

/**
 * @brief Uses AC ransac to robustly estimate the similarity between two sets of points.
 * 
 * @param[in] x1 The first 3xN matrix of euclidean points.
 * @param[in] x2 The second 3xN matrix of euclidean points.
 * @param[out] S The scale factor.
 * @param[out] t The 3x1 translation.
 * @param[out] R The 3x3 rotation.
 * @param[out] vec_inliers The vector containing the indices of inliers points.
 * @param[in] refine Enable/Disable refining of the found transformation.
 * @return true if the found transformation is a similarity
 * @see FindRTS()
 */
static bool ACRansac_FindRTS(const Mat &x1,
                             const Mat &x2,
                             double &S, 
                             Vec3 &t, 
                             Mat3 &R, 
                             std::vector<std::size_t> &vec_inliers,
                             bool refine = false)
{
  assert(3 == x1.rows());
  assert(3 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());
  
  const std::size_t numIterations = 1024;
  const double dPrecision = std::numeric_limits<double>::infinity();
  
  Mat4 RTS;

  typedef geometry::RTSSolver SolverType;
  typedef ACKernelAdaptor_PointsRegistrationSRT<
          SolverType,
          geometry::RTSSquaredResidualError> KernelType;

  KernelType kernel = KernelType(x1, x2);
  // Robust estimation of the Projection matrix and its precision
  const std::pair<double, double> ACRansacOut =
          robust::ACRANSAC(kernel, vec_inliers, numIterations, &RTS, dPrecision, true);
  
  const bool good = decomposeRTS(RTS, S, t, R);
  
  // return if it is not good or refine is not required
  if(!good || !refine)
    return good;

  const std::size_t nbInliers = vec_inliers.size();
  //only refine the inliers
  Mat inliers1 = Mat3X(3, nbInliers);
  Mat inliers2 = Mat3X(3, nbInliers);

  for(std::size_t i = 0; i < nbInliers; ++i)
  {
    inliers1.col(i) = x1.col(vec_inliers[i]);
    inliers2.col(i) = x2.col(vec_inliers[i]);
  }

  geometry::Refine_RTS(inliers1, inliers2, S, t, R);
  
  return good;

}

/**
 * @brief Uses AC ransac to robustly estimate the similarity between two sets of 
 * points. Just a wrapper that output the similarity in matrix form
 * @param[in] x1 The first 3xN matrix of euclidean points.
 * @param[in] x2 The second 3xN matrix of euclidean points.
 * @param[out] RTS The 4x4 similarity matrix. 
 * @param[out] vec_inliers The inliers used to estimate the similarity.
 * @param  refine Enable/Disable refining of the found transformation.
 * @return true if the found transformation is a similarity
 * @see geometry::FindRTS()
 * @see geometry::ACRansac_FindRTS()
 */
static bool ACRansac_FindRTS(const Mat &x1, 
                             const Mat &x2, 
                             Mat4 &RTS, 
                             std::vector<std::size_t> &vec_inliers,
                             bool refine = false)
{
  double S;
  Vec3 t; 
  Mat3 R; 
  const bool good = ACRansac_FindRTS(x1, x2, S, t, R, vec_inliers, refine);
  
  if(good)
    composeRTS(S, t, R, RTS);
  
  return good;
}

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_REGISTRATION_7DOF_H_
