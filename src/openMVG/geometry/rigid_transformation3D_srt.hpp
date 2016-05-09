// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_3D_REGISTRATION_7DOF_H_
#define OPENMVG_GEOMETRY_3D_REGISTRATION_7DOF_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/lm.hpp"

namespace openMVG
{
namespace geometry
{

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
static bool FindRTS( const Mat &x1,
                     const Mat &x2,
                     double * S,
                     Vec3 * t,
                     Mat3 * R )
{
  if ( x1.cols() < 3 || x2.cols() < 3 )
  {
    return false;
  }

  assert( 3 == x1.rows() );
  assert( 3 <= x1.cols() );
  assert( x1.rows() == x2.rows() );
  assert( x1.cols() == x2.cols() );

  // Get the transformation via Umeyama's least squares algorithm. This returns
  // a matrix of the form:
  // [ s * R t]
  // [ 0 1]
  // from which we can extract the scale, rotation, and translation.
  const Eigen::Matrix4d transform = Eigen::umeyama( x1, x2, true );

  // Check critical cases
  *R = transform.topLeftCorner<3, 3>();
  if ( R->determinant() < 0 )
  {
    return false;
  }
  *S = pow( R->determinant(), 1.0 / 3.0 );
  // Check for degenerate case (if all points have the same value...)
  if ( *S < std::numeric_limits<double>::epsilon() )
  {
    return false;
  }

  // Extract transformation parameters
  *S = pow( R->determinant(), 1.0 / 3.0 );
  *R /= *S;
  *t = transform.topRightCorner<3, 1>();

  return true;
}

/**
* @brief Eigen Levemberg-Marquardt functor to refine translation, Rotation and Scale parameter.
*/
struct lm_SRTRefine_functor : Functor<double>
{
  /**
  * @brief Constructor
  * @param inputs Number of inputs (nb elements to refine)
  * @param values Number of samples
  * @param x1 Input samples for first dataset
  * @param x2 Input samples for second dataset
  * @param S Scale
  * @param R Rotation
  * @param t Translation
  */
  lm_SRTRefine_functor( int inputs, int values,
                        const Mat &x1, const Mat &x2,
                        const double &S, const Mat3 & R, const Vec &t ): Functor<double>( inputs, values ),
    x1_( x1 ), x2_( x2 ), t_( t ), R_( R ), S_( S ) { }

  /**
  * @brief Computes error given a sample
  * @param x a Sample
  * @param[out] fvec Error for each values
  */
  int operator()( const Vec &x, Vec &fvec ) const
  {
    // convert x to rotation matrix and a translation vector and a Scale factor
    // x = {tx,ty,tz,anglex,angley,anglez,S}
    Vec3 transAdd = x.block<3, 1>( 0, 0 );
    Vec3 rot = x.block<3, 1>( 3, 0 );
    double Sadd = x( 6 );

    //Build the rotation matrix
    Mat3 Rcor =
      ( Eigen::AngleAxis<double>( rot( 0 ), Vec3::UnitX() )
        * Eigen::AngleAxis<double>( rot( 1 ), Vec3::UnitY() )
        * Eigen::AngleAxis<double>( rot( 2 ), Vec3::UnitZ() ) ).toRotationMatrix();

    const Mat3 nR  = R_ * Rcor;
    const Vec3 nt = t_ + transAdd;
    const double nS = S_ + Sadd;

    // Evaluate re-projection errors
    Vec3 proj;
    for ( Mat::Index i = 0; i < x1_.cols(); ++i )
    {
      proj = x2_.col( i ) -  ( nS *  nR * ( x1_.col( i ) ) + nt );
      fvec[i * 3]   = proj( 0 );
      fvec[i * 3 + 1] = proj( 1 );
      fvec[i * 3 + 2] = proj( 2 );
    }
    return 0;
  }

  Mat x1_, x2_;
  Vec3 t_;
  Mat3 R_;
  double S_;
};


/**
* @brief Eigen LM functor to refine Rotation.
*/
struct lm_RRefine_functor : Functor<double>
{
  /**
  * @brief Constructor
  * @param inputs Number of inputs (elements to refine)
  * @param values Number of samples
  * @param x1 Input samples for first dataset
  * @param x2 Input samples for second dataset
  * @param S Scale
  * @param R Rotation
  * @param t Translation
  */
  lm_RRefine_functor( int inputs, int values,
                      const Mat &x1, const Mat &x2,
                      const double &S, const Mat3 & R, const Vec &t ): Functor<double>( inputs, values ),
    x1_( x1 ), x2_( x2 ), t_( t ), R_( R ), S_( S ) { }

  /**
   * @brief Computes error given a sample
   * @param x a Sample
   * @param[out] fvec Error for each values
   */
  int operator()( const Vec &x, Vec &fvec ) const
  {
    // convert x to rotation matrix
    // x = {anglex,angley,anglez}
    Vec3 rot = x.block<3, 1>( 0, 0 );

    //Build the rotation matrix
    Mat3 Rcor =
      ( Eigen::AngleAxis<double>( rot( 0 ), Vec3::UnitX() )
        * Eigen::AngleAxis<double>( rot( 1 ), Vec3::UnitY() )
        * Eigen::AngleAxis<double>( rot( 2 ), Vec3::UnitZ() ) ).toRotationMatrix();

    const Mat3 nR  = R_ * Rcor;
    const Vec3 nt = t_;
    const double nS = S_;

    // Evaluate re-projection errors
    Vec3 proj;
    for ( Mat::Index i = 0; i < x1_.cols(); ++i )
    {
      proj = x2_.col( i ) -  ( nS *  nR * ( x1_.col( i ) ) + nt );
      fvec[i * 3]   = proj( 0 );
      fvec[i * 3 + 1] = proj( 1 );
      fvec[i * 3 + 2] = proj( 2 );
    }
    return 0;
  }

  Mat x1_, x2_;
  Vec3 t_;
  Mat3 R_;
  double S_;
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
static void Refine_RTS( const Mat &x1,
                        const Mat &x2,
                        double * S,
                        Vec3 * t,
                        Mat3 * R )
{
  {
    lm_SRTRefine_functor functor( 7, 3 * x1.cols(), x1, x2, *S, *R, *t );

    Eigen::NumericalDiff<lm_SRTRefine_functor> numDiff( functor );

    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lm_SRTRefine_functor> > lm( numDiff );
    lm.parameters.maxfev = 1000;
    Vec xlm = Vec::Zero( 7 ); // The deviation vector {tx,ty,tz,rotX,rotY,rotZ,S}

    lm.minimize( xlm );

    Vec3 transAdd = xlm.block<3, 1>( 0, 0 );
    Vec3 rot = xlm.block<3, 1>( 3, 0 );
    double SAdd = xlm( 6 );

    //Build the rotation matrix
    Mat3 Rcor =
      ( Eigen::AngleAxis<double>( rot( 0 ), Vec3::UnitX() )
        * Eigen::AngleAxis<double>( rot( 1 ), Vec3::UnitY() )
        * Eigen::AngleAxis<double>( rot( 2 ), Vec3::UnitZ() ) ).toRotationMatrix();

    *R = ( *R ) * Rcor;
    *t = ( *t ) + transAdd;
    *S = ( *S ) + SAdd;
  }

  // Refine rotation
  {
    lm_RRefine_functor functor( 3, 3 * x1.cols(), x1, x2, *S, *R, *t );

    Eigen::NumericalDiff<lm_RRefine_functor> numDiff( functor );

    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lm_RRefine_functor> > lm( numDiff );
    lm.parameters.maxfev = 1000;
    Vec xlm = Vec::Zero( 3 ); // The deviation vector {rotX,rotY,rotZ}

    lm.minimize( xlm );

    Vec3 rot = xlm.block<3, 1>( 0, 0 );

    //Build the rotation matrix
    Mat3 Rcor =
      ( Eigen::AngleAxis<double>( rot( 0 ), Vec3::UnitX() )
        * Eigen::AngleAxis<double>( rot( 1 ), Vec3::UnitY() )
        * Eigen::AngleAxis<double>( rot( 2 ), Vec3::UnitZ() ) ).toRotationMatrix();

    *R = ( *R ) * Rcor;
  }
}

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_REGISTRATION_7DOF_H_
