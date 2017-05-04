// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

namespace openMVG{
namespace geometry{

bool FindRTS
(
  const Mat &x1,
  const Mat &x2,
  double * S,
  Vec3 * t,
  Mat3 * R
)
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

lm_SRTRefine_functor::lm_SRTRefine_functor( int inputs, int values,
                      const Mat &x1, const Mat &x2,
                      const double &S, const Mat3 & R, const Vec &t ): Functor<double>( inputs, values ),
  x1_( x1 ), x2_( x2 ), t_( t ), R_( R ), S_( S ) { }

int lm_SRTRefine_functor::operator()( const Vec &x, Vec &fvec ) const
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


lm_RRefine_functor::lm_RRefine_functor( int inputs, int values,
                    const Mat &x1, const Mat &x2,
                    const double &S, const Mat3 & R, const Vec &t ): Functor<double>( inputs, values ),
  x1_( x1 ), x2_( x2 ), t_( t ), R_( R ), S_( S ) { }

int lm_RRefine_functor::operator()( const Vec &x, Vec &fvec ) const
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

void Refine_RTS
(
  const Mat &x1,
  const Mat &x2,
  double * S,
  Vec3 * t,
  Mat3 * R
)
{
  {
    lm_SRTRefine_functor functor( 7, 3 * x1.cols(), x1, x2, *S, *R, *t );

    Eigen::NumericalDiff<lm_SRTRefine_functor> numDiff( functor );

    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lm_SRTRefine_functor> > lm( numDiff );
    lm.parameters.maxfev = 1000;
    Vec xlm = Vec::Zero( 7 ); // The deviation vector {tx,ty,tz,rotX,rotY,rotZ,S}

    lm.minimize( xlm );

    const Vec3 transAdd = xlm.block<3, 1>( 0, 0 );
    const Vec3 rot = xlm.block<3, 1>( 3, 0 );
    const double SAdd = xlm( 6 );

    //Build the rotation matrix
    const Mat3 Rcor =
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

    const Vec3 rot = xlm.block<3, 1>( 0, 0 );

    //Build the rotation matrix
    const Mat3 Rcor =
      ( Eigen::AngleAxis<double>( rot( 0 ), Vec3::UnitX() )
        * Eigen::AngleAxis<double>( rot( 1 ), Vec3::UnitY() )
        * Eigen::AngleAxis<double>( rot( 2 ), Vec3::UnitZ() ) ).toRotationMatrix();

    *R = ( *R ) * Rcor;
  }
}

} // namespace geometry
} // namespace openMVG
