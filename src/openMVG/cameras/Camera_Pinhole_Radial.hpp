
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_PINHOLE_RADIAL_K_HPP
#define OPENMVG_CAMERA_PINHOLE_RADIAL_K_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Common.hpp"

#include <vector>

namespace openMVG
{
namespace cameras
{

/**
* @brief Namespace handling radial distortion computation
*/
namespace radial_distortion
{


/**
* @brief Solve by bisection the p' radius such that Square(disto(radius(p'))) = r^2
* @param params Parameters of the distortion
* @param r2 Target radius
* @param functor Functor used to paramatrize the distortion
* @param epsilon Error driven threshold
* @return Best radius
*/
template <class Disto_Functor>
double bisection_Radius_Solve(
  const std::vector<double> & params, // radial distortion parameters
  double r2, // targeted radius
  Disto_Functor & functor,
  double epsilon = 1e-8 // criteria to stop the bisection
)
{
  // Guess plausible upper and lower bound
  double lowerbound = r2, upbound = r2;
  while ( functor( params, lowerbound ) > r2 )
  {
    lowerbound /= 1.05;
  }
  while ( functor( params, upbound ) < r2 )
  {
    upbound *= 1.05;
  }

  // Perform a bisection until epsilon accuracy is not reached
  while ( epsilon < upbound - lowerbound )
  {
    const double mid = .5 * ( lowerbound + upbound );
    if ( functor( params, mid ) > r2 )
    {
      upbound = mid;
    }
    else
    {
      lowerbound = mid;
    }
  }
  return .5 * ( lowerbound + upbound );
}

} // namespace radial_distortion

/**
 * @brief Implement a Pinhole camera with a 1 radial distortion coefficient.
 * \f$ x_d = x_u (1 + K_1 r^2 ) \f$
 */
class Pinhole_Intrinsic_Radial_K1 : public Pinhole_Intrinsic
{
  protected:
    /// center of distortion is applied by the Intrinsics class
    std::vector<double> _params; // K1

  public:

    /**
    * @brief Constructor
    * @param w Width of the image
    * @param h Height of the image
    * @param focal Focal (in pixel) of the camera
    * @param ppx Principal point on X-Axis
    * @param ppy Principal point on Y-Axis
    * @param k1 Distortion coefficient
    */
    Pinhole_Intrinsic_Radial_K1(
      int w = 0, int h = 0,
      double focal = 0.0, double ppx = 0, double ppy = 0,
      double k1 = 0.0 )
      : Pinhole_Intrinsic( w, h, focal, ppx, ppy )
    {
      _params.resize( 1 );
      _params[0] = k1;
    }

    /**
    * @brief Tell from which type the embed camera is
    * @retval PINHOLE_CAMERA_RADIAL1
    */
    EINTRINSIC getType() const
    {
      return PINHOLE_CAMERA_RADIAL1;
    }

    /**
    * @brief Does the camera model handle a distortion field?
    * @retval true if intrinsic holds distortion
    * @retval false if intrinsic does not hold distortion
    */
    virtual bool have_disto() const
    {
      return true;
    }


    /**
    * @brief Add the distortion field to a point (that is in normalized camera frame)
    * @param p Point before distortion computation (in normalized camera frame)
    * @return point with distortion
    */
    virtual Vec2 add_disto( const Vec2 & p ) const
    {

      const double k1 = _params[0];

      const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
      const double r_coeff = ( 1. + k1 * r2 );

      return ( p * r_coeff );
    }


    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    */
    virtual Vec2 remove_disto( const Vec2& p ) const
    {
      // Compute the radius from which the point p comes from thanks to a bisection
      // Minimize disto(radius(p')^2) == actual Squared(radius(p))

      const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
      const double radius = ( r2 == 0 ) ?
                            1. :
                            ::sqrt( radial_distortion::bisection_Radius_Solve( _params, r2, distoFunctor ) / r2 );
      return radius * p;
    }


    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    virtual std::vector<double> getParams() const
    {
      std::vector<double> params = Pinhole_Intrinsic::getParams();
      params.push_back( _params[0] );
      return params;
    }


    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    virtual bool updateFromParams( const std::vector<double> & params )
    {
      if ( params.size() == 4 )
      {
        *this = Pinhole_Intrinsic_Radial_K1(
                  _w, _h,
                  params[0], params[1], params[2], // focal, ppx, ppy
                  params[3] ); //K1
        return true;
      }
      else
      {
        return false;
      }
    }


    /**
    * @brief Return the un-distorted pixel (with removed distortion)
    * @param p Input distorted pixel
    * @return Point without distortion
    */
    virtual Vec2 get_ud_pixel( const Vec2& p ) const
    {
      return cam2ima( remove_disto( ima2cam( p ) ) );
    }


    /**
    * @brief Return the distorted pixel (with added distortion)
    * @param p Input pixel
    * @return Distorted pixel
    */
    virtual Vec2 get_d_pixel( const Vec2& p ) const
    {
      return cam2ima( add_disto( ima2cam( p ) ) );
    }


    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    void save( Archive & ar ) const
    {
      Pinhole_Intrinsic::save( ar );
      ar( cereal::make_nvp( "disto_k1", _params ) );
    }


    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar )
    {
      Pinhole_Intrinsic::load( ar );
      ar( cereal::make_nvp( "disto_k1", _params ) );
    }

  private:


    /**
    * @brief Functor to solve Square(disto(radius(p'))) = r^2
    * @param params List of parameters (only the first one is used)
    * @param r2 square distance (relative to center)
    * @return distance
    */
    static double distoFunctor( const std::vector<double> & params, double r2 )
    {
      const double k1 = params[0];
      return r2 * Square( 1. + r2 * k1 );
    }
};

/**
* @brief Implement a Pinhole camera with a 3 radial distortion coefficients.
* \f$ x_d = x_u (1 + K_1 r^2 + K_2 r^4 + K_3 r^6) \f$
*/
class Pinhole_Intrinsic_Radial_K3 : public Pinhole_Intrinsic
{
  protected:
    // center of distortion is applied by the Intrinsics class
    /// K1, K2, K3
    std::vector<double> _params;

  public:

    /**
    * @brief Constructor
    * @param w Width of image
    * @param h Height of image
    * @param focal Focal (in pixel) of the camera
    * @param ppx Principal point on X-Axis
    * @param ppy Principal point on Y-Axis
    * @param k1 First radial distortion coefficient
    * @param k2 Second radial distortion coefficient
    * @param k3 Third radial distortion coefficient
    */
    Pinhole_Intrinsic_Radial_K3(
      int w = 0, int h = 0,
      double focal = 0.0, double ppx = 0, double ppy = 0,
      double k1 = 0.0, double k2 = 0.0, double k3 = 0.0 )
      : Pinhole_Intrinsic( w, h, focal, ppx, ppy )
    {
      _params.resize( 3 );
      _params[0] = k1;
      _params[1] = k2;
      _params[2] = k3;
    }

    /**
    * @brief Tell from which type the embed camera is
    * @retval PINHOLE_CAMERA_RADIAL3
    */
    EINTRINSIC getType() const
    {
      return PINHOLE_CAMERA_RADIAL3;
    }

    /**
    * @brief Does the camera model handle a distortion field?
    * @retval true
    */
    virtual bool have_disto() const
    {
      return true;
    }


    /**
    * @brief Add the distortion field to a point (that is in normalized camera frame)
    * @param p Point before distortion computation (in normalized camera frame)
    * @return point with distortion
    */
    virtual Vec2 add_disto( const Vec2 & p ) const
    {

      const double k1 = _params[0], k2 = _params[1], k3 = _params[2];

      const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
      const double r4 = r2 * r2;
      const double r6 = r4 * r2;
      const double r_coeff = ( 1. + k1 * r2 + k2 * r4 + k3 * r6 );

      return ( p * r_coeff );
    }


    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    */
    virtual Vec2 remove_disto( const Vec2& p ) const
    {
      // Compute the radius from which the point p comes from thanks to a bisection
      // Minimize disto(radius(p')^2) == actual Squared(radius(p))

      const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
      const double radius = ( r2 == 0 ) ? //1. : ::sqrt(bisectionSolve(_params, r2) / r2);
                            1. :
                            ::sqrt( radial_distortion::bisection_Radius_Solve( _params, r2, distoFunctor ) / r2 );
      return radius * p;
    }


    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    virtual std::vector<double> getParams() const
    {
      std::vector<double> params = Pinhole_Intrinsic::getParams();
      params.push_back( _params[0] );
      params.push_back( _params[1] );
      params.push_back( _params[2] );
      return params;
    }


    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    virtual bool updateFromParams( const std::vector<double> & params )
    {
      if ( params.size() == 6 )
      {
        *this = Pinhole_Intrinsic_Radial_K3(
                  _w, _h,
                  params[0], params[1], params[2], // focal, ppx, ppy
                  params[3], params[4], params[5] ); // K1, K2, K3
        return true;
      }
      else
      {
        return false;
      }
    }


    /**
    * @brief Return the un-distorted pixel (with removed distortion)
    * @param p Input distorted pixel
    * @return Point without distortion
    */
    virtual Vec2 get_ud_pixel( const Vec2& p ) const
    {
      return cam2ima( remove_disto( ima2cam( p ) ) );
    }


    /**
    * @brief Return the distorted pixel (with added distortion)
    * @param p Input pixel
    * @return Distorted pixel
    */
    virtual Vec2 get_d_pixel( const Vec2& p ) const
    {
      return cam2ima( add_disto( ima2cam( p ) ) );
    }


    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    void save( Archive & ar ) const
    {
      Pinhole_Intrinsic::save( ar );
      ar( cereal::make_nvp( "disto_k3", _params ) );
    }


    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar )
    {
      Pinhole_Intrinsic::load( ar );
      ar( cereal::make_nvp( "disto_k3", _params ) );
    }

  private:


    /**
    * @brief Functor to solve Square(disto(radius(p'))) = r^2
    * @param params List of parameters (only the first one is used)
    * @param r2 square distance (relative to center)
    * @return distance
    */
    static double distoFunctor( const std::vector<double> & params, double r2 )
    {
      const double k1 = params[0], k2 = params[1], k3 = params[2];
      return r2 * Square( 1. + r2 * ( k1 + r2 * ( k2 + r2 * k3 ) ) );
    }
};

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Radial_K1, "pinhole_radial_k1" );
CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Radial_K3, "pinhole_radial_k3" );

#endif // #ifndef OPENMVG_CAMERA_PINHOLE_RADIAL_K_HPP

