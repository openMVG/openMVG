
// Copyright (c) 2015 Sida Li.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_PINHOLE_BROWN_HPP
#define OPENMVG_CAMERA_PINHOLE_BROWN_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Common.hpp"

#include <vector>

namespace openMVG
{
namespace cameras
{

/**
* @brief Implement a Pinhole camera with a 3 radial distortion coefficients and 2 tangential distortion coefficients.
* \f$ x_d = x_u (1 + K_1 r^2 + K_2 r^4 + K_3 r^6) + (T_2 (r^2 + 2 x_u^2) + 2 T_1 x_u y_u) \f$
* \f$ y_d = y_u (1 + K_1 r^2 + K_2 r^4 + K_3 r^6) + (T_1 (r^2 + 2 y_u^2) + 2 T_2 x_u y_u) \f$
*/
class Pinhole_Intrinsic_Brown_T2 : public Pinhole_Intrinsic
{
  protected:

    /// center of distortion is applied by the Intrinsics class
    std::vector<double> _params; // K1, K2, K3, T1, T2

  public:

    /**
    * @brief Constructor
    * @param w Width of image
    * @param h Height of image
    * @param focal Focal distance (in pixel)
    * @param ppx Principal point on X-axis
    * @param ppy Principal point on Y-axis
    * @param k1 First radial distortion coefficient
    * @param k2 Second radial distortion coefficient
    * @param k3 Third radial distortion coefficient
    * @param t1 First tangential distortion coefficient
    * @param t2 Second tangential distortion coefficient
    */
    Pinhole_Intrinsic_Brown_T2(
      int w = 0, int h = 0,
      double focal = 0.0, double ppx = 0, double ppy = 0,
      double k1 = 0.0, double k2 = 0.0, double k3 = 0.0,
      double t1 = 0.0, double t2 = 0.0 )
      : Pinhole_Intrinsic( w, h, focal, ppx, ppy )
    {
      _params = {k1, k2, k3, t1, t2};
    }

    /**
    * @brief Get type of the intrinsic
    * @retval PINHOLE_CAMERA_BROWN
    */
    EINTRINSIC getType() const
    {
      return PINHOLE_CAMERA_BROWN;
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
      return ( p + distoFunction( _params, p ) );
    }

    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    * @note numerical approximation based on
    * Heikkila J (2000) Geometric Camera Calibration Using Circular Control Points.
    * IEEE Trans. Pattern Anal. Mach. Intell., 22:1066-1077
    */
    virtual Vec2 remove_disto( const Vec2 & p ) const
    {
      const double epsilon = 1e-8; //criteria to stop the iteration
      Vec2 p_u = p;

      while( ( add_disto( p_u ) - p ).lpNorm<1>() > epsilon ) //manhattan distance between the two points
      {
        p_u = p - distoFunction( _params, p_u );
      }

      return p_u;
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
      params.push_back( _params[3] );
      params.push_back( _params[4] );
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
      if ( params.size() == 8 )
      {
        *this = Pinhole_Intrinsic_Brown_T2(
                  _w, _h,
                  params[0], params[1], params[2], // focal, ppx, ppy
                  params[3], params[4], params[5], // K1, K2, K3
                  params[6], params[7] );          // T1, T2
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
      ar( cereal::make_nvp( "disto_t2", _params ) );
    }


    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar )
    {
      Pinhole_Intrinsic::load( ar );
      ar( cereal::make_nvp( "disto_t2", _params ) );
    }

  private:


    /**
    * @brief Functor to calculate distortion offset accounting for both radial and tangential distortion
    * @param params List of parameters to define a Brown camera
    * @param p Input point
    * @return Transformed point
    */
    static Vec2 distoFunction( const std::vector<double> & params, const Vec2 & p )
    {
      const double k1 = params[0], k2 = params[1], k3 = params[2], t1 = params[3], t2 = params[4];
      const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
      const double r4 = r2 * r2;
      const double r6 = r4 * r2;
      const double k_diff = ( k1 * r2 + k2 * r4 + k3 * r6 );
      const double t_x = t2 * ( r2 + 2 * p( 0 ) * p( 0 ) ) + 2 * t1 * p( 0 ) * p( 1 );
      const double t_y = t1 * ( r2 + 2 * p( 1 ) * p( 1 ) ) + 2 * t2 * p( 0 ) * p( 1 );
      Vec2 d( p( 0 ) * k_diff + t_x, p( 1 ) * k_diff + t_y );
      return d;
    }
};


} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Brown_T2, "pinhole_brown_t2" );

#endif // #ifndef OPENMVG_CAMERA_PINHOLE_BROWN_HPP
