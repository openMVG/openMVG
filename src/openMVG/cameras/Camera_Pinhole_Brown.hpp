// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Sida Li, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_HPP

#include <vector>

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

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
  using class_type = Pinhole_Intrinsic_Brown_T2;

  protected:

    /// center of distortion is applied by the Intrinsics class
    std::vector<double> params_; // K1, K2, K3, T1, T2

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
      params_ = {k1, k2, k3, t1, t2};
    }

    /**
    * @brief Get type of the intrinsic
    * @retval PINHOLE_CAMERA_BROWN
    */
    EINTRINSIC getType() const override
    {
      return PINHOLE_CAMERA_BROWN;
    }

    /**
    * @brief Does the camera model handle a distortion field?
    * @retval true
    */
    bool have_disto() const override
    {
      return true;
    }

    /**
    * @brief Add the distortion field to a point (that is in normalized camera frame)
    * @param p Point before distortion computation (in normalized camera frame)
    * @return point with distortion
    */
    Vec2 add_disto( const Vec2 & p ) const override
    {
      return ( p + distoFunction( params_, p ) );
    }

    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    * @note numerical approximation based on
    * Heikkila J (2000) Geometric Camera Calibration Using Circular Control Points.
    * IEEE Trans. Pattern Anal. Mach. Intell., 22:1066-1077
    */
    Vec2 remove_disto( const Vec2 & p ) const override
    {
      const double epsilon = 1e-10; //criteria to stop the iteration
      Vec2 p_u = p;

      Vec2 d = distoFunction(params_, p_u);
      while ((p_u + d - p).lpNorm<1>() > epsilon) //manhattan distance between the two points
      {
        p_u = p - d;
        d = distoFunction(params_, p_u);
      }

      return p_u;
    }

    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    std::vector<double> getParams() const override
    {
      std::vector<double> params = Pinhole_Intrinsic::getParams();
      params.insert(params.end(), std::begin(params_), std::end(params_));
      return params;
    }

    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    bool updateFromParams( const std::vector<double> & params ) override
    {
      if ( params.size() == 8 )
      {
        *this = Pinhole_Intrinsic_Brown_T2(
                  w_, h_,
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
    * @brief Return the list of parameter indexes that must be held constant
    * @param parametrization The given parametrization
    */
    std::vector<int> subsetParameterization
    (
      const Intrinsic_Parameter_Type & parametrization) const override
    {
      std::vector<int> constant_index;
      const int param = static_cast<int>(parametrization);
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH)
          || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(0);
      }
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT)
          || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(1);
        constant_index.push_back(2);
      }
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_DISTORTION)
          || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(3);
        constant_index.push_back(4);
        constant_index.push_back(5);
        constant_index.push_back(6);
        constant_index.push_back(7);
      }
      return constant_index;
    }

    /**
    * @brief Return the un-distorted pixel (with removed distortion)
    * @param p Input distorted pixel
    * @return Point without distortion
    */
    Vec2 get_ud_pixel( const Vec2& p ) const override
    {
      return cam2ima( remove_disto( ima2cam( p ) ) );
    }

    /**
    * @brief Return the distorted pixel (with added distortion)
    * @param p Input pixel
    * @return Distorted pixel
    */
    Vec2 get_d_pixel( const Vec2& p ) const override
    {
      return cam2ima( add_disto( ima2cam( p ) ) );
    }

    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    inline void save( Archive & ar ) const;

    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    inline void load( Archive & ar );

    /**
    * @brief Clone the object
    * @return A clone (copy of the stored object)
    */
    IntrinsicBase * clone( void ) const override
    {
      return new class_type( *this );
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

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_HPP
