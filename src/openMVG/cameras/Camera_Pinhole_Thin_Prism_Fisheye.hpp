// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong,Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This is an adaptation of the Fisheye distortion model implemented in OpenCV
// https://github.com/Itseez/opencv/blob/master/modules/calib3d/src/fisheye.cpp

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_THIN_PRISM_FISHEYE_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_THIN_PRISM_FISHEYE_HPP

#include <vector>

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"

namespace openMVG
{
namespace cameras
{

/**
* @brief Implement a thin_prism_fisheye camera model
* reference: https://www.eth3d.net/documentation#camera-models
*/
class Pinhole_Intrinsic_Thin_Prism_Fisheye : public Pinhole_Intrinsic
{
  using class_type = Pinhole_Intrinsic_Thin_Prism_Fisheye;

  protected:

    /// center of distortion is applied by the Intrinsics class
    std::vector<double> params_; // k1, k2, p1, p2, k3, k4, sx1, sx2


  public:

    /**
    * @brief Constructor
    * @param w Width of image plane
    * @param h Height of image plane
    * @param focal Focal distance in pixel
    * @param ppx Principal point on X-axis
    * @param ppy Principal point on Y-axis
    * @param k1 Distortion coefficient
    * @param k2 Distortion coefficient
    * @param p1 Distortion coefficient
    * @param p2 Distortion coefficient
    * @param k3 Distortion coefficient
    * @param k4 Distortion coefficient
    * @param sx1 Distortion coefficient
    * @param sx2 Distortion coefficient
    */
    Pinhole_Intrinsic_Thin_Prism_Fisheye(
      int w = 0, int h = 0,
      double focal = 0.0, double ppx = 0, double ppy = 0,
      double k1 = 0.0, double k2 = 0.0, 
      double p1 = 0.0, double p2 = 0.0,
      double k3 = 0.0, double k4 = 0.0,
      double sx1 = 0.0, double sx2 = 0.0 )
      : Pinhole_Intrinsic( w, h, focal, ppx, ppy ),
        params_({k1, k2, p1, p2, k3, k4, sx1, sx2})
    {
    }

    ~Pinhole_Intrinsic_Thin_Prism_Fisheye() override = default;

    /**
    * @brief Tell from which type the embed camera is
    * @retval PINHOLE_CAMERA_FISHEYE
    */
    EINTRINSIC getType() const override
    {
      return PINHOLE_CAMERA_THIN_PRISM_FISHEYE;
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
      if ( params.size() == 7 )
      {
        *this = Pinhole_Intrinsic_Fisheye(
                  w_, h_,
                  params[0], params[1], params[2], // focal, ppx, ppy
                  params[3], params[4], // k1, k2
                  params[5], params[6], // p1, p2
                  params[7], params[8], // k3, k4
                  params[9], params[10], // sx1, sx2
                   ); , 
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
        constant_index.push_back(8);
        constant_index.push_back(9);
        constant_index.push_back(10);
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
      const double eps = 1e-8;
      const double k1 = params[0], k2 = params[1];
      const double p1 = params[2], p1 = params[3];
      const double k3 = params[4], k4 = params[5];
      const double sx1 = params[6], sx2 = params[7];
      const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
      const double r  = sqrt(r2);
      if(r<eps)
      {
        return {p( 0 ), p( 1 )};
      }
      const double theta = atan(sqrt(r));
      const double theta_divide_r = theta/r;
      const double u_d = theta_divide_r * p( 0 );
      const double v_d = theta_divide_r * p( 1 );
      const double r4 = r2 * r2;
      const double r6 = r4 * r2;
      const double r8 = r4 * r4;
      const double t_r = 1 + k1*r2 + k2*r4 + k3*r6 + k4*r8;
      const double u_n = u_d * t_r + 2 * p1 * u_d * v_d + p_2 * (r2 + 2 * u_d * u_d) + sx1 * r2;
      const double v_n = v_d * t_r + 2 * p2 * u_d * v_d + p_1 * (r2 + 2 * v_d * v_d) + sx2 * r2;

      return { u_n, v_n};
    }
};


} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_THIN_PRISM_FISHEYE_HPP
