// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.
// Copyright (c) 2015 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_SAMPLE_HPP
#define OPENMVG_IMAGE_SAMPLE_HPP

#include <limits>
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"

namespace openMVG
{
namespace image
{

/**
** Sampling functors
** These functors computes weight associated to each pixels

  For a (relative) sampling position x (\in [0,1]) between two (consecutives) points :

    A  .... x ......... B

  w[0] is the weight associated to A
  w[1] is the weight associated to B

  Note: The following functors generalize the sampling to more than two neighbors
  ** They all contains the neighbor_width variable that specify the number of neighbors used for sampling
*/


/**
 ** Nearest sampling (ie: find the nearest pixel to a specified position)
 **/
struct SamplerNearest
{
  public:
    // Nearest sampling is only between two pixels
    static const int neighbor_width = 2;

    /**
     ** @brief Computes weight associated to neighboring pixels
     ** @author Romuald Perrot <perrot.romuald_AT_gmail.com>
     ** @param x Sampling position
     ** @param[out] weigth Sampling factors associated to the neighboring
     ** @note weight must be at least neighbor_width length
     **/
    void operator()( const double x , double * const weight ) const
    {
      weight[0] = ( x < 0.5 ) ? 1.0 : 0.0;
      weight[1] = ( x >= 0.5 ) ? 1.0 : 0.0;
    }
};


/**
 ** Linear sampling (ie: linear interpolation between two pixels)
 **/
struct SamplerLinear
{
  public:
    // Linear sampling is between two pixels
    static const int neighbor_width = 2;

    /**
     ** @brief Computes weight associated to neighboring pixels
     ** @author Romuald Perrot <perrot.romuald_AT_gmail.com>
     ** @param x Sampling position
     ** @param[out] weigth Sampling factors associated to the neighboring
     ** @note weight must be at least neighbor_width length
     **/
    void operator()( const double x , double * const weigth ) const
    {
      weigth[0] = 1.0 - x;
      weigth[1] = x;
    }
};


/**
 ** Cubic interpolation between 4 pixels
 **
 ** Interpolation weight is for A,B,C and D pixels given a x position as illustrated as follow :
 **
 ** A      B    x C      D
 **
 ** @ref : Cubic Convolution Interpolation for Digital Image Processing , R. Keys, eq(4)
 **/
struct SamplerCubic
{
  public:
    // Cubic interpolation is between 4 pixels
    static const int neighbor_width = 4;

    /**
     ** @brief Constructor
     ** @param sharpness_coef Sharpness coefficient used to control sharpness of the cubic curve
     ** @note sharpnedd_coef must be between -0.75 to -0.5
     ** -0.5 gives better mathematically result (ie: approximation at 3 order precision)
     **/
    explicit SamplerCubic( const double sharpness_coef = -0.5 )
      : sharpness_( sharpness_coef )
    {

    }

    /**
     ** @brief Computes weight associated to neighboring pixels
     ** @author Romuald Perrot <perrot.romuald_AT_gmail.com>
     ** @param x Sampling position
     ** @param[out] weigth Sampling factors associated to the neighboring
     ** @note weight must be at least neighbor_width length
     **/
    void operator()( const double x , double * const weigth ) const
    {
      // remember :
      // A      B    x  C       D

      // weigth[0] -> weight for A
      // weight[1] -> weight for B
      // weight[2] -> weight for C
      // weight[3] -> weigth for D

      weigth[0] = CubicInter12( sharpness_ , x + 1.0 );
      weigth[1] = CubicInter01( sharpness_ , x );
      weigth[2] = CubicInter01( sharpness_ , 1.0 - x );
      weigth[3] = CubicInter12( sharpness_ , 2.0 - x );
    }

  private:

    /**
     * @brief Cubic interpolation for x \in [0; 1 [
     * @param sharpness Cubic parameter
     * @param x Position where to sample the cubic spline
     */
    static double CubicInter01( const double sharpness , const double x )
    {
      // A = sharpness
      //  f(x) = ( A + 2 ) * x ^ 3 - ( A + 3 ) * x ^ 2 + 1
      return ( ( sharpness + 2.0 ) * x - ( sharpness + 3.0 ) ) * x * x + 1.0;
    }

    /**
     * @brief Cubic interpolation for x \in [1; 2 [
     * @param sharpness Cubic parameter
     * @param x Position where to sample the cubic spline
     */
    static double CubicInter12( const double sharpness , const double x )
    {

      // A = sharpness
      // f(x) = A * x^3 - 5 * A * x^2 + 8 * A * x - 4 * a

      return ( ( sharpness * x - 5.0 * sharpness ) * x + 8.0 * sharpness ) * x - 4.0 * sharpness;
    }

    /// Sharpness coefficient
    double sharpness_;
};

/**
 ** Sampler spline16 -> Interpolation on 4 points used for 2d ressampling (16 = 4x4 sampling)
 ** Cubic interpolation with 0-derivative at edges (ie at A and D points)
 ** See Helmut Dersch for more details
 **
 ** Some refs :
 **  -   http://forum.doom9.org/archive/index.php/t-147117.html
 **  -   http://avisynth.nl/index.php/Resampling
 **  -   http://www.ipol.im/pub/art/2011/g_lmii/
 **
 ** The idea is to consider 3 cubic splines (f1,f2,f3) in the sampling interval :
 **
 ** A   f1    B    f2    C     f3     D
 **
 ** with curves defined as follow :
 ** f1(x) = a1 x^3 + b1 x^2 + c1 x + d1
 ** f2(x) = a2 x^3 + b2 x^2 + c2 x + d2
 ** f3(x) = a3 x^3 + b2 x^2 + c3 x + d3
 **
 ** We want to compute spline coefs for A,B,C,D assuming that:
 ** y0 = coef[A] = f1(-1)
 ** y1 = coef[B] = f1(0) = f2(0)
 ** y2 = coef[C] = f2(1) = f3(1)
 ** y3 = coef[D] = f3(2)
 **
 ** coef are computed using the following constraints :
 ** Curve is continuous, ie:
 ** f1(0)  = f2(0)
 ** f2(1)  = f3(1)
 ** First derivative are equals, ie:
 ** f1'(0) = f2'(0)
 ** f2'(1) = f3'(1)
 ** Second derivative are equals, ie:
 ** f1''(0) = f2''(0)
 ** f2''(1) = f3''(0)
 ** Curve is, at boundary, with second derivative set to zero (it's a constraint introduced by Dersch), ie:
 ** f1''(-1) = 0
 ** f3''(2) = 0
 **
 ** Then, you can solve for (a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3)
 **
 ** for ex, for curve f2 you find :
 **
 ** d2 = y1                                      // easy since y1 = f2(0)
 ** c2 = - 7/15 y0 - 1/5 y1 + 4/5 y2 - 2/15 y3
 ** b2 = 4/5 y0 - 9/5 y1 + 6/5 y2 - 1/5 y3
 ** a2 = - 1/3 y0 + y1 - y2 + 1/3 y3
 **
 **
 ** When you have coefs, you just have to express your curve as a linear combinaison of the control points, fort ex
 ** with f2 :
 **
 **
 ** f2(x) = w0(x) * y0 + w1(x) + y1 + w2(x) * y2 + w3(x) * y3
 **
 ** with :
 **
 ** w0(x) = - 1/3 * x^3 + 4/5 * x^2 - 7/15 * x
 ** w1(x) = x^3 - 9/5 * x^2 - 1/5 * x + 1
 ** w2(x) = -x^3 + 6/5 * x^2 + 4/5 * x
 ** w3(x) = 1/3 * x^3 - 1/5 * x^2 - 2/15 * x
 **
 ** substituing boundary conditions gives the correct coeficients for y0,y1,y2,y3 giving the final sampling scheme
 **/
struct SamplerSpline16
{
  public:
    static const int neighbor_width = 4;


    /**
     ** @brief Computes weight associated to neighboring pixels
     ** @author Romuald Perrot <perrot.romuald_AT_gmail.com>
     ** @param x Sampling position
     ** @param[out] weigth Sampling factors associated to the neighboring
     ** @note weight must be at least neighbor_width length
     **/
    void operator()( const double x , double * const weigth ) const
    {
      weigth[0] = ( ( -1.0 / 3.0 * x + 4.0 / 5.0     ) * x -   7.0 / 15.0 ) * x;
      weigth[1] = ( ( x - 9.0 / 5.0 ) * x -   1.0 / 5.0     ) * x + 1.0;
      weigth[2] = ( ( 6.0 / 5.0 - x     ) * x +   4.0 / 5.0 ) * x;
      weigth[3] = ( ( 1.0 / 3.0  * x - 1.0 / 5.0 ) * x -   2.0 / 15.0 ) * x;
    }
};

/**
 ** Sampler spline 36
 ** Same as spline 16 but on 6 neighbors (used for 6x6 frame)
 **/
struct SamplerSpline36
{
  public:
    static const int neighbor_width = 6;

    /**
     ** @brief Computes weight associated to neighboring pixels
     ** @author Romuald Perrot <perrot.romuald_AT_gmail.com>
     ** @param x Sampling position
     ** @param[out] weigth Sampling factors associated to the neighboring
     ** @note weight must be at least neighbor_width length
     **/
    void operator()( const double x , double * const weigth ) const
    {
      weigth[0] = ( (    1.0 / 11.0  * x -  45.0 / 209.0 ) * x +  26.0 / 209.0  ) * x;
      weigth[1] = ( ( -  6.0 / 11.0  * x + 270.0 / 209.0 ) * x - 156.0 / 209.0  ) * x;
      weigth[2] = ( (   13.0 / 11.0  * x - 453.0 / 209.0 ) * x -   3.0 / 209.0  ) * x + 1.0;
      weigth[3] = ( ( - 13.0 / 11.0  * x + 288.0 / 209.0 ) * x + 168.0 / 209.0  ) * x;
      weigth[4] = ( (    6.0 / 11.0  * x -  72.0 / 209.0 ) * x -  42.0 / 209.0  ) * x;
      weigth[5] = ( ( -  1.0 / 11.0  * x +  12.0 / 209.0 ) * x +   7.0 / 209.0  ) * x;
    }
};

/**
 ** Sampler spline 64
 ** Same as spline 16 but on 8 neighbors (used for 8x8 frame)
 **/
struct SamplerSpline64
{
  public:
    static const int neighbor_width = 8;

    /**
     ** @brief Computes weight associated to neighboring pixels
     ** @author Romuald Perrot <perrot.romuald_AT_gmail.com>
     ** @param x Sampling position
     ** @param[out] weigth Sampling factors associated to the neighboring
     ** @note weight must be at least neighbor_width length
     **/
    void operator()( const double x , double * const weigth ) const
    {
      weigth[0] = ( ( - 1.0 / 41.0 * x +  168.0 / 2911.0 ) * x -   97.0 / 2911.0 ) * x;
      weigth[1] = ( (  6.0 / 41.0 * x - 1008.0 / 2911.0 ) * x +  582.0 / 2911.0 ) * x;
      weigth[2] = ( ( -24.0 / 41.0 * x + 4032.0 / 2911.0 ) * x - 2328.0 / 2911.0 ) * x;
      weigth[3] = ( ( 49.0 / 41.0 * x - 6387.0 / 2911.0 ) * x -    3.0 / 2911.0 ) * x + 1.0;
      weigth[4] = ( ( -49.0 / 41.0 * x + 4050.0 / 2911.0 ) * x + 2340.0 / 2911.0 ) * x;
      weigth[5] = ( ( 24.0 / 41.0 * x - 1080.0 / 2911.0 ) * x -  624.0 / 2911.0 ) * x;
      weigth[6] = ( ( - 6.0 / 41.0 * x +  270.0 / 2911.0 ) * x +  156.0 / 2911.0 ) * x;
      weigth[7] = ( (  1.0 / 41.0 * x -   45.0 / 2911.0 ) * x -   26.0 / 2911.0 ) * x;
    }
};


/**
* @brief Template class used to convert a pixel from a type to pixel in double format
*/
template<typename T>
struct RealPixel
{
  using base_type = T;
  using real_type = double;

  /**
  * @brief Cast pixel value to real
  * @param val Input value
  * @return casted value
  */
  static real_type convert_to_real( const base_type & val )
  {
    return static_cast<real_type>( val );
  }

  /**
  * @brief Cast pixel value to base_type
  * @param val Input value
  * @return casted value
  */
  static base_type convert_from_real( const real_type & val )
  {
    return static_cast<base_type>( val );
  }
};

/**
* @brief overloading for unsigned char
*/
template<>
struct RealPixel<unsigned char>
{
  using base_type = unsigned char;
  using real_type = double;

  /**
  * @brief Cast pixel value to real
  * @param val Input value
  * @return casted value
  */
  static real_type convert_to_real( const base_type & val )
  {
    return static_cast<real_type>( val );
  }

  /**
  * @brief Cast pixel value to base_type
  * @param val Input value
  * @return casted value
  */
  static base_type convert_from_real( const real_type & val )
  {
    // handle out of range values.
    return ( val < 0.0 ) ?
           0 :
           ( val > static_cast<real_type>( std::numeric_limits<base_type>::max() ) ?
             std::numeric_limits<base_type>::max() :
             static_cast<base_type>( val + 0.5 ) );
  }
};


/**
* @brief overloading for float
*/
template<>
struct RealPixel<float>
{
  using base_type = float;
  using real_type = double;

  /**
  * @brief Convert pixel value to real_type
  * @param val Input value
  * @return casted value
  */
  static real_type convert_to_real( const base_type & val )
  {
    return val;
  }

  /**
  * @brief Convert pixel value to base_type
  * @param val Input type
  * @return casted value
  */
  static base_type convert_from_real( const real_type & val )
  {
    return val;
  }
};


/**
* @brief overloading for Rgb
*/
template<typename T>
struct RealPixel<Rgb<T>>
{
  using base_type = Rgb<T>;
  using real_type = Rgb<double>;

  /**
  * @brief Convert pixel value to real_type
  * @param val Input value
  * @return casted value
  */
  static real_type convert_to_real( const base_type & val )
  {
    return real_type( val.template cast<double>() );
  }

  /**
  * @brief Convert pixel value to base_type
  * @param val Input type
  * @return casted value
  */
  static base_type convert_from_real( const real_type & val )
  {
    return base_type( RealPixel<T>::convert_from_real( val.r() ) ,
                      RealPixel<T>::convert_from_real( val.g() ) ,
                      RealPixel<T>::convert_from_real( val.b() ) );
  }
};


/**
* @brief overloading for rgba
*/
template<typename T>
struct RealPixel<Rgba<T>>
{
  using base_type = Rgba<T>;
  using real_type = Rgba<double>;

  /**
  * @brief Convert pixel value to real_type
  * @param val Input value
  * @return casted value
  */
  static real_type convert_to_real( const base_type & val )
  {
    return real_type( val.template cast<double>() );
  }

  /**
  * @brief Convert pixel value to base_type
  * @param val Input type
  * @return casted value
  */
  static base_type convert_from_real( const real_type & val )
  {
    return base_type( RealPixel<T>::convert_from_real( val.r() ) ,
                      RealPixel<T>::convert_from_real( val.g() ) ,
                      RealPixel<T>::convert_from_real( val.b() ) ,
                      RealPixel<T>::convert_from_real( val.a() ) );
  }
};

/**
 ** Generic sampling of image using a sampling function
 **/
template<typename SamplerFunc>
struct Sampler2d
{
    explicit Sampler2d( const SamplerFunc & sampler = SamplerFunc() )
      : sampler_( sampler ) ,
        half_width_( SamplerFunc::neighbor_width / 2 )
    {

    }

    /**
     ** Sample image at a specified position
     ** @param src Input image
     ** @param y Y-coordinate of sampling
     ** @param x X-coordinate of sampling
     ** @return Sampled value
     **/
    template <typename T>
    T operator()( const Image<T> & src , const float y , const float x ) const
    {
      const int im_width = src.Width();
      const int im_height = src.Height();

      // Get sampler coefficients
      double coefs_x[ SamplerFunc::neighbor_width ];
      double coefs_y[ SamplerFunc::neighbor_width ];

      // Compute difference between exact pixel location and sample
      const double dx = static_cast<double>( x ) - floor( x );
      const double dy = static_cast<double>( y ) - floor( y );

      // Get sampler weights
      sampler_( dx , coefs_x );
      sampler_( dy , coefs_y );

      typename RealPixel<T>::real_type res( 0 );

      // integer position of sample (x,y)
      const int grid_x = static_cast<int>( floor( x ) );
      const int grid_y = static_cast<int>( floor( y ) );

      // Sample a grid around specified grid point
      double total_weight = 0.0;
      for (int i = 0; i < SamplerFunc::neighbor_width; ++i )
      {
        // Get current i value
        // +1 for correct scheme (draw it to be conviced)
        const int cur_i = grid_y + 1 + i - half_width_;

        // handle out of range
        if (cur_i < 0 || cur_i >= im_height )
        {
          continue;
        }

        for (int j = 0; j < SamplerFunc::neighbor_width; ++j )
        {
          // Get current j value
          // +1 for the same reason
          const int cur_j = grid_x + 1 + j - half_width_;

          // handle out of range
          if (cur_j < 0 || cur_j >= im_width )
          {
            continue;
          }


          // sample input image and weight according to sampler
          const double w = coefs_x[ j ] * coefs_y[ i ];
          const typename RealPixel<T>::real_type pix = RealPixel<T>::convert_to_real( src( cur_i , cur_j ) );
          const typename RealPixel<T>::real_type wp = pix * w;
          res += wp;

          total_weight += w;
        }
      }

      // If value too small, it should be so instable, so return the sampled value
      if (total_weight <= 0.2 )
      {
        return T();
      }

      if (total_weight != 1.0 )
      {
        res /= total_weight;
      }


      return RealPixel<T>::convert_from_real( res );
    }

  private:

    /// Sampler function used to resample input image
    SamplerFunc sampler_;

    /// Sampling window
    const int half_width_;
};

} // namespace image
} // namespace openMVG

#endif  // OPENMVG_IMAGE_SAMPLE_HPP
