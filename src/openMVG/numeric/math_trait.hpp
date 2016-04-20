// Copyright (c) 2014 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef _OPENMVG_NUMERIC_MATH_TRAIT_HPP_
#define _OPENMVG_NUMERIC_MATH_TRAIT_HPP_

#include <cmath>

namespace openMVG
{


  /**
  * @brief Type trait used to specialize math functions
  * @tparam T type used for computation 
  */
  template <typename T>
  class MathTrait
  {
  public:

    // Trigonometric functions
    static inline T cos( const T val ) ;
    static inline T sin( const T val ) ;
    static inline T tan( const T val ) ;

    // Inverse trigonometric functions
    static inline T acos( const T val ) ;
    static inline T asin( const T val ) ;
    static inline T atan( const T val ) ;
    static inline T atan2( const T y , const T x ) ;

    // Exponential functions
    static inline T exp( const T val ) ;
    static inline T log( const T val ) ;
    static inline T log10( const T val ) ;

    // Power functions
    static inline T pow( const T base , const T expo ) ;
    static inline T sqrt( const T val ) ;
    static inline T cbrt( const T val ) ;

    // Rounding functions
    static inline T floor( const T val ) ;
    static inline T ceil( const T val ) ;
    static inline T round( const T val ) ;

    // Absolute value
    static inline T abs( const T val ) ;
  } ;


  /**
   ** Float specialization
   **/
  template<>
  inline float MathTrait<float>::cos( const float val )
  {
    return cosf( val ) ;
  }

  template<>
  inline float MathTrait<float>::sin( const float val )
  {
    return sinf( val ) ;
  }

  template<>
  inline float MathTrait<float>::tan( const float val )
  {
    return tanf( val ) ;
  }

  template<>
  inline float MathTrait<float>::acos( const float val )
  {
    return acosf( val ) ;
  }

  template<>
  inline float MathTrait<float>::asin( const float val )
  {
    return asinf( val ) ;
  }

  template<>
  inline float MathTrait<float>::atan( const float val )
  {
    return atanf( val ) ;
  }

  template<>
  inline float MathTrait<float>::atan2( const float y , const float x )
  {
    return atan2f( y , x ) ;
  }

  template<>
  inline float MathTrait<float>::exp( const float val )
  {
    return expf( val ) ;
  }

  template<>
  inline float MathTrait<float>::log( const float val )
  {
    return logf( val ) ;
  }

  template<>
  inline float MathTrait<float>::log10( const float val )
  {
    return log10f( val ) ;
  }

  template<>
  inline float MathTrait<float>::pow( const float base , const float expo )
  {
    return powf( base , expo ) ;
  }

  template<>
  inline float MathTrait<float>::sqrt( const float val )
  {
    return sqrtf( val ) ;
  }

  template<>
  inline float MathTrait<float>::cbrt( const float val )
  {
#ifdef _MSC_VER
    return cbrt( val ) ;
#else
    return cbrtf( val ) ;
#endif
  }

  template<>
  inline float MathTrait<float>::floor( const float val )
  {
    return floorf( val ) ;
  }

  template<>
  inline float MathTrait<float>::ceil( const float val )
  {
    return ceilf( val ) ;
  }

  template<>
  inline float MathTrait<float>::round( const float val )
  {
    if( val >= 0.0f )
    {
      return floorf( val + 0.5f ) ;
    }
    else
    {
      return ceilf( val - 0.5f ) ;
    }
  }

  template<>
  inline float MathTrait<float>::abs( const float val )
  {
    return fabsf( val ) ;
  }


  /**
   ** Double specialization
  **/
  template<>
  inline double MathTrait<double>::cos( const double val )
  {
    return cos( val ) ;
  }

  template<>
  inline double MathTrait<double>::sin( const double val )
  {
    return sin( val ) ;
  }

  template<>
  inline double MathTrait<double>::tan( const double val )
  {
    return tan( val ) ;
  }

  template<>
  inline double MathTrait<double>::acos( const double val )
  {
    return acos( val ) ;
  }

  template<>
  inline double MathTrait<double>::asin( const double val )
  {
    return asin( val ) ;
  }

  template<>
  inline double MathTrait<double>::atan( const double val )
  {
    return atan( val ) ;
  }

  template<>
  inline double MathTrait<double>::atan2( const double y , const double x )
  {
    return atan2( y , x ) ;
  }

  template<>
  inline double MathTrait<double>::exp( const double val )
  {
    return exp( val ) ;
  }

  template<>
  inline double MathTrait<double>::log( const double val )
  {
    return log( val ) ;
  }

  template<>
  inline double MathTrait<double>::log10( const double val )
  {
    return log10( val ) ;
  }

  template<>
  inline double MathTrait<double>::pow( const double base , const double expo )
  {
    return pow( base , expo ) ;
  }

  template<>
  inline double MathTrait<double>::sqrt( const double val )
  {
    return sqrt( val ) ;
  }

  template<>
  inline double MathTrait<double>::cbrt( const double val )
  {
    return cbrt( val ) ;
  }

  template<>
  inline double MathTrait<double>::floor( const double val )
  {
    return floor( val ) ;
  }

  template<>
  inline double MathTrait<double>::ceil( const double val )
  {
    return ceil( val ) ;
  }

  template<>
  inline double MathTrait<double>::round( const double val )
  {
    if( val >= 0.0 )
    {
      return floor( val + 0.5 ) ;
    }
    else
    {
      return ceil( val - 0.5 ) ;
    }
  }

  template<>
  inline double MathTrait<double>::abs( const double val )
  {
    return fabs( val ) ;
  }


  /**
   ** long double specialization
  **/
  template<>
  inline long double MathTrait<long double>::cos( const long double val )
  {
    return cosl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::sin( const long double val )
  {
    return sinl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::tan( const long double val )
  {
    return tanl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::acos( const long double val )
  {
    return acosl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::asin( const long double val )
  {
    return asinl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::atan( const long double val )
  {
    return atanl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::atan2( const long double y , const long double x )
  {
    return atan2l( y , x ) ;
  }

  template<>
  inline long double MathTrait<long double>::exp( const long double val )
  {
    return expl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::log( const long double val )
  {
    return logl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::log10( const long double val )
  {
    return log10l( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::pow( const long double base , const long double expo )
  {
    return powl( base , expo ) ;
  }

  template<>
  inline long double MathTrait<long double>::sqrt( const long double val )
  {
    return sqrtl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::cbrt( const long double val )
  {
#ifdef _MSC_VER
    return cbrt( val ) ;
#else
    return cbrtl( val ) ;
#endif
  }

  template<>
  inline long double MathTrait<long double>::floor( const long double val )
  {
    return floorl( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::ceil( const long double val )
  {
    return ceill( val ) ;
  }

  template<>
  inline long double MathTrait<long double>::round( const long double val )
  {
    if( val >= 0.0l )
    {
      return floorl( val + 0.5l ) ;
    }
    else
    {
      return ceill( val - 0.5l ) ;
    }
  }

  template<>
  inline long double MathTrait<long double>::abs( const long double val )
  {
    return fabsl( val ) ;
  }
} // namespace openMVG 

#endif