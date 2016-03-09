// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_PIXELTYPES_HPP
#define OPENMVG_IMAGE_PIXELTYPES_HPP

#include "openMVG/numeric/numeric.h"

namespace openMVG
{
namespace image
{

/**
 * Red Greeen Blue template pixel type
 * @tparam T type of each channel
 */
template <typename T>
class Rgb : public Eigen::Matrix<T, 3, 1, 0, 3, 1>
{
    /// Full internal type
    typedef Eigen::Matrix<T, 3, 1, 0, 3, 1> Base;

    /// Color component type
    typedef T TBase;
  public:

    //------------------------------
    //-- construction method
    /**
    * @brief Full constructor
    * @param red Red value
    * @param green value
    * @param blue value
    */
    inline Rgb( T red, T green, T blue )
      : Base( red, green, blue )
    {

    }

    /**
    * @brief Copy constructor
    * @param val Source RGB value
    */
    explicit inline Rgb( const Base& val )
      : Base( val )
    {

    }

    /**
    * @brief Single value construction
    * @param val Value to set in each channel
    * @note This is equivalent to Rgb( val , val , val )
    */
    explicit inline Rgb( const T val = 0 )
      : Base( val, val, val )
    {

    }
    //-- construction method
    //------------------------------

    //------------------------------
    //-- accessors/getters methods
    /**
    * @brief Get readonly Red channel value
    * @return Red channel value
    */
    inline const T& r() const
    {
      return ( *this )( 0 );
    }

    /**
    * @brief Get rw Red channel value
    * @return Red channel value
    */
    inline T& r()
    {
      return ( *this )( 0 );
    }

    /**
    * @brief Get readonly Green channel value
    * @return Green channel value
    */
    inline const T& g() const
    {
      return ( *this )( 1 );
    }

    /**
    * @brief Get rw Green channel value
    * @return Green channel value
    */
    inline T& g()
    {
      return ( *this )( 1 );
    }

    /**
    * @brief Get readonly Blue channel value
    * @return Blue channel value
    */
    inline const T& b() const
    {
      return ( *this )( 2 );
    }

    /**
    * @brief Get rw Blue channel value
    * @return Blue channel value
    */
    inline T& b()
    {
      return ( *this )( 2 );
    }

    /**
    * @brief Get approximate Gray value using rec601 Luma conversion
    * @return Luminance value of pixel
    */
    inline operator T() const
    {
      return T( 0.299 * r() + 0.587 * g() + 0.114 * b() );
    }

    //-- accessors/getters methods
    //------------------------------

    /**
    * @brief stream operator
    * @param os Stream in which rgb value is outputed
    * @param col Color to store into the stream
    * @return stream after output
    */
    friend std::ostream& operator<<( std::ostream& os, const Rgb& col )
    {
      os << " {" ;
      for( int i = 0; i < 2; ++i )
      {
        os << col( i ) << ",";
      }
      os << col( 2 ) << "} ";
      return os;
    }

    /**
    * @brief scalar division
    * @param val Scalar divisor factor
    * @return Rgb color after scalar division
    * @note This does not modify the Rgb value (ie: only return a modified copy)
    */
    template <typename Z>
    inline Rgb operator /( const Z& val ) const
    {
      return Rgb( T( ( Z )( ( *this )( 0 ) ) / val ),
                  T( ( Z )( ( *this )( 1 ) ) / val ),
                  T( ( Z )( ( *this )( 2 ) ) / val ) );
    }

    /**
    * @brief scalar multiplication
    * @param val Scale multiplication factor
    * @return Rgb color after scalar multiplication
    * @note This does not modify the Rgb value (ie: only return a modified copy)
    */
    template <typename Z>
    inline Rgb operator *( const Z& val ) const
    {
      return Rgb( T( ( Z )( *this )( 0 ) * val ),
                  T( ( Z )( *this )( 1 ) * val ),
                  T( ( Z )( *this )( 2 ) * val ) );
    }
};

/// Instantiation for unsigned char color component
typedef Rgb<unsigned char> RGBColor;
/// Instantiation for float color component
typedef Rgb<float> RGBfColor;

/**
* @brief RGBA templated pixel type
*/
template <typename T>
class Rgba : public Eigen::Matrix<T, 4, 1, 0, 4, 1>
{
    /// Full internal type
    typedef Eigen::Matrix<T, 4, 1, 0, 4, 1> Base;

    /// Color component type
    typedef T TBase;

  public:

    //------------------------------
    //-- construction method
    /**
    * @brief Full constructor
    * @param red Red component value
    * @param green Green component value
    * @param blue component value
    * @param alpha component value
    */
    inline Rgba( const T red, const T green, const T blue, const T alpha = static_cast<T>( 1 ) )
      : Base( red, green, blue, alpha )
    {

    }

    /**
    * @brief Copy constructor from internal data
    * @param val Source RGBA value
    */
    explicit inline Rgba( const Base& val )
      : Base( val )
    {

    }

    /**
    * @brief RGB constructor with default alpha value
    * @param val Value to set in each RGB component
    * @note This is equivalent to RGBA( val , val , val , 1 )
    */
    explicit inline Rgba( const T val = 0 )
      : Base( val, val, val, static_cast<T>( 1 ) )
    {

    }

    /**
    * @brief Copy constructor
    * @param val Source RGBA value
    */
    inline Rgba( const RGBColor & val )
      : Base( val.r(), val.g(), val.b(), static_cast<T>( 1 ) )
    {

    }
    //-- construction method
    //------------------------------

    //------------------------------
    //-- accessors/getters methods
    /**
    * @brief Get readonly Red channel value
    * @return Red channel value
    */
    inline const T& r() const
    {
      return ( *this )( 0 );
    }
    /**
    * @brief Get rw Red channel value
    * @return Red channel value
    */
    inline T& r()
    {
      return ( *this )( 0 );
    }
    /**
    * @brief Get readonly Green channel value
    * @return Green channel value
    */
    inline const T& g() const
    {
      return ( *this )( 1 );
    }
    /**
    * @brief Get rw Green channel value
    * @return Green channel value
    */
    inline T& g()
    {
      return ( *this )( 1 );
    }
    /**
    * @brief Get readonly Blue channel value
    * @return Blue channel value
    */
    inline const T& b() const
    {
      return ( *this )( 2 );
    }
    /**
    * @brief Get rw Blue channel value
    * @return Blue channel value
    */
    inline T& b()
    {
      return ( *this )( 2 );
    }
    /**
    * @brief Get readonly Alpha channel value
    * @return Alpha channel value
    */
    inline const T& a() const
    {
      return ( *this )( 3 );
    }
    /**
    * @brief Get rw Alpha channel value
    * @return Alpha channel value
    */
    inline T& a()
    {
      return ( *this )( 3 );
    }

    /**
    * @brief Get approximate Gray value using rec601 Luma conversion
    * @return Luminance value of pixel
    */
    inline operator T() const
    {
      return T( 0.299 * r() + 0.587 * g() + 0.114 * b() );
    }
    //-- accessors/getters methods
    //------------------------------
    /**
    * @brief stream operator
    * @param os Stream in which rgb value is outputed
    * @param col Color to store into the stream
    * @return stream after output
    */
    friend std::ostream& operator<<( std::ostream& os, const Rgba& col )
    {
      os << " {" ;
      for( int i = 0; i < 3; ++i )
      {
        os << col( i ) << ",";
      }
      os << col( 3 ) << "} ";
      return os;
    }

    /**
    * @brief scalar division
    * @param val Scalar divisor factor
    * @return Rgba color after scalar division
    * @note This does not modify the Rgba value (ie: only return a modified copy)
    */
    template<class Z>
    inline Rgba operator /( const Z& val ) const
    {
      return Rgba( T( ( Z )( *this )( 0 ) / val ),
                   T( ( Z )( *this )( 1 ) / val ),
                   T( ( Z )( *this )( 2 ) / val ),
                   T( ( Z )( *this )( 3 ) / val ) );
    }

    /**
    * @brief scalar multiplication
    * @param val Scale multiplication factor
    * @return Rgba color after scalar multiplication
    * @note This does not modify the Rgba value (ie: only return a modified copy)
    */
    template<class Z>
    inline Rgba operator *( const Z& val ) const
    {
      return Rgba( T( ( Z )( *this )( 0 ) * val ),
                   T( ( Z )( *this )( 1 ) * val ),
                   T( ( Z )( *this )( 2 ) * val ),
                   T( ( Z )( *this )( 3 ) * val ) );
    }
};

/// Type used to handle RGBA color in unsigned char format for each component
typedef Rgba<unsigned char> RGBAColor;

const RGBColor WHITE( 255, 255, 255 );
const RGBColor BLACK( 0, 0, 0 );
const RGBColor BLUE( 0, 0, 255 );
const RGBColor RED( 255, 0, 0 );
const RGBColor GREEN( 0, 255, 0 );
const RGBColor YELLOW( 255, 255, 0 );
const RGBColor CYAN( 0, 255, 255 );
const RGBColor MAGENTA( 255, 0, 255 );

} // namespace image
} // namespace openMVG

#endif // OPENMVG_IMAGE_PIXELTYPES_HPP

