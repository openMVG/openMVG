// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVERTER_HPP
#define OPENMVG_IMAGE_IMAGE_CONVERTER_HPP

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"

namespace openMVG
{
namespace image
{

/**
* @brief Convert from one color type to another one
* @param valin Input color value
* @param out Output color value after conversion
* @tparam Tin Input color type
* @tparam Tout Output color type
*/
template<typename Tin, typename Tout>
inline void Convert( const Tin& valin, Tout& out )
{
  out = static_cast<Tout>( valin );
}

/**
* @brief Convert RGBColor to unsigned char color
* @param valin Input color value
* @param valOut Output color value after conversion
*/
template<>
inline void Convert<RGBColor, unsigned char>(
  const RGBColor& valin, unsigned char& valOut )
{
  valOut = static_cast<unsigned char>( 0.3 * valin.r() + 0.59 * valin.g() + 0.11 * valin.b() );
}

/**
* @brief Convert RGBAColor to unsigned char color
* @param valin Input color value
* @param valOut Output color value after conversion
*/
template<>
inline void Convert<RGBAColor, unsigned char>(
  const RGBAColor& valin, unsigned char& valOut )
{
  valOut = static_cast<unsigned char>(
             ( valin.a() / 255.f ) *
             ( 0.3 * valin.r() + 0.59 * valin.g() + 0.11 * valin.b() ) );
}

/**
* @brief Convert RGBAColor to RGBColor
* @param valin Input color value
* @param valOut Output color value after conversion
*/
template<>
inline void Convert<RGBAColor, RGBColor>(
  const RGBAColor& valin, RGBColor& valOut )
{
  valOut = RGBColor(
             static_cast<unsigned char> ( ( valin.a() / 255.f ) * valin.r() ),
             static_cast<unsigned char> ( ( valin.a() / 255.f ) * valin.g() ),
             static_cast<unsigned char> ( ( valin.a() / 255.f ) * valin.b() ) );
}

/**
* @brief Convert image from one internal pixel type to another one
* @param imaIn Input image
* @param imaOut Output image
* @tparam ImageIn Type of input image
* @tparam ImageOut Type of output image
*/
template<typename ImageIn, typename ImageOut>
void ConvertPixelType( const ImageIn& imaIn, ImageOut *imaOut )
{
  ( *imaOut ) = ImageOut( imaIn.Width(), imaIn.Height() );
  // Convert each input pixel to destination pixel
  for (int j = 0; j < imaIn.Height(); ++j )
    for (int i = 0; i < imaIn.Width(); ++i )
    {
      Convert( imaIn( j, i ), ( *imaOut )( j, i ) );
    }
}

//--------------------------------------------------------------------------
// RGB ( unsigned char or int ) to Float
//--------------------------------------------------------------------------

/**
* @brief Convert RGB color stored as an int to RGB color stored as a float
* @param valIn Input color
* @param valOut Output color
* @param scaling factor applied to input color component
* @tparam Tin Input color type
* @tparam[out] Tout Output color type
* @todo Use SFINAE to ensure input type is an intergral one
* @todo Why not using RGBfColor as output type ?
*/
template< typename Tin, typename Tout >
inline void convertRGB2Float(
  const Tin& valIn,
  Tout& valOut,
  float factor = 1.0f / 255.f )
{
  for (int channel = 0; channel < 3; ++channel )
  {
    valOut( channel ) = ( float )( ( int )( valIn( channel ) ) * factor );
  }
}

/**
* @brief Convert RGB image stored as an int to RGB image stored as a float
* @param imaIn Input image
* @param[out] imaOut Output image
* @param scaling factor applied to input image color component
* @tparam ImageIn Input image type
* @todo Use SFINAE to ensure input type is an intergral one
*/
template< typename ImageIn >
void rgb2Float( const ImageIn& imaIn,
                Image< RGBfColor > *imaOut, float factor = 1.0f / 255.f )
{
  assert( imaIn.Depth() == 3 );
  ( *imaOut ).resize( imaIn.Width(), imaIn.Height() );
  // Convert each int RGB to float RGB values
  for (int j = 0; j < imaIn.Height(); ++j )
    for (int i = 0; i < imaIn.Width(); ++i )
    {
      convertRGB2Float( imaIn( j, i ), ( *imaOut )( j, i ), factor );
    }
}

//--------------------------------------------------------------------------
// Float to RGB ( unsigned char or int )
//--------------------------------------------------------------------------

/**
* @brief Convert float RGB color to Int RGB color
* @param valin Input color
* @param[out] Output color
* @param factor scaling factor applied to input color components
*/
inline
void convertFloatToInt
(
  const RGBfColor& valIn,
  RGBColor& valOut,
  float factor = 255.f
)
{
  for (int channel = 0; channel < 3; ++channel )
  {
    valOut( channel ) = ( int )( valIn( channel ) * factor );
  }
}

/**
* @brief Convert RGB image stored as float components to RGB image stored as int components
* @param imaIn Input image
* @param[out] imaOut Output image
* @param factor scaling factor applied to each input color component
*/
inline void rgbFloat2rgbInt(
  const Image< RGBfColor >& imaIn,
  Image< RGBColor > *imaOut,
  float factor = 255.f )
{
  assert( imaIn.Depth() == 3 );
  ( *imaOut ).resize( imaIn.Width(), imaIn.Height() );
  // Convert each int RGB to float RGB values
  for (int j = 0; j < imaIn.Height(); ++j )
  {
    for (int i = 0; i < imaIn.Width(); ++i )
    {
      convertFloatToInt( imaIn( j, i ), ( *imaOut )( j, i ), factor  );
    }
  }
}

} // namespace image
} // namespace openMVG

#endif  // OPENMVG_IMAGE_IMAGE_CONVERTER_HPP
