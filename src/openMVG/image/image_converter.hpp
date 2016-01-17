// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVERTER_HPP
#define OPENMVG_IMAGE_IMAGE_CONVERTER_HPP

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"

namespace openMVG{
namespace image {

template<typename T>
// The factor comes from http://www.easyrgb.com/
// RGB to XYZ : Y is the luminance channel
// var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722
inline T Rgb2Gray(const T r,const T g, const T b) {
  return r * 0.2126 + g * 0.7152 + b * 0.0722;
}

template<typename Tin, typename Tout>
inline void Convert(const Tin& valin, Tout& out) {
  out = static_cast<Tout>(valin);
}

template<>
inline void Convert<RGBColor, unsigned char>(
  const RGBColor& valin, unsigned char& valOut)
{
  valOut = static_cast<unsigned char>(0.3 * valin.r() + 0.59 * valin.g() + 0.11 * valin.b());
}

template<>
inline void Convert<RGBAColor, unsigned char>(
  const RGBAColor& valin, unsigned char& valOut)
{
  valOut = static_cast<unsigned char>(
    (valin.a()/255.f) *
    (0.3 * valin.r() + 0.59 * valin.g() + 0.11 * valin.b()));
}

template<>
inline void Convert<RGBAColor, RGBColor>(
  const RGBAColor& valin, RGBColor& valOut)
{
  valOut = RGBColor(
    static_cast<unsigned char> ((valin.a()/255.f) * valin.r()),
    static_cast<unsigned char> ((valin.a()/255.f) * valin.g()),
    static_cast<unsigned char> ((valin.a()/255.f) * valin.b()));
}

template<typename ImageIn, typename ImageOut>
void ConvertPixelType(const ImageIn& imaIn, ImageOut *imaOut)
{
  (*imaOut) = ImageOut(imaIn.Width(), imaIn.Height());
  // Convert each input pixel to destination pixel
  for(int j = 0; j < imaIn.Height(); ++j)
    for(int i = 0; i < imaIn.Width(); ++i)
      Convert(imaIn(j,i), (*imaOut)(j,i));
}

//--------------------------------------------------------------------------
// RGB ( unsigned char or int ) to Float
//--------------------------------------------------------------------------

template< typename Tin, typename Tout >
inline void convertRGB2Float(
    const Tin& valIn,
    Tout& valOut,
    float factor = 1.0f / 255.f)
{
  for( int channel = 0; channel < 3; ++channel )
    valOut(channel) = (float)((int)(valIn(channel)) * factor);
}

template< typename ImageIn >
void rgb2Float( const ImageIn& imaIn,
                Image< RGBfColor > *imaOut, float factor = 1.0f / 255.f )
{
  assert( imaIn.Depth() == 3 );
  (*imaOut).resize(imaIn.Width(), imaIn.Height());
  // Convert each int RGB to float RGB values
  for( int j = 0; j < imaIn.Height(); ++j )
    for( int i = 0; i < imaIn.Width(); ++i )
      convertRGB2Float( imaIn( j, i ), ( *imaOut )( j, i ), factor );
}

//--------------------------------------------------------------------------
// Float to RGB ( unsigned char or int )
//--------------------------------------------------------------------------

static inline
void convertFloatToInt
(
  const RGBfColor& valIn,
  RGBColor& valOut,
  float factor = 255.f
)
{
  for( int channel = 0; channel < 3; ++channel )
    valOut(channel) = (int)(valIn(channel) * factor);
}

static void rgbFloat2rgbInt(
        const Image< RGBfColor >& imaIn,
        Image< RGBColor > *imaOut,
        float factor = 255.f )
{
  assert( imaIn.Depth() == 3 );
  (*imaOut).resize(imaIn.Width(), imaIn.Height());
  // Convert each int RGB to float RGB values
  for( int j = 0; j < imaIn.Height(); ++j )
    for( int i = 0; i < imaIn.Width(); ++i )
      convertFloatToInt( imaIn( j, i ), (*imaOut)( j, i ), factor  );
}

} // namespace image
} // namespace openMVG

#endif  // OPENMVG_IMAGE_IMAGE_CONVERTER_HPP
