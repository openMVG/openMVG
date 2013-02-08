// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVERTER_HPP
#define OPENMVG_IMAGE_IMAGE_CONVERTER_HPP

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"

namespace openMVG{

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

template<typename ImageIn, typename ImageOut>
void Rgb2Gray(const ImageIn& imaIn, ImageOut *imaOut)
{
  assert( imaIn.Depth() == 3 );
  (*imaOut) = ImageOut(imaIn.Width(), imaIn.Height());
  // Convert each RGB pixel into Grey value (luminance)
  for(int j = 0; j < imaIn.Height(); ++j)
    for(int i = 0; i < imaIn.Width(); ++i)
      Convert(imaIn(j,i), (*imaOut)(j,i));
}

} // namespace openMVG

#endif  // OPENMVG_IMAGE_IMAGE_CONVERTER_HPP
