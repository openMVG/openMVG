// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONCAT_H_
#define OPENMVG_IMAGE_IMAGE_CONCAT_H_

#include "openMVG/image/image_container.hpp"

/// Horizontal concatenation of images
template < class Image >
void ConcatH(const Image & imageA, const Image & imageB, Image & Out)
{
  // Compute new dimensions.
  int ww = imageA.Width() + imageB.Width();

  Out = Image(ww, std::max(imageA.Height(), imageB.Height()));

  // Fill with original data from imageA.
  for(size_t i = 0; i < imageA.Width(); ++i)
    for(size_t j = 0; j < imageA.Height(); ++j)
      Out(j,i) = imageA(j,i);

  // Fill with original data from imageB with the imageA Width offset.
  const size_t offset = imageA.Width();
  for(size_t i = 0; i < imageB.Width(); ++i)
    for(size_t j = 0; j < imageB.Height(); ++j)
      Out(j,i+offset) = imageB(j,i);
}

/// Vertical concatenation of images
template < class Image >
void ConcatV(const Image & imageA, const Image & imageB, Image & Out)
{
  // Compute new dimensions.
  int hh = imageA.Height() + imageB.Height();

  Out = Image(max(imageA.Width(), imageB.Width()), hh);

  // Fill with original data from imageA.
  for(size_t i = 0; i < imageA.Width(); ++i)
    for(size_t j = 0; j < imageA.Height(); ++j)
      Out(j,i) = imageA(j,i);

  // Fill with original data from imageB with the imageA Height offset.
  const size_t offset = imageA.Height();
  for(size_t i = 0; i < imageB.Width(); ++i)
    for(size_t j = 0; j < imageB.Height(); ++j)
      Out(j+offset,i) = imageB(j,i);
}

#endif // OPENMVG_IMAGE_IMAGE_CONCAT_H_
