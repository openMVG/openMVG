// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_INTEGRAL_HPP
#define OPENMVG_IMAGE_IMAGE_INTEGRAL_HPP

#include "openMVG/image/image_container.hpp"

#include <iostream>
namespace openMVG
{
namespace image
{

/**
* @brief Create a summed area table:
* -  http://en.wikipedia.org/wiki/Summed_Area_Table
* - "Summed-area tables for texture mapping". SIGGRAPH84. Crow, Franklin.
* @param[in] image Input image.
* @param[out] columns Output integral image.
* note: in order to avoid sum overflow please choose carefully the
* TIntegralImage::Tpixel of the integral image. It must be larger than the
* input bit depth of the input image.
* => i.e:
* (uint8_t input -< uint32_t integral image) => only work up to a 4K image
* If we have a fully staurated image we can handle only the following width, height:
* std::sqrt(std::numeric_limits<uint32_t>::max()/255) => 4104
* If you use a larger image you must use a Image<uint64_t> integral_image.
**/
template <typename TImage, typename TIntegralImage>
inline void IntegralImage
(
  const TImage &image,
  TIntegralImage *integral_image
)
{
  typedef typename TIntegralImage::Tpixel Tpixel;
  integral_image->resize(image.Width(), image.Height());

  // Split first row from the rest to avoid an if in the inner loop.
  Tpixel row_sum = Tpixel(0);
  for (int c = 0; c < image.Width(); ++c)
  {
    row_sum += Tpixel(image(0, c));
    (*integral_image)(0, c) = row_sum;
  }

  // Each pixel is a sum of all the pixels to the left and up of that pixel.
  for (int r = 1; r < image.Height(); ++r)
  {
    row_sum = Tpixel(0);
    for (int c = 0; c < image.Width(); ++c)
    {
      row_sum += Tpixel(image(r, c));
      (*integral_image)(r, c) = row_sum + (*integral_image)(r - 1, c);
    }
  }
}

/**
* @brief Compute the sum of the pixels in a rectangular region:
* -  http://en.wikipedia.org/wiki/Summed_Area_Table
* - "Summed-area tables for texture mapping". SIGGRAPH84. Crow, Franklin.
* @param[in] left_corner Input top left corner that defines the rectangular region.
* @param[in] right_corner Input bottom right corner that defines the rectangular region.
* @param[in] integral_image Output integral image.
* @return The sum of pixel intensitied in the depicted area.
**/
template <typename TIntegralImage>
typename TIntegralImage::Tpixel BlockSum
(
  const Eigen::Vector2i & left_corner,
  const Eigen::Vector2i & right_corner,
  const TIntegralImage & integral_image
)
{
  typedef typename TIntegralImage::Tpixel Tpixel;
  assert(right_corner.x() > left_corner.x());
  assert(right_corner.y() > left_corner.y());
  // A ---- B
  // |      |
  // C----- D
  // Sum in the area is equal to: D + A − (C + B).
  const auto & A = integral_image(left_corner.y() - 1, left_corner.x() - 1);
  const auto & B = integral_image(left_corner.y() - 1, right_corner.x());
  const auto & C = integral_image(right_corner.y()   , left_corner.x() - 1);
  const auto & D = integral_image(right_corner.y()   , right_corner.x());
  return (D + A - C - B);
}

} // namespace image
} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_INTEGRAL_HPP
