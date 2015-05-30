// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_SAMPLE_HPP
#define OPENMVG_IMAGE_SAMPLE_HPP

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/pixel_types.hpp"

namespace openMVG {
namespace image {

/// Nearest neighbor interpolation.
template<typename T>
inline T sampleNearest(const Image<T>& image,
                       float y, float x) {
  return image(static_cast<int>(round(y)), static_cast<int>(round(x)));
}

static inline void LinearInitAxis(float fx, int width,
  int *x1, int *x2,
  float *dx1, float *dx2) {

  const int ix = static_cast<int>(fx);
  if (ix < 0) {
    *x1 = *x2 = 0;
    *dx1 = 1;
    *dx2 = 0;
  } else if (ix > width-2) {
    *x1 = *x2 = width-1;
    *dx1 = 1;
    *dx2 = 0;
  } else {
    *x1 = ix;
    *x2 = *x1 + 1;
    *dx1 = *x2 - fx;
    *dx2 = 1 - *dx1;
  }
}

/// Linear interpolation.
template<typename T>
inline T SampleLinear(const Image<T>& image, float y, float x) {
  int x1, y1, x2, y2;
  float dx1, dy1, dx2, dy2;

  LinearInitAxis(y, image.Height(), &y1, &y2, &dy1, &dy2);
  LinearInitAxis(x, image.Width(),  &x1, &x2, &dx1, &dx2);

  const T im11 = image(y1, x1);
  const T im12 = image(y1, x2);
  const T im21 = image(y2, x1);
  const T im22 = image(y2, x2);

  return T(( im11 * dx1 + im12 * dx2) * dy1 +
           ( im21 * dx1 + im22 * dx2) * dy2 );
}

/// Linear interpolation.
/// RGBColor specialization.
template<>
inline RGBColor SampleLinear<RGBColor>(const Image<RGBColor>& image, float y, float x) {
  int x1, y1, x2, y2;
  float dx1, dy1, dx2, dy2;

  LinearInitAxis(y, image.Height(), &y1, &y2, &dy1, &dy2);
  LinearInitAxis(x, image.Width(),  &x1, &x2, &dx1, &dx2);

  const Rgb<float> im11((image(y1, x1).cast<float>()));
  const Rgb<float> im12((image(y1, x2).cast<float>()));
  const Rgb<float> im21((image(y2, x1).cast<float>()));
  const Rgb<float> im22((image(y2, x2).cast<float>()));

  return RGBColor(
      (
        ( im11 * dx1 + im12 * dx2) * dy1 +
        ( im21 * dx1 + im22 * dx2) * dy2
      )
      .cast<unsigned char>()
    );
}

} // namespace image
} // namespace openMVG

#endif  // OPENMVG_IMAGE_SAMPLE_HPP
