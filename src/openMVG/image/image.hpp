// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_HPP
#define OPENMVG_IMAGE_IMAGE_HPP

// Get rid of the specific MSVC compiler warnings.
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
# define _CRT_SECURE_NO_WARNINGS
#endif


#include "openMVG/numeric/numeric.h"

#include "openMVG/image/image_concat.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/image_converter.hpp"
#include "openMVG/image/image_convolution.hpp"
#include "openMVG/image/image_convolution_base.hpp"
#include "openMVG/image/image_diffusion.hpp"
#include "openMVG/image/image_drawing.hpp"
#include "openMVG/image/image_filtering.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_resampling.hpp"
#include "openMVG/image/pixel_types.hpp"
#include "openMVG/image/sample.hpp"


namespace openMVG
{
  /**
  * @brief namespace containing various image manipulation classes and functions
  */
  namespace image
  {

    template class Image<unsigned char>;
    template class Image<float>;
    template class Image<double>;
    template class Image<RGBColor>;
    template class Image<RGBAColor>;

  }
}


#endif /* OPENMVG_IMAGE_IMAGE_HPP */
