
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_UNDISTORT_IMAGE_HPP
#define OPENMVG_CAMERA_UNDISTORT_IMAGE_HPP

#include "openMVG/image/image.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"


namespace openMVG
{
namespace cameras
{


/**
* @brief  Undistort an image according a given camera & it's distortion model
* @param imageIn Input image
* @param cam Input intrinsic parameter used to undistort image
* @param[out] image_ud Output undistorted image
* @param fillcolor color used to fill pixels where no input pixel is found
*/
template <typename Image>
void UndistortImage(
  const Image& imageIn,
  const IntrinsicBase * cam,
  Image & image_ud,
  typename Image::Tpixel fillcolor = typename Image::Tpixel( 0 ) )
{
  if ( !cam->have_disto() ) // no distortion, perform a direct copy
  {
    image_ud = imageIn;
  }
  else // There is distortion
  {
    image_ud.resize( imageIn.Width(), imageIn.Height(), true, fillcolor );
    const image::Sampler2d<image::SamplerLinear> sampler;
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
#endif
    for ( int j = 0; j < imageIn.Height(); ++j )
      for ( int i = 0; i < imageIn.Width(); ++i )
      {
        const Vec2 undisto_pix( i, j );
        // compute coordinates with distortion
        const Vec2 disto_pix = cam->get_d_pixel( undisto_pix );
        // pick pixel if it is in the image domain
        if ( imageIn.Contains( disto_pix( 1 ), disto_pix( 0 ) ) )
        {
          image_ud( j, i ) = sampler( imageIn, disto_pix( 1 ), disto_pix( 0 ) );
        }
      }
  }
}

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERA_UNDISTORT_IMAGE_HPP

