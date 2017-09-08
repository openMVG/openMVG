// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_UNDISTORT_IMAGE_HPP
#define OPENMVG_CAMERAS_CAMERA_UNDISTORT_IMAGE_HPP

#include <limits>
#include <cmath>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"

namespace openMVG
{
namespace cameras
{

/**
* @brief  Undistort an image according a given camera & its distortion model
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

/**
* @brief  Undistort an image according a given camera & its distortion model
* @param imageIn Input image
* @param cam Input intrinsic parameter used to undistort image
* @param[out] image_ud Output undistorted image
* @param fillcolor color used to fill pixels where no input pixel is found
* @param max_ud_width Maximum width of the undistorted image
* @param max_ud_height Maximum height of the undistorted image
* @note This function produces an image with size that try to fit the image plane
*/
template <typename Image>
void UndistortImageResized(
  const Image& imageIn,
  const IntrinsicBase * cam,
  Image & image_ud,
  typename Image::Tpixel fillcolor = typename Image::Tpixel( 0 ),
  const uint32_t max_ud_width = 10000,
  const uint32_t max_ud_height = 10000)
{
  if ( !cam->have_disto() ) // no distortion, perform a direct copy
  {
    image_ud = imageIn;
  }
  else // There is distortion
  {
    // 1 - Compute size of the Undistorted image
    int min_x = std::numeric_limits<int>::max();
    int min_y = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::lowest();
    int max_y = std::numeric_limits<int>::lowest();

    for (int id_row = 0; id_row < imageIn.Height(); ++id_row )
    {
      for (int id_col = 0; id_col < imageIn.Width(); ++id_col )
      {
        const Vec2 dist_pix( id_col , id_row );
        const Vec2 undist_pix = cam->get_ud_pixel( dist_pix );

        const int x = static_cast<int>( undist_pix[0] );
        const int y = static_cast<int>( undist_pix[1] );

        min_x = std::min( x , min_x );
        min_y = std::min( y , min_y );
        max_x = std::max( x , max_x );
        max_y = std::max( y , max_y );
      }
    }

    // Ensure size is at least 1 pixel (width and height)
    const int computed_size_x = std::max( 1 , max_x - min_x + 1 );
    const int computed_size_y = std::max( 1 , max_y - min_y + 1 );

    // Compute real size (ensure we do not have infinite size)
    const uint32_t real_size_x = std::min( max_ud_width , (uint32_t)computed_size_x );
    const uint32_t real_size_y = std::min( max_ud_height , (uint32_t)computed_size_y );

    // 2 - Compute inverse projection to fill the output image
    image_ud.resize( real_size_x , real_size_y , true, fillcolor );
    const image::Sampler2d<image::SamplerLinear> sampler;

#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
#endif
    for ( int j = 0; j < real_size_y; ++j )
      for ( int i = 0; i < real_size_x; ++i )
      {
        const Vec2 undisto_pix( i + min_x , j + min_y );
        // compute coordinates with distortion
        const Vec2 disto_pix = cam->get_d_pixel( undisto_pix );

        // pick pixel if it is in the image domain
        if ( imageIn.Contains( disto_pix( 1 ), disto_pix( 0 ) ) &&
             (! std::isnan( disto_pix[0] )) &&
             (! std::isnan( disto_pix[1] )) &&
             (! std::isinf( disto_pix[0] )) &&
             (! std::isinf( disto_pix[1] )) )
        {
          image_ud( j, i ) = sampler( imageIn, disto_pix( 1 ), disto_pix( 0 ) );
        }
      }
  }
}


} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_UNDISTORT_IMAGE_HPP
