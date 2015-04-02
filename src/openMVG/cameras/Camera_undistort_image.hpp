
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_UNDISTORT_IMAGE_HPP
#define OPENMVG_CAMERA_UNDISTORT_IMAGE_HPP

#include "openMVG/image/image.hpp"

namespace openMVG{

/// Undistort an image according a given camera & it's distortion model
template <typename Image>
void UndistortImage(
  const Image& imageIn,
  const IntrinsicBase * cam,
  Image & image_ud,
  typename Image::Tpixel fillcolor = typename Image::Tpixel(0))
{
  if (!cam->have_disto()) // no distortion, perform a direct copy
  {
    image_ud = imageIn;
  }
  else // There is distortion
  {
    image_ud.resize(imageIn.Width(), imageIn.Height(), true, fillcolor);
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int j = 0; j < imageIn.Height(); ++j)
    for (int i = 0; i < imageIn.Width(); ++i)
    {
      const Vec2 distorted(i,j);
      // compute coordinates with distortion
      const Vec2 undistorted = cam->cam2ima( cam->add_disto(cam->ima2cam(distorted)) );
      // pick pixel if it is in the image domain
      if ( imageIn.Contains(undistorted(1), undistorted(0)) )
        image_ud( j, i ) = SampleLinear(imageIn, undistorted(1), undistorted(0));
    }
  }
}

} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERA_UNDISTORT_IMAGE_HPP

