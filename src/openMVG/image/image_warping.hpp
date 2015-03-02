
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_HOMOGRAPHY_WARP
#define OPENMVG_IMAGE_HOMOGRAPHY_WARP

namespace openMVG{
namespace image{

/// Apply inplace homography transform for the given point (x,y).
/// Return true if H is orientation preserving around the point.
bool ApplyH_AndCheckOrientation(const Mat3 &H, double &x, double &y)
{
  Vec3 X(x, y, 1.0);
  X = H*X;
  X /= X(2);
  x = X(0);
  y = X(1);
  return (X(2) * H(2,2) > 0.0);
}

/// Warp an image im given a homography H with a backward approach
/// H must be already have been resized accordingly
template <class Image>
void Warp(const Image &im, const Mat3 & H, Image &out)
{
  const int wOut = static_cast<int>(out.Width());
  const int hOut = static_cast<int>(out.Height());

  for (int j = 0; j < hOut; ++j)
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for
#endif
    for (int i = 0; i < wOut; ++i)
    {
      double xT = i, yT = j;
      if (ApplyH_AndCheckOrientation(H, xT, yT)
          && im.Contains(yT,xT))
        out(j,i) = SampleLinear(im, (float)yT, (float)xT);
    }
}

}; // namespace image
}; // namespace openMVG

#endif // OPENMVG_IMAGE_HOMOGRAPHY_WARP
