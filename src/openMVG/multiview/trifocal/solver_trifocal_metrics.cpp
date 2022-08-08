// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\author Pierre MOULON
//\author Gabriel ANDRADE Rio de Janeiro State U.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"
#include "openMVG/multiview/triangulation.hpp"

namespace openMVG {
namespace trifocal {
  
double NormalizedSquaredPointReprojectionOntoOneViewError::
Error(
  const trifocal_model_t &tt,
  const Vec &bearing_0, // x,y,tangentialx,tangentialy
  const Vec &bearing_1,
  const Vec &bearing_2) 
{
  // Return the cost related to this model and those sample data point
  // Ideal algorithm:
  // 1) reconstruct the 3D points and orientations
  // 2) project the 3D points and orientations on all images_
  // 3) compute error 
  // 
  // In practice we ignore the directions and only reproject to one third view
  // 3x3: each column is x,y,1
  Mat3 bearing;
  bearing << bearing_0.head(2).homogeneous(),
             bearing_1.head(2).homogeneous(), 
             bearing_2.head(2).homogeneous();
  
  Vec4 triangulated_homg;
  unsigned third_view = 0;
  // pick the wider baseline. TODO: measure all pairwise translation distances
  if (tt[1].col(3).squaredNorm() > tt[2].col(3).squaredNorm()) {
    // TODO(trifocal future) compare to triangulation from the three views at once
    TriangulateDLT(tt[0], bearing.col(0), tt[1], bearing.col(1), &triangulated_homg);
    third_view = 2;
  } else {
    TriangulateDLT(tt[0], bearing.col(0), tt[2], bearing.col(2), &triangulated_homg);
    third_view = 1;
  }

  // Computing the projection of triangulated points using projection.hpp
  // For prototyping and speed, for now we will only project to the third view
  // and report only one error
  // TODO: it is a good idea to filter the inliers after a robust estimation
  // using a more complete (and heavier) error metric.
  Vec2 p_reprojected = (tt[third_view]*triangulated_homg).hnormalized();
  return (p_reprojected - bearing.col(third_view).head(2)).squaredNorm();
}

} // namespace trifocal
} // namespace OpenMVG
