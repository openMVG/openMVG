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

//: clump arg from minus 1 to 1; also assert abs(arg) not much bigger than 1
inline double clump_to_acos(double x)
{ 
  if (x > 1.0 || x < -1.0) {
    assert(std::fabs(std::fabs(x)-1) < 1e-5);
    if (x > 1.0)
      return 1.0;
    if (x < -1.0)
      return -1.0;
  }
  return x;
}

  
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


// Meant to be run by the 3 points given to trifocal solver
bool  NormalizedSquaredPointReprojectionOntoOneViewError::
Check(
  const trifocal_model_t &tt,
  const Vec &bearing_0, // x,y,tangentialx,tangentialy
  const Vec &bearing_1,
  const Vec &bearing_2) 
{
  // Ideal algorithm:
  // 1) reconstruct the 3D points and orientations
  // 2) compute depths of the 3D points on all views
  // 3) make sure depths are > 1
  // 4) make sure tangent at 3rd point match
  
  Mat3 bearing;
  bearing << bearing_0.head(2).homogeneous(),
             bearing_1.head(2).homogeneous(), 
             bearing_2.head(2).homogeneous();

  Mat3 t;
  t       << bearing_0.tail(2).homogeneous(),
             bearing_1.tail(2).homogeneous(), 
             bearing_2.tail(2).homogeneous();

  t(2,0) = t(2,1) = t(2,2) = 0;
  
  Vec4 triangulated_homg;
  Vec3 Trec;
  unsigned third_view = 0;
  // pick the wider baseline. TODO: measure all pairwise translation distances
  if (tt[1].col(3).squaredNorm() > tt[2].col(3).squaredNorm()) {
    // TODO(trifocal future) compare to triangulation from the three views at once
    TriangulateDLT(tt[0], bearing.col(0), tt[1], bearing.col(1), &triangulated_homg);
//    Trec = t.col(0).cross(bearing.col(0)).cross(t.col(1).cross(bearing.col(1)));
    third_view = 2;
  } else {
    TriangulateDLT(tt[0], bearing.col(0), tt[2], bearing.col(2), &triangulated_homg);
//    Trec = t.col(0).cross(bearing.col(0)).cross(t.col(2).cross(bearing.col(2)));
    third_view = 1;
  }

  // Computing the projection of triangulated points using projection.hpp
  // For prototyping and speed, for now we will only project to the third view
  // and report only one error
  // TODO: it is a good idea to filter the inliers after a robust estimation
  // using a more complete (and heavier) error metric.
  Vec3 p_third_view = tt[third_view]*triangulated_homg;
  Vec2 p_reprojected = p_third_view.hnormalized();

  OPENMVG_LOG_INFO << "P reproj " << p_reprojected;
  OPENMVG_LOG_INFO << "P third " <<bearing.col(third_view).head(2);
  double err = (p_reprojected - bearing.col(third_view).head(2)).squaredNorm();
  OPENMVG_LOG_INFO <<  "Solver 3rd point sq reprojection error: " << err << std::endl;
  if (err > 1e-3)
    return false;

  Vec3 p_second_view = tt[(third_view == 1)?2:1] * triangulated_homg;

  if (triangulated_homg.hnormalized()(2) <= 0. || p_third_view(2) <= 0. || p_second_view(2) <= 0.) {
    OPENMVG_LOG_INFO << "Internal Cheirality check FAIL" << std::endl;
    OPENMVG_LOG_INFO <<  triangulated_homg.hnormalized()(2)  << " , " <<  p_third_view(2) << " , " << p_second_view(2);
    return false;
  }
  OPENMVG_LOG_INFO << "Internal Cheirality check PASS" << std::endl;
  
#if 0 // code below is working, just have to check
  tproj = (Trec - Trec.dot(tt[third_view].(2,seq(0,2))/* z axis of third camera relative to 1st, last line of R*/)*bearing.col(third_view));

  // compute angle between tproj t.col(3)
  double angular_error = std::acos(clump_to_acos(tproj.dot(t.col(3))));

  // about 30 degrees tolerance
  if (angular_error > 0.52) {
    OPENMVG_LOG_INFO << "Internal 3rd view reprojection angle check FAIL" << std::endl;
    return false;
  } else {
    OPENMVG_LOG_INFO << "Internal 3rd view reprojection angle check pass" << std::endl;
  }
#endif
  return true;
}


} // namespace trifocal
} // namespace OpenMVG
