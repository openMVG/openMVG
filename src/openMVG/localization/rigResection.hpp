/* 
 * File:   rigResection.hpp
 * Author: sgaspari
 *
 * Created on January 2, 2016, 5:51 PM
 */

#pragma once

#include <openMVG/types.hpp>
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>
#include <openMVG/geometry/pose3.hpp>

#include <vector>

namespace openMVG{
namespace localization{

#if HAVE_OPENGV

bool rigResection(const std::vector<openMVG::Mat2X> &pts2d, 
                  const std::vector<openMVG::Mat3X> &pts3d,
                  const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                  const std::vector<geometry::Pose3 > &vec_subPoses,
                  geometry::Pose3 &rigPose,
                  std::vector<std::vector<std::size_t> > &inliers,
                  double threshold = 1e-6,
                  size_t maxIterations = 100, 
                  bool verbosity = true);

#endif

}
}
