// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 nomoko AG, Srivathsan Murali<srivathsan@nomoko.camera>

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef _NOMOKO_TYPES_H_
#define _NOMOKO_TYPES_H_

namespace nomoko {
  /* brief:
     Contains the camera parameters
     currently contains only a intrinsic matrix for a pinhole camera
     */
  struct Camera {
    Eigen::Matrix3f K;
    unsigned int width;
    unsigned int height;
  }; // struct Camera

  /* brief:
     Contains the information for each view;
     rot       -> Rotation matrix for the pose
     trans     -> Translation vector for the pose
     cameraID  -> index for the camera used of the vector of Cameras
     filepath  -> filepath to the image
     */
  struct View {
    Eigen::Matrix3f rot;
    Eigen::Vector3f trans;
    size_t cameraId;
    std::string filename;
    std::vector<size_t> viewPoints;
  }; // struct View

  /* brief:
     Contains the information for each point in the sparse point cloud
     pos       -> 3D position of the point
     color     -> RGB values for the point
     viewList  -> List of views that see this point
     */
  struct Point {
    Eigen::Vector3f pos;
    std::vector<size_t> viewList;
  };
} // namespace nomoko
#endif // _NOMOKO_TYPES_H_
