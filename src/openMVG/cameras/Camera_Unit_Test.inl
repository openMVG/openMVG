// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <random>

// Macro used to perform unit test on the OpenMVG Camera
//
// - Generate random point inside the image domain
// - For each point:
//   - Check bijection between transformation between camera and image domain
//   - Check bijection of the distortion function
//   - Check bijection of the bearing vector and its projection
#define Test_camera(cam) \
{ \
 \
  std::default_random_engine gen; \
  std::uniform_int_distribution<> \
    rand_x(0, cam.w()-1), \
    rand_y(0, cam.h()-1); \
 \
  static const int nb_sampled_pts = cam.w() * cam.h() * 0.10; \
  static const double epsilon = 1e-4; \
 \
  for (int i = 0; i < nb_sampled_pts; ++i) \
  { \
    /* generate random point inside the image domain */ \
    const Vec2 ptImage = {rand_x(gen), rand_y(gen)}; \
 \
    /* Check bijection between transformation between camera and image domain */ \
    EXPECT_MATRIX_NEAR( ptImage, cam.cam2ima(cam.ima2cam(ptImage)), epsilon); \
 \
    /* Check bijection of the distortion function */ \
    const Vec2 ptCamera = cam.ima2cam(ptImage); \
    EXPECT_MATRIX_NEAR( ptCamera, cam.remove_disto(cam.add_disto(ptCamera)), epsilon); \
    EXPECT_MATRIX_NEAR( ptImage, cam.cam2ima(cam.remove_disto(cam.add_disto(ptCamera))), epsilon); \
 \
    if (cam.have_disto()) \
    { /* Assert that distortion field is not null and it has moved the initial provided point */ \
      EXPECT_FALSE( (cam.add_disto(ptCamera) == cam.remove_disto(cam.add_disto(ptCamera))) ); \
    } \
 \
    /* Check bijection of the bearing vector and its projection */ \
    EXPECT_MATRIX_NEAR( \
      ptImage, \
      /* (ptImage -> bearing -> projected image point -> remove the distortion) */ \
      cam.cam2ima(cam.remove_disto(cam.ima2cam(cam.project(geometry::Pose3(), cam(ptImage))))), \
      epsilon); \
  } \
}

