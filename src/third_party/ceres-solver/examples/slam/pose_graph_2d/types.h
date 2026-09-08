// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: vitus@google.com (Michael Vitus)
//
// Defines the types used in the 2D pose graph SLAM formulation. Each vertex of
// the graph has a unique integer ID with a position and orientation. There are
// delta transformation constraints between two vertices.

#ifndef CERES_EXAMPLES_POSE_GRAPH_2D_TYPES_H_
#define CERES_EXAMPLES_POSE_GRAPH_2D_TYPES_H_

#include <fstream>

#include "Eigen/Core"
#include "normalize_angle.h"

namespace ceres::examples {

// The state for each vertex in the pose graph.
struct Pose2d {
  double x;
  double y;
  double yaw_radians;

  // The name of the data type in the g2o file format.
  static std::string name() { return "VERTEX_SE2"; }
};

inline std::istream& operator>>(std::istream& input, Pose2d& pose) {
  input >> pose.x >> pose.y >> pose.yaw_radians;
  // Normalize the angle between -pi to pi.
  pose.yaw_radians = NormalizeAngle(pose.yaw_radians);
  return input;
}

// The constraint between two vertices in the pose graph. The constraint is the
// transformation from vertex id_begin to vertex id_end.
struct Constraint2d {
  int id_begin;
  int id_end;

  double x;
  double y;
  double yaw_radians;

  // The inverse of the covariance matrix for the measurement. The order of the
  // entries are x, y, and yaw.
  Eigen::Matrix3d information;

  // The name of the data type in the g2o file format.
  static std::string name() { return "EDGE_SE2"; }
};

inline std::istream& operator>>(std::istream& input, Constraint2d& constraint) {
  input >> constraint.id_begin >> constraint.id_end >> constraint.x >>
      constraint.y >> constraint.yaw_radians >> constraint.information(0, 0) >>
      constraint.information(0, 1) >> constraint.information(0, 2) >>
      constraint.information(1, 1) >> constraint.information(1, 2) >>
      constraint.information(2, 2);

  // Set the lower triangular part of the information matrix.
  constraint.information(1, 0) = constraint.information(0, 1);
  constraint.information(2, 0) = constraint.information(0, 2);
  constraint.information(2, 1) = constraint.information(1, 2);

  // Normalize the angle between -pi to pi.
  constraint.yaw_radians = NormalizeAngle(constraint.yaw_radians);
  return input;
}

}  // namespace ceres::examples

#endif  // CERES_EXAMPLES_POSE_GRAPH_2D_TYPES_H_
