// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2016 Google Inc. All rights reserved.
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

#ifndef EXAMPLES_CERES_TYPES_H_
#define EXAMPLES_CERES_TYPES_H_

#include <istream>
#include <map>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace ceres {
namespace examples {

struct Pose3d {
  Eigen::Vector3d p;
  Eigen::Quaterniond q;

  // The name of the data type in the g2o file format.
  static std::string name() {
    return "VERTEX_SE3:QUAT";
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

std::istream& operator>>(std::istream& input, Pose3d& pose) {
  input >> pose.p.x() >> pose.p.y() >> pose.p.z() >> pose.q.x() >>
      pose.q.y() >> pose.q.z() >> pose.q.w();
  // Normalize the quaternion to account for precision loss due to
  // serialization.
  pose.q.normalize();
  return input;
}

typedef std::map<int, Pose3d, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, Pose3d> > >
    MapOfPoses;

// The constraint between two vertices in the pose graph. The constraint is the
// transformation from vertex id_begin to vertex id_end.
struct Constraint3d {
  int id_begin;
  int id_end;

  // The transformation that represents the pose of the end frame E w.r.t. the
  // begin frame B. In other words, it transforms a vector in the E frame to
  // the B frame.
  Pose3d t_be;

  // The inverse of the covariance matrix for the measurement. The order of the
  // entries are x, y, z, delta orientation.
  Eigen::Matrix<double, 6, 6> information;

  // The name of the data type in the g2o file format.
  static std::string name() {
    return "EDGE_SE3:QUAT";
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

std::istream& operator>>(std::istream& input, Constraint3d& constraint) {
  Pose3d& t_be = constraint.t_be;
  input >> constraint.id_begin >> constraint.id_end >> t_be;

  for (int i = 0; i < 6 && input.good(); ++i) {
    for (int j = i; j < 6 && input.good(); ++j) {
      input >> constraint.information(i, j);
      if (i != j) {
        constraint.information(j, i) = constraint.information(i, j);
      }
    }
  }
  return input;
}

typedef std::vector<Constraint3d, Eigen::aligned_allocator<Constraint3d> >
    VectorOfConstraints;

}  // namespace examples
}  // namespace ceres

#endif  // EXAMPLES_CERES_TYPES_H_
