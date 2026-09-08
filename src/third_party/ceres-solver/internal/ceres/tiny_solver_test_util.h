
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
// Author: mierle@gmail.com (Keir Mierle)

#ifndef CERES_INTERNAL_TINY_SOLVER_TEST_UTIL_H_
#define CERES_INTERNAL_TINY_SOLVER_TEST_UTIL_H_

namespace ceres {

template <typename T>
bool EvaluateResidualsAndJacobians(const T* parameters,
                                   T* residuals,
                                   T* jacobian) {
  T x = parameters[0];
  T y = parameters[1];
  T z = parameters[2];

  residuals[0] = x + static_cast<T>(2) * y + static_cast<T>(4) * z;
  residuals[1] = y * z;

  if (jacobian) {
    jacobian[0 * 2 + 0] = static_cast<T>(1);
    jacobian[0 * 2 + 1] = static_cast<T>(0);

    jacobian[1 * 2 + 0] = static_cast<T>(2);
    jacobian[1 * 2 + 1] = z;

    jacobian[2 * 2 + 0] = static_cast<T>(4);
    jacobian[2 * 2 + 1] = y;
  }
  return true;
}

}  // namespace ceres

#endif  // CERES_INTERNAL_TINY_SOLVER_TEST_UTIL_H_
