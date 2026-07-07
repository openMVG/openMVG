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
// Author: richie.stebbing@gmail.com (Richard Stebbing)
//         sameeragarwal@google.com (Sameer Agarwal)
//
// Based on examples/ellipse_approximation.cc

#include <cmath>
#include <utility>
#include <vector>

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// Data generated with the following Python code.
//   import numpy as np
//   np.random.seed(1337)
//   t = np.linspace(0.0, 2.0 * np.pi, 212, endpoint=False)
//   t += 2.0 * np.pi * 0.01 * np.random.randn(t.size)
//   theta = np.deg2rad(15)
//   a, b = np.cos(theta), np.sin(theta)
//   R = np.array([[a, -b],
//                 [b, a]])
//   Y = np.dot(np.c_[4.0 * np.cos(t), np.sin(t)], R.T)

const int kYRows = 212;
const int kYCols = 2;
// clang-format off
const double kYData[kYRows * kYCols] = {
  +3.871364e+00, +9.916027e-01,
  +3.864003e+00, +1.034148e+00,
  +3.850651e+00, +1.072202e+00,
  +3.868350e+00, +1.014408e+00,
  +3.796381e+00, +1.153021e+00,
  +3.857138e+00, +1.056102e+00,
  +3.787532e+00, +1.162215e+00,
  +3.704477e+00, +1.227272e+00,
  +3.564711e+00, +1.294959e+00,
  +3.754363e+00, +1.191948e+00,
  +3.482098e+00, +1.322725e+00,
  +3.602777e+00, +1.279658e+00,
  +3.585433e+00, +1.286858e+00,
  +3.347505e+00, +1.356415e+00,
  +3.220855e+00, +1.378914e+00,
  +3.558808e+00, +1.297174e+00,
  +3.403618e+00, +1.343809e+00,
  +3.179828e+00, +1.384721e+00,
  +3.054789e+00, +1.398759e+00,
  +3.294153e+00, +1.366808e+00,
  +3.247312e+00, +1.374813e+00,
  +2.988547e+00, +1.404247e+00,
  +3.114508e+00, +1.392698e+00,
  +2.899226e+00, +1.409802e+00,
  +2.533256e+00, +1.414778e+00,
  +2.654773e+00, +1.415909e+00,
  +2.565100e+00, +1.415313e+00,
  +2.976456e+00, +1.405118e+00,
  +2.484200e+00, +1.413640e+00,
  +2.324751e+00, +1.407476e+00,
  +1.930468e+00, +1.378221e+00,
  +2.329017e+00, +1.407688e+00,
  +1.760640e+00, +1.360319e+00,
  +2.147375e+00, +1.396603e+00,
  +1.741989e+00, +1.358178e+00,
  +1.743859e+00, +1.358394e+00,
  +1.557372e+00, +1.335208e+00,
  +1.280551e+00, +1.295087e+00,
  +1.429880e+00, +1.317546e+00,
  +1.213485e+00, +1.284400e+00,
  +9.168172e-01, +1.232870e+00,
  +1.311141e+00, +1.299839e+00,
  +1.231969e+00, +1.287382e+00,
  +7.453773e-01, +1.200049e+00,
  +6.151587e-01, +1.173683e+00,
  +5.935666e-01, +1.169193e+00,
  +2.538707e-01, +1.094227e+00,
  +6.806136e-01, +1.187089e+00,
  +2.805447e-01, +1.100405e+00,
  +6.184807e-01, +1.174371e+00,
  +1.170550e-01, +1.061762e+00,
  +2.890507e-01, +1.102365e+00,
  +3.834234e-01, +1.123772e+00,
  +3.980161e-04, +1.033061e+00,
  -3.651680e-01, +9.370367e-01,
  -8.386351e-01, +7.987201e-01,
  -8.105704e-01, +8.073702e-01,
  -8.735139e-01, +7.878886e-01,
  -9.913836e-01, +7.506100e-01,
  -8.784011e-01, +7.863636e-01,
  -1.181440e+00, +6.882566e-01,
  -1.229556e+00, +6.720191e-01,
  -1.035839e+00, +7.362765e-01,
  -8.031520e-01, +8.096470e-01,
  -1.539136e+00, +5.629549e-01,
  -1.755423e+00, +4.817306e-01,
  -1.337589e+00, +6.348763e-01,
  -1.836966e+00, +4.499485e-01,
  -1.913367e+00, +4.195617e-01,
  -2.126467e+00, +3.314900e-01,
  -1.927625e+00, +4.138238e-01,
  -2.339862e+00, +2.379074e-01,
  -1.881736e+00, +4.322152e-01,
  -2.116753e+00, +3.356163e-01,
  -2.255733e+00, +2.754930e-01,
  -2.555834e+00, +1.368473e-01,
  -2.770277e+00, +2.895711e-02,
  -2.563376e+00, +1.331890e-01,
  -2.826715e+00, -9.000818e-04,
  -2.978191e+00, -8.457804e-02,
  -3.115855e+00, -1.658786e-01,
  -2.982049e+00, -8.678322e-02,
  -3.307892e+00, -2.902083e-01,
  -3.038346e+00, -1.194222e-01,
  -3.190057e+00, -2.122060e-01,
  -3.279086e+00, -2.705777e-01,
  -3.322028e+00, -2.999889e-01,
  -3.122576e+00, -1.699965e-01,
  -3.551973e+00, -4.768674e-01,
  -3.581866e+00, -5.032175e-01,
  -3.497799e+00, -4.315203e-01,
  -3.565384e+00, -4.885602e-01,
  -3.699493e+00, -6.199815e-01,
  -3.585166e+00, -5.061925e-01,
  -3.758914e+00, -6.918275e-01,
  -3.741104e+00, -6.689131e-01,
  -3.688331e+00, -6.077239e-01,
  -3.810425e+00, -7.689015e-01,
  -3.791829e+00, -7.386911e-01,
  -3.789951e+00, -7.358189e-01,
  -3.823100e+00, -7.918398e-01,
  -3.857021e+00, -8.727074e-01,
  -3.858250e+00, -8.767645e-01,
  -3.872100e+00, -9.563174e-01,
  -3.864397e+00, -1.032630e+00,
  -3.846230e+00, -1.081669e+00,
  -3.834799e+00, -1.102536e+00,
  -3.866684e+00, -1.022901e+00,
  -3.808643e+00, -1.139084e+00,
  -3.868840e+00, -1.011569e+00,
  -3.791071e+00, -1.158615e+00,
  -3.797999e+00, -1.151267e+00,
  -3.696278e+00, -1.232314e+00,
  -3.779007e+00, -1.170504e+00,
  -3.622855e+00, -1.270793e+00,
  -3.647249e+00, -1.259166e+00,
  -3.655412e+00, -1.255042e+00,
  -3.573218e+00, -1.291696e+00,
  -3.638019e+00, -1.263684e+00,
  -3.498409e+00, -1.317750e+00,
  -3.304143e+00, -1.364970e+00,
  -3.183001e+00, -1.384295e+00,
  -3.202456e+00, -1.381599e+00,
  -3.244063e+00, -1.375332e+00,
  -3.233308e+00, -1.377019e+00,
  -3.060112e+00, -1.398264e+00,
  -3.078187e+00, -1.396517e+00,
  -2.689594e+00, -1.415761e+00,
  -2.947662e+00, -1.407039e+00,
  -2.854490e+00, -1.411860e+00,
  -2.660499e+00, -1.415900e+00,
  -2.875955e+00, -1.410930e+00,
  -2.675385e+00, -1.415848e+00,
  -2.813155e+00, -1.413363e+00,
  -2.417673e+00, -1.411512e+00,
  -2.725461e+00, -1.415373e+00,
  -2.148334e+00, -1.396672e+00,
  -2.108972e+00, -1.393738e+00,
  -2.029905e+00, -1.387302e+00,
  -2.046214e+00, -1.388687e+00,
  -2.057402e+00, -1.389621e+00,
  -1.650250e+00, -1.347160e+00,
  -1.806764e+00, -1.365469e+00,
  -1.206973e+00, -1.283343e+00,
  -8.029259e-01, -1.211308e+00,
  -1.229551e+00, -1.286993e+00,
  -1.101507e+00, -1.265754e+00,
  -9.110645e-01, -1.231804e+00,
  -1.110046e+00, -1.267211e+00,
  -8.465274e-01, -1.219677e+00,
  -7.594163e-01, -1.202818e+00,
  -8.023823e-01, -1.211203e+00,
  -3.732519e-01, -1.121494e+00,
  -1.918373e-01, -1.079668e+00,
  -4.671988e-01, -1.142253e+00,
  -4.033645e-01, -1.128215e+00,
  -1.920740e-01, -1.079724e+00,
  -3.022157e-01, -1.105389e+00,
  -1.652831e-01, -1.073354e+00,
  +4.671625e-01, -9.085886e-01,
  +5.940178e-01, -8.721832e-01,
  +3.147557e-01, -9.508290e-01,
  +6.383631e-01, -8.591867e-01,
  +9.888923e-01, -7.514088e-01,
  +7.076339e-01, -8.386023e-01,
  +1.326682e+00, -6.386698e-01,
  +1.149834e+00, -6.988221e-01,
  +1.257742e+00, -6.624207e-01,
  +1.492352e+00, -5.799632e-01,
  +1.595574e+00, -5.421766e-01,
  +1.240173e+00, -6.684113e-01,
  +1.706612e+00, -5.004442e-01,
  +1.873984e+00, -4.353002e-01,
  +1.985633e+00, -3.902561e-01,
  +1.722880e+00, -4.942329e-01,
  +2.095182e+00, -3.447402e-01,
  +2.018118e+00, -3.768991e-01,
  +2.422702e+00, -1.999563e-01,
  +2.370611e+00, -2.239326e-01,
  +2.152154e+00, -3.205250e-01,
  +2.525121e+00, -1.516499e-01,
  +2.422116e+00, -2.002280e-01,
  +2.842806e+00, +9.536372e-03,
  +3.030128e+00, +1.146027e-01,
  +2.888424e+00, +3.433444e-02,
  +2.991609e+00, +9.226409e-02,
  +2.924807e+00, +5.445844e-02,
  +3.007772e+00, +1.015875e-01,
  +2.781973e+00, -2.282382e-02,
  +3.164737e+00, +1.961781e-01,
  +3.237671e+00, +2.430139e-01,
  +3.046123e+00, +1.240014e-01,
  +3.414834e+00, +3.669060e-01,
  +3.436591e+00, +3.833600e-01,
  +3.626207e+00, +5.444311e-01,
  +3.223325e+00, +2.336361e-01,
  +3.511963e+00, +4.431060e-01,
  +3.698380e+00, +6.187442e-01,
  +3.670244e+00, +5.884943e-01,
  +3.558833e+00, +4.828230e-01,
  +3.661807e+00, +5.797689e-01,
  +3.767261e+00, +7.030893e-01,
  +3.801065e+00, +7.532650e-01,
  +3.828523e+00, +8.024454e-01,
  +3.840719e+00, +8.287032e-01,
  +3.848748e+00, +8.485921e-01,
  +3.865801e+00, +9.066551e-01,
  +3.870983e+00, +9.404873e-01,
  +3.870263e+00, +1.001884e+00,
  +3.864462e+00, +1.032374e+00,
  +3.870542e+00, +9.996121e-01,
  +3.865424e+00, +1.028474e+00
};
// clang-format on

ConstMatrixRef kY(kYData, kYRows, kYCols);

class PointToLineSegmentContourCostFunction : public CostFunction {
 public:
  // This class needs to have an Eigen aligned operator new as it contains
  // fixed-size Eigen types.
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  PointToLineSegmentContourCostFunction(const int num_segments,
                                        Eigen::Vector2d y)
      : num_segments_(num_segments), y_(std::move(y)) {
    // The first parameter is the preimage position.
    mutable_parameter_block_sizes()->push_back(1);
    // The next parameters are the control points for the line segment contour.
    for (int i = 0; i < num_segments_; ++i) {
      mutable_parameter_block_sizes()->push_back(2);
    }
    set_num_residuals(2);
  }

  bool Evaluate(const double* const* x,
                double* residuals,
                double** jacobians) const final {
    // Convert the preimage position `t` into a segment index `i0` and the
    // line segment interpolation parameter `u`. `i1` is the index of the next
    // control point.
    const double t = ModuloNumSegments(*x[0]);
    CHECK_GE(t, 0.0);
    CHECK_LT(t, num_segments_);
    const int i0 = floor(t), i1 = (i0 + 1) % num_segments_;
    const double u = t - i0;

    // Linearly interpolate between control points `i0` and `i1`.
    residuals[0] = y_[0] - ((1.0 - u) * x[1 + i0][0] + u * x[1 + i1][0]);
    residuals[1] = y_[1] - ((1.0 - u) * x[1 + i0][1] + u * x[1 + i1][1]);

    if (jacobians == nullptr) {
      return true;
    }

    if (jacobians[0] != nullptr) {
      jacobians[0][0] = x[1 + i0][0] - x[1 + i1][0];
      jacobians[0][1] = x[1 + i0][1] - x[1 + i1][1];
    }
    for (int i = 0; i < num_segments_; ++i) {
      if (jacobians[i + 1] != nullptr) {
        MatrixRef(jacobians[i + 1], 2, 2).setZero();
        if (i == i0) {
          jacobians[i + 1][0] = -(1.0 - u);
          jacobians[i + 1][3] = -(1.0 - u);
        } else if (i == i1) {
          jacobians[i + 1][0] = -u;
          jacobians[i + 1][3] = -u;
        }
      }
    }
    return true;
  }

  static CostFunction* Create(const int num_segments,
                              const Eigen::Vector2d& y) {
    return new PointToLineSegmentContourCostFunction(num_segments, y);
  }

 private:
  inline double ModuloNumSegments(const double t) const {
    return t - num_segments_ * floor(t / num_segments_);
  }

  const int num_segments_;
  const Eigen::Vector2d y_;
};

class EuclideanDistanceFunctor {
 public:
  explicit EuclideanDistanceFunctor(const double sqrt_weight)
      : sqrt_weight_(sqrt_weight) {}

  template <typename T>
  bool operator()(const T* x0, const T* x1, T* residuals) const {
    residuals[0] = sqrt_weight_ * (x0[0] - x1[0]);
    residuals[1] = sqrt_weight_ * (x0[1] - x1[1]);
    return true;
  }

  static CostFunction* Create(const double sqrt_weight) {
    return new AutoDiffCostFunction<EuclideanDistanceFunctor, 2, 2, 2>(
        new EuclideanDistanceFunctor(sqrt_weight));
  }

 private:
  const double sqrt_weight_;
};

TEST(DynamicSparsity, StaticAndDynamicSparsityProduceSameSolution) {
  // Skip test if there is no sparse linear algebra library that
  // supports dynamic sparsity.
  if (!IsSparseLinearAlgebraLibraryTypeAvailable(SUITE_SPARSE) &&
      !IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    return;
  }

  // Problem configuration.
  const int num_segments = 151;
  const double regularization_weight = 1e-2;

  // `X` is the matrix of control points which make up the contour of line
  // segments. The number of control points is equal to the number of line
  // segments because the contour is closed.
  //
  // Initialize `X` to points on the unit circle.
  Vector w(num_segments + 1);
  w.setLinSpaced(num_segments + 1, 0.0, 2.0 * constants::pi);
  w.conservativeResize(num_segments);
  Matrix X(num_segments, 2);
  X.col(0) = w.array().cos();
  X.col(1) = w.array().sin();

  // Each data point has an associated preimage position on the line segment
  // contour. For each data point we initialize the preimage positions to
  // the index of the closest control point.
  const int num_observations = kY.rows();
  Vector t(num_observations);
  for (int i = 0; i < num_observations; ++i) {
    (X.rowwise() - kY.row(i)).rowwise().squaredNorm().minCoeff(&t[i]);
  }

  Problem problem;

  // For each data point add a residual which measures its distance to its
  // corresponding position on the line segment contour.
  std::vector<double*> parameter_blocks(1 + num_segments);
  parameter_blocks[0] = nullptr;
  for (int i = 0; i < num_segments; ++i) {
    parameter_blocks[i + 1] = X.data() + 2 * i;
  }
  for (int i = 0; i < num_observations; ++i) {
    parameter_blocks[0] = &t[i];
    problem.AddResidualBlock(
        PointToLineSegmentContourCostFunction::Create(num_segments, kY.row(i)),
        nullptr,
        parameter_blocks);
  }

  // Add regularization to minimize the length of the line segment contour.
  for (int i = 0; i < num_segments; ++i) {
    problem.AddResidualBlock(
        EuclideanDistanceFunctor::Create(sqrt(regularization_weight)),
        nullptr,
        X.data() + 2 * i,
        X.data() + 2 * ((i + 1) % num_segments));
  }

  Solver::Options options;
  options.max_num_iterations = 100;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  // Only SuiteSparse & EigenSparse currently support dynamic sparsity.
  options.sparse_linear_algebra_library_type =
#if !defined(CERES_NO_SUITESPARSE)
      ceres::SUITE_SPARSE;
#elif defined(CERES_USE_EIGEN_SPARSE)
      ceres::EIGEN_SPARSE;
#endif

  // First, solve `X` and `t` jointly with dynamic_sparsity = true.
  Matrix X0 = X;
  Vector t0 = t;
  options.dynamic_sparsity = false;
  Solver::Summary static_summary;
  Solve(options, &problem, &static_summary);
  EXPECT_EQ(static_summary.termination_type, CONVERGENCE)
      << static_summary.FullReport();

  X = X0;
  t = t0;
  options.dynamic_sparsity = true;
  Solver::Summary dynamic_summary;
  Solve(options, &problem, &dynamic_summary);
  EXPECT_EQ(dynamic_summary.termination_type, CONVERGENCE)
      << dynamic_summary.FullReport();

  EXPECT_NEAR(static_summary.final_cost,
              dynamic_summary.final_cost,
              std::numeric_limits<double>::epsilon())
      << "Static: \n"
      << static_summary.FullReport() << "\nDynamic: \n"
      << dynamic_summary.FullReport();
}

}  // namespace ceres::internal
