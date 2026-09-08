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
// Author: sameeragarwal@google.com (Sameer Agarwal)
//
// This example illustrates the use of the EvaluationCallback, which can be used
// to perform high performance computation of the residual and Jacobians outside
// Ceres (in this case using Eigen's vectorized code) and then the CostFunctions
// just copy these computed residuals and Jacobians appropriately and pass them
// to Ceres Solver.
//
// The results of running this example should be identical to the results
// obtained by running curve_fitting.cc. The only difference between the two
// examples is how the residuals and Jacobians are computed.
//
// The observant reader will note that both here and curve_fitting.cc instead of
// creating one ResidualBlock for each observation one can just do one
// ResidualBlock/CostFunction for the entire problem. The reason for keeping one
// residual per observation is that it is what is needed if and when we need to
// introduce a loss function which is what we do in robust_curve_fitting.cc

#include <iostream>

#include "Eigen/Core"
#include "ceres/ceres.h"
#include "glog/logging.h"

// Data generated using the following octave code.
//   randn('seed', 23497);
//   m = 0.3;
//   c = 0.1;
//   x=[0:0.075:5];
//   y = exp(m * x + c);
//   noise = randn(size(x)) * 0.2;
//   y_observed = y + noise;
//   data = [x', y_observed'];

const int kNumObservations = 67;
// clang-format off
const double data[] = {
  0.000000e+00, 1.133898e+00,
  7.500000e-02, 1.334902e+00,
  1.500000e-01, 1.213546e+00,
  2.250000e-01, 1.252016e+00,
  3.000000e-01, 1.392265e+00,
  3.750000e-01, 1.314458e+00,
  4.500000e-01, 1.472541e+00,
  5.250000e-01, 1.536218e+00,
  6.000000e-01, 1.355679e+00,
  6.750000e-01, 1.463566e+00,
  7.500000e-01, 1.490201e+00,
  8.250000e-01, 1.658699e+00,
  9.000000e-01, 1.067574e+00,
  9.750000e-01, 1.464629e+00,
  1.050000e+00, 1.402653e+00,
  1.125000e+00, 1.713141e+00,
  1.200000e+00, 1.527021e+00,
  1.275000e+00, 1.702632e+00,
  1.350000e+00, 1.423899e+00,
  1.425000e+00, 1.543078e+00,
  1.500000e+00, 1.664015e+00,
  1.575000e+00, 1.732484e+00,
  1.650000e+00, 1.543296e+00,
  1.725000e+00, 1.959523e+00,
  1.800000e+00, 1.685132e+00,
  1.875000e+00, 1.951791e+00,
  1.950000e+00, 2.095346e+00,
  2.025000e+00, 2.361460e+00,
  2.100000e+00, 2.169119e+00,
  2.175000e+00, 2.061745e+00,
  2.250000e+00, 2.178641e+00,
  2.325000e+00, 2.104346e+00,
  2.400000e+00, 2.584470e+00,
  2.475000e+00, 1.914158e+00,
  2.550000e+00, 2.368375e+00,
  2.625000e+00, 2.686125e+00,
  2.700000e+00, 2.712395e+00,
  2.775000e+00, 2.499511e+00,
  2.850000e+00, 2.558897e+00,
  2.925000e+00, 2.309154e+00,
  3.000000e+00, 2.869503e+00,
  3.075000e+00, 3.116645e+00,
  3.150000e+00, 3.094907e+00,
  3.225000e+00, 2.471759e+00,
  3.300000e+00, 3.017131e+00,
  3.375000e+00, 3.232381e+00,
  3.450000e+00, 2.944596e+00,
  3.525000e+00, 3.385343e+00,
  3.600000e+00, 3.199826e+00,
  3.675000e+00, 3.423039e+00,
  3.750000e+00, 3.621552e+00,
  3.825000e+00, 3.559255e+00,
  3.900000e+00, 3.530713e+00,
  3.975000e+00, 3.561766e+00,
  4.050000e+00, 3.544574e+00,
  4.125000e+00, 3.867945e+00,
  4.200000e+00, 4.049776e+00,
  4.275000e+00, 3.885601e+00,
  4.350000e+00, 4.110505e+00,
  4.425000e+00, 4.345320e+00,
  4.500000e+00, 4.161241e+00,
  4.575000e+00, 4.363407e+00,
  4.650000e+00, 4.161576e+00,
  4.725000e+00, 4.619728e+00,
  4.800000e+00, 4.737410e+00,
  4.875000e+00, 4.727863e+00,
  4.950000e+00, 4.669206e+00,
};
// clang-format on

// This implementation of the EvaluationCallback interface also stores the
// residuals and Jacobians that the CostFunction copies their values from.
class MyEvaluationCallback : public ceres::EvaluationCallback {
 public:
  // m and c are passed by reference so that we have access to their values as
  // they evolve over time through the course of optimization.
  MyEvaluationCallback(const double& m, const double& c) : m_(m), c_(c) {
    x_ = Eigen::VectorXd::Zero(kNumObservations);
    y_ = Eigen::VectorXd::Zero(kNumObservations);
    residuals_ = Eigen::VectorXd::Zero(kNumObservations);
    jacobians_ = Eigen::MatrixXd::Zero(kNumObservations, 2);
    for (int i = 0; i < kNumObservations; ++i) {
      x_[i] = data[2 * i];
      y_[i] = data[2 * i + 1];
    }
    PrepareForEvaluation(true, true);
  }

  void PrepareForEvaluation(bool evaluate_jacobians,
                            bool new_evaluation_point) final {
    if (new_evaluation_point) {
      ComputeResidualAndJacobian(evaluate_jacobians);
      jacobians_are_stale_ = !evaluate_jacobians;
    } else {
      if (evaluate_jacobians && jacobians_are_stale_) {
        ComputeResidualAndJacobian(evaluate_jacobians);
        jacobians_are_stale_ = false;
      }
    }
  }

  const Eigen::VectorXd& residuals() const { return residuals_; }
  const Eigen::MatrixXd& jacobians() const { return jacobians_; }
  bool jacobians_are_stale() const { return jacobians_are_stale_; }

 private:
  void ComputeResidualAndJacobian(bool evaluate_jacobians) {
    residuals_ = -(m_ * x_.array() + c_).exp();
    if (evaluate_jacobians) {
      jacobians_.col(0) = residuals_.array() * x_.array();
      jacobians_.col(1) = residuals_;
    }
    residuals_ += y_;
  }

  const double& m_;
  const double& c_;
  Eigen::VectorXd x_;
  Eigen::VectorXd y_;
  Eigen::VectorXd residuals_;
  Eigen::MatrixXd jacobians_;

  // jacobians_are_stale_ keeps track of whether the jacobian matrix matches the
  // residuals or not, we only compute it if we know that Solver is going to
  // need access to it.
  bool jacobians_are_stale_ = true;
};

// As the name implies this CostFunction does not do any computation, it just
// copies the appropriate residual and Jacobian from the matrices stored in
// MyEvaluationCallback.
class CostAndJacobianCopyingCostFunction
    : public ceres::SizedCostFunction<1, 1, 1> {
 public:
  CostAndJacobianCopyingCostFunction(
      int index, const MyEvaluationCallback& evaluation_callback)
      : index_(index), evaluation_callback_(evaluation_callback) {}
  ~CostAndJacobianCopyingCostFunction() override = default;

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    residuals[0] = evaluation_callback_.residuals()(index_);
    if (!jacobians) return true;

    // Ensure that we are not using stale Jacobians.
    CHECK(!evaluation_callback_.jacobians_are_stale());

    if (jacobians[0] != nullptr)
      jacobians[0][0] = evaluation_callback_.jacobians()(index_, 0);
    if (jacobians[1] != nullptr)
      jacobians[1][0] = evaluation_callback_.jacobians()(index_, 1);
    return true;
  }

 private:
  int index_ = -1;
  const MyEvaluationCallback& evaluation_callback_;
};

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  const double initial_m = 0.0;
  const double initial_c = 0.0;
  double m = initial_m;
  double c = initial_c;

  MyEvaluationCallback evaluation_callback(m, c);
  ceres::Problem::Options problem_options;
  problem_options.evaluation_callback = &evaluation_callback;
  ceres::Problem problem(problem_options);
  for (int i = 0; i < kNumObservations; ++i) {
    problem.AddResidualBlock(
        new CostAndJacobianCopyingCostFunction(i, evaluation_callback),
        nullptr,
        &m,
        &c);
  }

  ceres::Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";
  std::cout << "Initial m: " << initial_m << " c: " << initial_c << "\n";
  std::cout << "Final   m: " << m << " c: " << c << "\n";
  return 0;
}
