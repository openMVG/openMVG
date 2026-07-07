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
// Author: strandmark@google.com (Petter Strandmark)
//
// Denoising using Fields of Experts and the Ceres minimizer.
//
// Note that for good denoising results the weighting between the data term
// and the Fields of Experts term needs to be adjusted. This is discussed
// in [1]. This program assumes Gaussian noise. The noise model can be changed
// by substituting another function for QuadraticCostFunction.
//
// [1] S. Roth and M.J. Black. "Fields of Experts." International Journal of
//     Computer Vision, 82(2):205--229, 2009.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "ceres/ceres.h"
#include "fields_of_experts.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "pgm_image.h"

DEFINE_string(input, "", "File to which the output image should be written");
DEFINE_string(foe_file, "", "FoE file to use");
DEFINE_string(output, "", "File to which the output image should be written");
DEFINE_double(sigma, 20.0, "Standard deviation of noise");
DEFINE_string(trust_region_strategy,
              "levenberg_marquardt",
              "Options are: levenberg_marquardt, dogleg.");
DEFINE_string(dogleg,
              "traditional_dogleg",
              "Options are: traditional_dogleg,"
              "subspace_dogleg.");
DEFINE_string(linear_solver,
              "sparse_normal_cholesky",
              "Options are: "
              "sparse_normal_cholesky and cgnr.");
DEFINE_string(preconditioner,
              "jacobi",
              "Options are: "
              "identity, jacobi, subset");
DEFINE_string(sparse_linear_algebra_library,
              "suite_sparse",
              "Options are: suite_sparse, cx_sparse and eigen_sparse");
DEFINE_double(eta,
              1e-2,
              "Default value for eta. Eta determines the "
              "accuracy of each linear solve of the truncated newton step. "
              "Changing this parameter can affect solve performance.");
DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 10, "Number of iterations.");
DEFINE_bool(nonmonotonic_steps,
            false,
            "Trust region algorithm can use"
            " nonmonotic steps.");
DEFINE_bool(inner_iterations,
            false,
            "Use inner iterations to non-linearly "
            "refine each successful trust region step.");
DEFINE_bool(mixed_precision_solves, false, "Use mixed precision solves.");
DEFINE_int32(max_num_refinement_iterations,
             0,
             "Iterative refinement iterations");
DEFINE_bool(line_search,
            false,
            "Use a line search instead of trust region "
            "algorithm.");
DEFINE_double(subset_fraction,
              0.2,
              "The fraction of residual blocks to use for the"
              " subset preconditioner.");

namespace ceres::examples {
namespace {

// This cost function is used to build the data term.
//
//   f_i(x) = a * (x_i - b)^2
//
class QuadraticCostFunction : public ceres::SizedCostFunction<1, 1> {
 public:
  QuadraticCostFunction(double a, double b) : sqrta_(std::sqrt(a)), b_(b) {}
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    const double x = parameters[0][0];
    residuals[0] = sqrta_ * (x - b_);
    if (jacobians != nullptr && jacobians[0] != nullptr) {
      jacobians[0][0] = sqrta_;
    }
    return true;
  }

 private:
  double sqrta_, b_;
};

// Creates a Fields of Experts MAP inference problem.
void CreateProblem(const FieldsOfExperts& foe,
                   const PGMImage<double>& image,
                   Problem* problem,
                   PGMImage<double>* solution) {
  // Create the data term
  CHECK_GT(CERES_GET_FLAG(FLAGS_sigma), 0.0);
  const double coefficient =
      1 / (2.0 * CERES_GET_FLAG(FLAGS_sigma) * CERES_GET_FLAG(FLAGS_sigma));
  for (int index = 0; index < image.NumPixels(); ++index) {
    ceres::CostFunction* cost_function = new QuadraticCostFunction(
        coefficient, image.PixelFromLinearIndex(index));
    problem->AddResidualBlock(
        cost_function, nullptr, solution->MutablePixelFromLinearIndex(index));
  }

  // Create Ceres cost and loss functions for regularization. One is needed for
  // each filter.
  std::vector<ceres::LossFunction*> loss_function(foe.NumFilters());
  std::vector<ceres::CostFunction*> cost_function(foe.NumFilters());
  for (int alpha_index = 0; alpha_index < foe.NumFilters(); ++alpha_index) {
    loss_function[alpha_index] = foe.NewLossFunction(alpha_index);
    cost_function[alpha_index] = foe.NewCostFunction(alpha_index);
  }

  // Add FoE regularization for each patch in the image.
  for (int x = 0; x < image.width() - (foe.Size() - 1); ++x) {
    for (int y = 0; y < image.height() - (foe.Size() - 1); ++y) {
      // Build a vector with the pixel indices of this patch.
      std::vector<double*> pixels;
      const std::vector<int>& x_delta_indices = foe.GetXDeltaIndices();
      const std::vector<int>& y_delta_indices = foe.GetYDeltaIndices();
      for (int i = 0; i < foe.NumVariables(); ++i) {
        double* pixel = solution->MutablePixel(x + x_delta_indices[i],
                                               y + y_delta_indices[i]);
        pixels.push_back(pixel);
      }
      // For this patch with coordinates (x, y), we will add foe.NumFilters()
      // terms to the objective function.
      for (int alpha_index = 0; alpha_index < foe.NumFilters(); ++alpha_index) {
        problem->AddResidualBlock(
            cost_function[alpha_index], loss_function[alpha_index], pixels);
      }
    }
  }
}

void SetLinearSolver(Solver::Options* options) {
  CHECK(StringToLinearSolverType(CERES_GET_FLAG(FLAGS_linear_solver),
                                 &options->linear_solver_type));
  CHECK(StringToPreconditionerType(CERES_GET_FLAG(FLAGS_preconditioner),
                                   &options->preconditioner_type));
  CHECK(StringToSparseLinearAlgebraLibraryType(
      CERES_GET_FLAG(FLAGS_sparse_linear_algebra_library),
      &options->sparse_linear_algebra_library_type));
  options->use_mixed_precision_solves =
      CERES_GET_FLAG(FLAGS_mixed_precision_solves);
  options->max_num_refinement_iterations =
      CERES_GET_FLAG(FLAGS_max_num_refinement_iterations);
}

void SetMinimizerOptions(Solver::Options* options) {
  options->max_num_iterations = CERES_GET_FLAG(FLAGS_num_iterations);
  options->minimizer_progress_to_stdout = true;
  options->num_threads = CERES_GET_FLAG(FLAGS_num_threads);
  options->eta = CERES_GET_FLAG(FLAGS_eta);
  options->use_nonmonotonic_steps = CERES_GET_FLAG(FLAGS_nonmonotonic_steps);
  if (CERES_GET_FLAG(FLAGS_line_search)) {
    options->minimizer_type = ceres::LINE_SEARCH;
  }

  CHECK(StringToTrustRegionStrategyType(
      CERES_GET_FLAG(FLAGS_trust_region_strategy),
      &options->trust_region_strategy_type));
  CHECK(
      StringToDoglegType(CERES_GET_FLAG(FLAGS_dogleg), &options->dogleg_type));
  options->use_inner_iterations = CERES_GET_FLAG(FLAGS_inner_iterations);
}

// Solves the FoE problem using Ceres and post-processes it to make sure the
// solution stays within [0, 255].
void SolveProblem(Problem* problem, PGMImage<double>* solution) {
  // These parameters may be experimented with. For example, ceres::DOGLEG tends
  // to be faster for 2x2 filters, but gives solutions with slightly higher
  // objective function value.
  ceres::Solver::Options options;
  SetMinimizerOptions(&options);
  SetLinearSolver(&options);
  options.function_tolerance = 1e-3;  // Enough for denoising.

  if (options.linear_solver_type == ceres::CGNR &&
      options.preconditioner_type == ceres::SUBSET) {
    std::vector<ResidualBlockId> residual_blocks;
    problem->GetResidualBlocks(&residual_blocks);

    // To use the SUBSET preconditioner we need to provide a list of
    // residual blocks (rows of the Jacobian). The denoising problem
    // has fairly general sparsity, and there is no apriori reason to
    // select one residual block over another, so we will randomly
    // subsample the residual blocks with probability subset_fraction.
    std::default_random_engine engine;
    std::uniform_real_distribution<> distribution(0, 1);  // rage 0 - 1
    for (auto residual_block : residual_blocks) {
      if (distribution(engine) <= CERES_GET_FLAG(FLAGS_subset_fraction)) {
        options.residual_blocks_for_subset_preconditioner.insert(
            residual_block);
      }
    }
  }

  ceres::Solver::Summary summary;
  ceres::Solve(options, problem, &summary);
  std::cout << summary.FullReport() << "\n";

  // Make the solution stay in [0, 255].
  for (int x = 0; x < solution->width(); ++x) {
    for (int y = 0; y < solution->height(); ++y) {
      *solution->MutablePixel(x, y) =
          std::min(255.0, std::max(0.0, solution->Pixel(x, y)));
    }
  }
}

}  // namespace
}  // namespace ceres::examples

int main(int argc, char** argv) {
  using namespace ceres::examples;
  GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  if (CERES_GET_FLAG(FLAGS_input).empty()) {
    std::cerr << "Please provide an image file name using -input.\n";
    return 1;
  }

  if (CERES_GET_FLAG(FLAGS_foe_file).empty()) {
    std::cerr << "Please provide a Fields of Experts file name using -foe_file."
                 "\n";
    return 1;
  }

  // Load the Fields of Experts filters from file.
  FieldsOfExperts foe;
  if (!foe.LoadFromFile(CERES_GET_FLAG(FLAGS_foe_file))) {
    std::cerr << "Loading \"" << CERES_GET_FLAG(FLAGS_foe_file)
              << "\" failed.\n";
    return 2;
  }

  // Read the images
  PGMImage<double> image(CERES_GET_FLAG(FLAGS_input));
  if (image.width() == 0) {
    std::cerr << "Reading \"" << CERES_GET_FLAG(FLAGS_input) << "\" failed.\n";
    return 3;
  }
  PGMImage<double> solution(image.width(), image.height());
  solution.Set(0.0);

  ceres::Problem problem;
  CreateProblem(foe, image, &problem, &solution);

  SolveProblem(&problem, &solution);

  if (!CERES_GET_FLAG(FLAGS_output).empty()) {
    CHECK(solution.WriteToFile(CERES_GET_FLAG(FLAGS_output)))
        << "Writing \"" << CERES_GET_FLAG(FLAGS_output) << "\" failed.";
  }

  return 0;
}
