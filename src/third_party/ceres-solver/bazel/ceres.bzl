# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2023 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Support for building Ceres Solver with a specific configuration.

CERES_SRCS = ["internal/ceres/" + filename for filename in [
    "accelerate_sparse.cc",
    "array_utils.cc",
    "block_evaluate_preparer.cc",
    "block_jacobi_preconditioner.cc",
    "block_jacobian_writer.cc",
    "block_random_access_dense_matrix.cc",
    "block_random_access_diagonal_matrix.cc",
    "block_random_access_matrix.cc",
    "block_random_access_sparse_matrix.cc",
    "block_sparse_matrix.cc",
    "block_structure.cc",
    "c_api.cc",
    "callbacks.cc",
    "canonical_views_clustering.cc",
    "cgnr_solver.cc",
    "compressed_col_sparse_matrix_utils.cc",
    "compressed_row_jacobian_writer.cc",
    "compressed_row_sparse_matrix.cc",
    "conditioned_cost_function.cc",
    "context.cc",
    "context_impl.cc",
    "coordinate_descent_minimizer.cc",
    "corrector.cc",
    "cost_function.cc",
    "covariance.cc",
    "covariance_impl.cc",
    "dense_cholesky.cc",
    "dense_normal_cholesky_solver.cc",
    "dense_qr.cc",
    "dense_qr_solver.cc",
    "dense_sparse_matrix.cc",
    "detect_structure.cc",
    "dogleg_strategy.cc",
    "dynamic_compressed_row_jacobian_writer.cc",
    "dynamic_compressed_row_sparse_matrix.cc",
    "dynamic_sparse_normal_cholesky_solver.cc",
    "eigensparse.cc",
    "evaluation_callback.cc",
    "evaluator.cc",
    "file.cc",
    "first_order_function.cc",
    "float_suitesparse.cc",
    "function_sample.cc",
    "gradient_checker.cc",
    "gradient_checking_cost_function.cc",
    "gradient_problem.cc",
    "gradient_problem_solver.cc",
    "implicit_schur_complement.cc",
    "inner_product_computer.cc",
    "is_close.cc",
    "iteration_callback.cc",
    "iterative_refiner.cc",
    "iterative_schur_complement_solver.cc",
    "levenberg_marquardt_strategy.cc",
    "line_search.cc",
    "line_search_direction.cc",
    "line_search_minimizer.cc",
    "line_search_preprocessor.cc",
    "linear_least_squares_problems.cc",
    "linear_operator.cc",
    "linear_solver.cc",
    "loss_function.cc",
    "low_rank_inverse_hessian.cc",
    "manifold.cc",
    "minimizer.cc",
    "normal_prior.cc",
    "parallel_for.cc",
    "parallel_invoke.cc",
    "parallel_utils.cc",
    "parameter_block_ordering.cc",
    "partitioned_matrix_view.cc",
    "polynomial.cc",
    "power_series_expansion_preconditioner.cc",
    "preconditioner.cc",
    "preprocessor.cc",
    "problem.cc",
    "problem_impl.cc",
    "program.cc",
    "reorder_program.cc",
    "residual_block.cc",
    "residual_block_utils.cc",
    "schur_complement_solver.cc",
    "schur_eliminator.cc",
    "schur_jacobi_preconditioner.cc",
    "schur_templates.cc",
    "scratch_evaluate_preparer.cc",
    "single_linkage_clustering.cc",
    "solver.cc",
    "solver_utils.cc",
    "sparse_cholesky.cc",
    "sparse_matrix.cc",
    "sparse_normal_cholesky_solver.cc",
    "stringprintf.cc",
    "subset_preconditioner.cc",
    "suitesparse.cc",
    "thread_pool.cc",
    "thread_token_provider.cc",
    "triplet_sparse_matrix.cc",
    "trust_region_minimizer.cc",
    "trust_region_preprocessor.cc",
    "trust_region_step_evaluator.cc",
    "trust_region_strategy.cc",
    "types.cc",
    "visibility_based_preconditioner.cc",
    "visibility.cc",
    "wall_time.cc",
]]

# TODO(rodrigoq): add support to configure Ceres into various permutations,
# like SuiteSparse or not, threading or not, glog or not, and so on.
# See https://github.com/ceres-solver/ceres-solver/issues/335.
def ceres_library(name,
                  restrict_schur_specializations=False):
    # The path to internal/ depends on whether Ceres is the main workspace or
    # an external repository.
    if native.repository_name() != '@':
        internal = 'external/%s/internal' % native.repository_name().lstrip('@')
    else:
        internal = 'internal'

    # The fixed-size Schur eliminator template instantiations incur a large
    # binary size penalty, and are slow to compile, so support disabling them.
    schur_eliminator_copts = []
    if restrict_schur_specializations:
        schur_eliminator_copts.append("-DCERES_RESTRICT_SCHUR_SPECIALIZATION")
        schur_sources = [
            "internal/ceres/generated/schur_eliminator_d_d_d.cc",
            "internal/ceres/generated/partitioned_matrix_view_d_d_d.cc",
        ]
    else:
        schur_sources = native.glob(["internal/ceres/generated/*.cc"])

    native.cc_library(
        name = name,

        # Internal sources, options, and dependencies.
        srcs = CERES_SRCS + schur_sources + native.glob([
            "include/ceres/internal/*.h",
        ]) + native.glob([
            "internal/ceres/*.h",
        ]),

        # These headers are made available to other targets.
        hdrs =
            native.glob(["include/ceres/*.h"]) + native.glob([
                "include/ceres/internal/*.h",
            ]) +

            # This is an empty config and export, since the
            # Bazel-based build does not generate a
            # config.h/export.h. This is fine, since Bazel properly
            # handles propagating -D defines to dependent targets.
            native.glob([
                "config/ceres/internal/config.h",
                "config/ceres/internal/export.h",
            ]),
        copts = [
            "-I" + internal,
            "-Wunused-parameter",
            "-Wno-sign-compare",
        ] + schur_eliminator_copts,

        # These include directories and defines are propagated to other targets
        # depending on Ceres.
        # TODO(keir): These defines are placeholders for now to facilitate getting
        # started with a Bazel build. However, these should become configurable as
        # part of a Skylark Ceres target macro.
        # https://github.com/ceres-solver/ceres-solver/issues/396
        defines = [
            "CERES_EXPORT=",
            "CERES_NO_ACCELERATE_SPARSE",
            "CERES_NO_CHOLMOD_PARTITION",
            "CERES_NO_CUDA",
            "CERES_NO_EIGEN_METIS",
            "CERES_NO_EXPORT=",
            "CERES_NO_LAPACK",
            "CERES_NO_SUITESPARSE",
            "CERES_USE_EIGEN_SPARSE",
        ],
        includes = [
            "config",
            "include",
        ],
        visibility = ["//visibility:public"],
        deps = [
            "@com_gitlab_libeigen_eigen//:eigen",
            "@com_github_google_glog//:glog",
        ],
    )
