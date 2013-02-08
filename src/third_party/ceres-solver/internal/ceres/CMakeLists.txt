# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2010, 2011, 2012 Google Inc. All rights reserved.
# http://code.google.com/p/ceres-solver/
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
# Author: keir@google.com (Keir Mierle)

SET(CERES_INTERNAL_SRC
    array_utils.cc
    block_evaluate_preparer.cc
    block_jacobi_preconditioner.cc
    block_jacobian_writer.cc
    block_random_access_dense_matrix.cc
    block_random_access_matrix.cc
    block_random_access_sparse_matrix.cc
    block_sparse_matrix.cc
    block_structure.cc
    canonical_views_clustering.cc
    cgnr_solver.cc
    compressed_row_jacobian_writer.cc
    compressed_row_sparse_matrix.cc
    conditioned_cost_function.cc
    conjugate_gradients_solver.cc
    coordinate_descent_minimizer.cc
    corrector.cc
    cxsparse.cc
    dense_normal_cholesky_solver.cc
    dense_qr_solver.cc
    dense_sparse_matrix.cc
    detect_structure.cc
    dogleg_strategy.cc
    evaluator.cc
    file.cc
    gradient_checking_cost_function.cc
    implicit_schur_complement.cc
    iterative_schur_complement_solver.cc
    levenberg_marquardt_strategy.cc
    line_search.cc
    line_search_direction.cc
    line_search_minimizer.cc
    linear_least_squares_problems.cc
    linear_operator.cc
    linear_solver.cc
    local_parameterization.cc
    loss_function.cc
    low_rank_inverse_hessian.cc
    minimizer.cc
    normal_prior.cc
    parameter_block_ordering.cc
    partitioned_matrix_view.cc
    polynomial.cc
    problem.cc
    problem_impl.cc
    program.cc
    residual_block.cc
    residual_block_utils.cc
    runtime_numeric_diff_cost_function.cc
    schur_complement_solver.cc
    schur_eliminator.cc
    scratch_evaluate_preparer.cc
    solver.cc
    solver_impl.cc
    sparse_matrix.cc
    sparse_normal_cholesky_solver.cc
    split.cc
    stringprintf.cc
    suitesparse.cc
    triplet_sparse_matrix.cc
    trust_region_minimizer.cc
    trust_region_strategy.cc
    types.cc
    visibility.cc
    visibility_based_preconditioner.cc
    wall_time.cc
)

# Also depend on the header files so that they appear in IDEs.
FILE(GLOB CERES_INTERNAL_HDRS *.h)

# Include the specialized schur solvers.
IF (${SCHUR_SPECIALIZATIONS})
  FILE(GLOB CERES_INTERNAL_SCHUR_FILES generated/*.cc)
ELSE (${SCHUR_SPECIALIZATIONS})
  # Only the fully dynamic solver. The build is much faster this way.
  FILE(GLOB CERES_INTERNAL_SCHUR_FILES generated/schur_eliminator_d_d_d.cc)
ENDIF (${SCHUR_SPECIALIZATIONS})

# For Android, use the internal Glog implementation.
IF (${MINIGLOG})
  ADD_LIBRARY(miniglog STATIC
              miniglog/glog/logging.cc)

  IF (UNIX)
    SET_TARGET_PROPERTIES(miniglog PROPERTIES COMPILE_FLAGS -fPIC)
  ENDIF(UNIX)

  # The Android logging library that defines e.g. __android_log_print is
  # creatively named "log".
#  TARGET_LINK_LIBRARIES(miniglog log)

  INSTALL(TARGETS miniglog
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
ENDIF (${MINIGLOG})

SET(CERES_LIBRARY_DEPENDENCIES ${GLOG_LIB})

IF (${CXSPARSE_FOUND})
  LIST(APPEND CERES_LIBRARY_DEPENDENCIES cxsparse)
ENDIF (${CXSPARSE_FOUND})

IF (${OPENMP_FOUND})
  IF (NOT MSVC)
    LIST(APPEND CERES_LIBRARY_DEPENDENCIES gomp)
  ENDIF (NOT MSVC)
ENDIF (${OPENMP_FOUND})

SET(CERES_LIBRARY_SOURCE
    ${PROTO_SRCS}
    ${PROTO_HDRS}
    ${CERES_INTERNAL_SRC}
    ${CERES_INTERNAL_HDRS}
    ${CERES_INTERNAL_SCHUR_FILES})

ADD_LIBRARY(ceres STATIC ${CERES_LIBRARY_SOURCE})
TARGET_LINK_LIBRARIES(ceres ${CERES_LIBRARY_DEPENDENCIES})

INSTALL(TARGETS ceres
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
