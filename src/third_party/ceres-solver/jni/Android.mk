# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2015 Google Inc. All rights reserved.
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
# Author: settinger@google.com (Scott Ettinger)
#         keir@google.com (Keir Mierle)
#
# Builds Ceres for Android, using the standard toolchain (not
# standalone). It uses LLVM's libc++ as the standard library. It is a
# modern BSD licensed implementation of the standard c++ library. We
# do this to avoid any licensing issues that may arise from using
# GCC's libstdc++ which is licensed under GPL3.
#
# Building
# --------
#
# You will have to specify the environment EIGEN_PATH to point to the
# Eigen sources when building. For example:
#
#   EIGEN_PATH=/home/keir/src/eigen-3.0.5 ndk-build -j
#
# It is also possible to specify CERES_EXTRA_DEFINES, in case you need
# to pass more definitions to the C compiler.
#
# Using the library
# -----------------
# Copy the static library:
#
#   ../obj/local/armeabi-v7a/libceres.a
#
# into your own project, then link it into your binary in your
# Android.mk file.
#
# Reducing binary size
# --------------------
# This build includes the Schur specializations, which increase the
# size of the binary. If you don't need them for your application,
# consider adding:
#
#   -DCERES_RESTRICT_SCHUR_SPECIALIZATION
#
# to the LOCAL_CFLAGS variable below.
#
# Changing the logging library
# ----------------------------
# Ceres Solver ships with a replacement for glog that provides a
# simple and small implementation that builds on Android. However, if
# you wish to supply a header only version yourself, then you may
# define CERES_GLOG_DIR to point to it.

LOCAL_PATH := $(call my-dir)

# Ceres requires at least NDK version r9d to compile.
ifneq ($(shell $(LOCAL_PATH)/assert_ndk_version.sh "r9d" $(NDK_ROOT)), true)
  $(error Ceres requires NDK version r9d or greater)
endif

EIGEN_PATH := $(EIGEN_PATH)
CERES_INCLUDE_PATHS := $(CERES_EXTRA_INCLUDES)
CERES_INCLUDE_PATHS += $(LOCAL_PATH)/../internal
CERES_INCLUDE_PATHS += $(LOCAL_PATH)/../internal/ceres
CERES_INCLUDE_PATHS += $(LOCAL_PATH)/../include
CERES_INCLUDE_PATHS += $(LOCAL_PATH)/../config

# Use the alternate glog implementation if provided by the user.
ifdef CERES_GLOG_DIR
  CERES_INCLUDE_PATHS += $(CERES_GLOG_DIR)
else
  CERES_INCLUDE_PATHS += $(LOCAL_PATH)/../internal/ceres/miniglog
endif
CERES_SRC_PATH := ../internal/ceres

include $(CLEAR_VARS)
LOCAL_C_INCLUDES := $(CERES_INCLUDE_PATHS)
LOCAL_C_INCLUDES += $(EIGEN_PATH)

LOCAL_CPP_EXTENSION := .cc
LOCAL_CFLAGS := $(CERES_EXTRA_DEFINES) \
                -DCERES_NO_LAPACK \
                -DCERES_NO_SUITESPARSE \
                -DCERES_NO_CXSPARSE \
                -DCERES_STD_UNORDERED_MAP


# If the user did not enable threads in CERES_EXTRA_DEFINES, then add
# CERES_NO_THREADS.
#
# TODO(sameeragarwal): Update comments here and in the docs to
# demonstrate how OpenMP can be used by the user.
ifeq (,$(findstring CERES_HAVE_PTHREAD, $(LOCAL_CFLAGS)))
  LOCAL_CFLAGS += -DCERES_NO_THREADS
endif

LOCAL_SRC_FILES := $(CERES_SRC_PATH)/array_utils.cc \
                   $(CERES_SRC_PATH)/blas.cc \
                   $(CERES_SRC_PATH)/block_evaluate_preparer.cc \
                   $(CERES_SRC_PATH)/block_jacobian_writer.cc \
                   $(CERES_SRC_PATH)/block_jacobi_preconditioner.cc \
                   $(CERES_SRC_PATH)/block_random_access_dense_matrix.cc \
                   $(CERES_SRC_PATH)/block_random_access_diagonal_matrix.cc \
                   $(CERES_SRC_PATH)/block_random_access_matrix.cc \
                   $(CERES_SRC_PATH)/block_random_access_sparse_matrix.cc \
                   $(CERES_SRC_PATH)/block_sparse_matrix.cc \
                   $(CERES_SRC_PATH)/block_structure.cc \
                   $(CERES_SRC_PATH)/callbacks.cc \
                   $(CERES_SRC_PATH)/canonical_views_clustering.cc \
                   $(CERES_SRC_PATH)/cgnr_solver.cc \
                   $(CERES_SRC_PATH)/compressed_row_jacobian_writer.cc \
                   $(CERES_SRC_PATH)/compressed_row_sparse_matrix.cc \
                   $(CERES_SRC_PATH)/conditioned_cost_function.cc \
                   $(CERES_SRC_PATH)/conjugate_gradients_solver.cc \
                   $(CERES_SRC_PATH)/coordinate_descent_minimizer.cc \
                   $(CERES_SRC_PATH)/corrector.cc \
                   $(CERES_SRC_PATH)/covariance.cc \
                   $(CERES_SRC_PATH)/covariance_impl.cc \
                   $(CERES_SRC_PATH)/dense_normal_cholesky_solver.cc \
                   $(CERES_SRC_PATH)/dense_qr_solver.cc \
                   $(CERES_SRC_PATH)/dense_sparse_matrix.cc \
                   $(CERES_SRC_PATH)/detect_structure.cc \
                   $(CERES_SRC_PATH)/dogleg_strategy.cc \
                   $(CERES_SRC_PATH)/dynamic_compressed_row_jacobian_writer.cc \
                   $(CERES_SRC_PATH)/dynamic_compressed_row_sparse_matrix.cc \
                   $(CERES_SRC_PATH)/evaluator.cc \
                   $(CERES_SRC_PATH)/file.cc \
                   $(CERES_SRC_PATH)/gradient_checker.cc \
                   $(CERES_SRC_PATH)/gradient_checking_cost_function.cc \
                   $(CERES_SRC_PATH)/gradient_problem.cc \
                   $(CERES_SRC_PATH)/gradient_problem_solver.cc \
                   $(CERES_SRC_PATH)/is_close.cc \
                   $(CERES_SRC_PATH)/implicit_schur_complement.cc \
                   $(CERES_SRC_PATH)/iterative_schur_complement_solver.cc \
                   $(CERES_SRC_PATH)/lapack.cc \
                   $(CERES_SRC_PATH)/levenberg_marquardt_strategy.cc \
                   $(CERES_SRC_PATH)/line_search.cc \
                   $(CERES_SRC_PATH)/line_search_direction.cc \
                   $(CERES_SRC_PATH)/line_search_minimizer.cc \
                   $(CERES_SRC_PATH)/linear_least_squares_problems.cc \
                   $(CERES_SRC_PATH)/linear_operator.cc \
                   $(CERES_SRC_PATH)/line_search_preprocessor.cc \
                   $(CERES_SRC_PATH)/linear_solver.cc \
                   $(CERES_SRC_PATH)/local_parameterization.cc \
                   $(CERES_SRC_PATH)/loss_function.cc \
                   $(CERES_SRC_PATH)/low_rank_inverse_hessian.cc \
                   $(CERES_SRC_PATH)/minimizer.cc \
                   $(CERES_SRC_PATH)/normal_prior.cc \
                   $(CERES_SRC_PATH)/parameter_block_ordering.cc \
                   $(CERES_SRC_PATH)/partitioned_matrix_view.cc \
                   $(CERES_SRC_PATH)/polynomial.cc \
                   $(CERES_SRC_PATH)/preconditioner.cc \
                   $(CERES_SRC_PATH)/preprocessor.cc \
                   $(CERES_SRC_PATH)/problem.cc \
                   $(CERES_SRC_PATH)/problem_impl.cc \
                   $(CERES_SRC_PATH)/program.cc \
                   $(CERES_SRC_PATH)/reorder_program.cc \
                   $(CERES_SRC_PATH)/residual_block.cc \
                   $(CERES_SRC_PATH)/residual_block_utils.cc \
                   $(CERES_SRC_PATH)/schur_complement_solver.cc \
                   $(CERES_SRC_PATH)/schur_eliminator.cc \
                   $(CERES_SRC_PATH)/schur_jacobi_preconditioner.cc \
                   $(CERES_SRC_PATH)/scratch_evaluate_preparer.cc \
                   $(CERES_SRC_PATH)/solver.cc \
                   $(CERES_SRC_PATH)/solver_utils.cc \
                   $(CERES_SRC_PATH)/sparse_matrix.cc \
                   $(CERES_SRC_PATH)/sparse_normal_cholesky_solver.cc \
                   $(CERES_SRC_PATH)/split.cc \
                   $(CERES_SRC_PATH)/stringprintf.cc \
                   $(CERES_SRC_PATH)/suitesparse.cc \
                   $(CERES_SRC_PATH)/triplet_sparse_matrix.cc \
                   $(CERES_SRC_PATH)/trust_region_minimizer.cc \
                   $(CERES_SRC_PATH)/trust_region_preprocessor.cc \
                   $(CERES_SRC_PATH)/trust_region_step_evaluator.cc \
                   $(CERES_SRC_PATH)/trust_region_strategy.cc \
                   $(CERES_SRC_PATH)/types.cc \
                   $(CERES_SRC_PATH)/visibility_based_preconditioner.cc \
                   $(CERES_SRC_PATH)/visibility.cc \
                   $(CERES_SRC_PATH)/wall_time.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_d_d_d.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_2_2.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_2_3.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_2_4.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_2_d.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_3_3.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_3_4.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_3_6.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_3_9.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_3_d.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_4_3.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_4_4.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_4_8.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_4_9.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_4_d.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_2_d_d.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_4_4_2.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_4_4_3.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_4_4_4.cc \
                   $(CERES_SRC_PATH)/generated/schur_eliminator_4_4_d.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_d_d_d.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_2_2.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_2_3.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_2_4.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_2_d.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_3_3.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_3_4.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_3_6.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_3_9.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_3_d.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_4_3.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_4_4.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_4_8.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_4_9.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_4_d.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_2_d_d.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_4_4_2.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_4_4_3.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_4_4_4.cc \
                   $(CERES_SRC_PATH)/generated/partitioned_matrix_view_4_4_d.cc

ifndef CERES_GLOG_DIR
LOCAL_SRC_FILES += $(CERES_SRC_PATH)/miniglog/glog/logging.cc
endif

LOCAL_MODULE := ceres
include $(BUILD_STATIC_LIBRARY)
