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
# Author: keir@google.com (Keir Mierle)

set(CERES_INTERNAL_SRC
    array_utils.cc
    blas.cc
    block_evaluate_preparer.cc
    block_jacobi_preconditioner.cc
    block_jacobian_writer.cc
    block_random_access_dense_matrix.cc
    block_random_access_diagonal_matrix.cc
    block_random_access_matrix.cc
    block_random_access_sparse_matrix.cc
    block_sparse_matrix.cc
    block_structure.cc
    c_api.cc
    canonical_views_clustering.cc
    cgnr_solver.cc
    callbacks.cc
    compressed_col_sparse_matrix_utils.cc
    compressed_row_jacobian_writer.cc
    compressed_row_sparse_matrix.cc
    conditioned_cost_function.cc
    conjugate_gradients_solver.cc
    coordinate_descent_minimizer.cc
    corrector.cc
    covariance.cc
    covariance_impl.cc
    cxsparse.cc
    dense_normal_cholesky_solver.cc
    dense_qr_solver.cc
    dense_sparse_matrix.cc
    detect_structure.cc
    dogleg_strategy.cc
    dynamic_compressed_row_jacobian_writer.cc
    dynamic_compressed_row_sparse_matrix.cc
    evaluator.cc
    file.cc
    gradient_checker.cc
    gradient_checking_cost_function.cc
    gradient_problem.cc
    gradient_problem_solver.cc
    implicit_schur_complement.cc
    is_close.cc
    iterative_schur_complement_solver.cc
    levenberg_marquardt_strategy.cc
    lapack.cc
    line_search.cc
    line_search_direction.cc
    line_search_minimizer.cc
    line_search_preprocessor.cc
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
    preconditioner.cc
    preprocessor.cc
    problem.cc
    problem_impl.cc
    program.cc
    reorder_program.cc
    residual_block.cc
    residual_block_utils.cc
    schur_complement_solver.cc
    schur_eliminator.cc
    schur_jacobi_preconditioner.cc
    scratch_evaluate_preparer.cc
    single_linkage_clustering.cc
    solver.cc
    solver_utils.cc
    sparse_matrix.cc
    sparse_normal_cholesky_solver.cc
    split.cc
    stringprintf.cc
    suitesparse.cc
    triplet_sparse_matrix.cc
    trust_region_preprocessor.cc
    trust_region_minimizer.cc
    trust_region_step_evaluator.cc
    trust_region_strategy.cc
    types.cc
    visibility.cc
    visibility_based_preconditioner.cc
    wall_time.cc
)

# Also depend on the header files so that they appear in IDEs.
file(GLOB CERES_INTERNAL_HDRS *.h)

# Include the specialized schur solvers.
if (SCHUR_SPECIALIZATIONS)
  file(GLOB CERES_INTERNAL_SCHUR_FILES generated/*.cc)
else (SCHUR_SPECIALIZATIONS)
  # Only the fully dynamic solver. The build is much faster this way.
  file(GLOB CERES_INTERNAL_SCHUR_FILES generated/*_d_d_d.cc)
endif (SCHUR_SPECIALIZATIONS)

# Build the list of dependencies for Ceres based on the current configuration.
if (NOT MINIGLOG AND GLOG_FOUND)
  list(APPEND CERES_LIBRARY_PUBLIC_DEPENDENCIES ${GLOG_LIBRARIES})
  if (GFLAGS_FOUND)
    # If glog & gflags are both found, we assume that glog was built with
    # gflags, as it is awkward to perform a try_compile() to verify this
    # when gflags is an imported target (as it is in newer versions).
    # As glog #includes gflags/gflags.h in glog/logging.h if compiled with
    # gflags, it is thus a public dependency for Ceres in this case.
    list(APPEND CERES_LIBRARY_PUBLIC_DEPENDENCIES ${GFLAGS_LIBRARIES})
  endif()
endif (NOT MINIGLOG AND GLOG_FOUND)

if (SUITESPARSE AND SUITESPARSE_FOUND)
  # Define version information for use in Solver::FullReport.
  add_definitions(-DCERES_SUITESPARSE_VERSION="${SUITESPARSE_VERSION}")
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${SUITESPARSE_LIBRARIES})
endif (SUITESPARSE AND SUITESPARSE_FOUND)

if (CXSPARSE AND CXSPARSE_FOUND)
  # Define version information for use in Solver::FullReport.
  add_definitions(-DCERES_CXSPARSE_VERSION="${CXSPARSE_VERSION}")
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${CXSPARSE_LIBRARIES})
endif (CXSPARSE AND CXSPARSE_FOUND)

if (BLAS_FOUND AND LAPACK_FOUND)
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${LAPACK_LIBRARIES})
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${BLAS_LIBRARIES})
endif (BLAS_FOUND AND LAPACK_FOUND)

if (OPENMP_FOUND)
  # OpenMP support in Clang requires a non-GNU OpenMP library.
  if (CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES gomp)
  endif()
  if (NOT MSVC)
    list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${CMAKE_THREAD_LIBS_INIT})
  endif()
endif (OPENMP_FOUND)

set(CERES_LIBRARY_SOURCE
    ${CERES_INTERNAL_SRC}
    ${CERES_INTERNAL_HDRS}
    ${CERES_INTERNAL_SCHUR_FILES})

# Primarily for Android, but optionally for others, compile the minimal
# glog implementation into Ceres.
if (MINIGLOG)
  list(APPEND CERES_LIBRARY_SOURCE miniglog/glog/logging.cc)
endif (MINIGLOG)

add_library(ceres ${CERES_LIBRARY_SOURCE})
set_target_properties(ceres PROPERTIES
  VERSION ${CERES_VERSION}
  SOVERSION ${CERES_VERSION_MAJOR})

include(AppendTargetProperty)
# Always build position-independent code (PIC), even when building Ceres as a
# static library so that shared libraries can link against it, not just
# executables (PIC does not apply on Windows).
if (NOT WIN32 AND NOT BUILD_SHARED_LIBS)
  # Use the explicit POSITION_INDEPENDENT_CODE target property on CMake versions
  # that support it (>= 2.8.9).  Otherwise, manually add the -fPIC flag as an
  # additional compile definitions for the target.
  if (CMAKE_VERSION VERSION_LESS "2.8.9")
    append_target_property(ceres COMPILE_FLAGS "-fPIC")
  else()
    # Use set_target_properties() not append_target_property() here as
    # POSITION_INDEPENDENT_CODE is a binary ON/OFF switch.
    set_target_properties(ceres PROPERTIES POSITION_INDEPENDENT_CODE ON)
  endif()
endif()

if (CMAKE_VERSION VERSION_LESS "2.8.12")
  # CMake version < 2.8.12 does not support target_compile_options(), warn
  # user that they will have to add compile flags to their own projects
  # manually if required.
  if (CXX11 AND COMPILER_HAS_CXX11_FLAG)
    message("-- Warning: Detected CMake version: ${CMAKE_VERSION} < 2.8.12, "
      "which is the minimum required for compile options to be included in an "
      "exported CMake target, but CERES_USE_CXX11 is enabled and the detected. "
      "compiler requires -std=c++11. The client is responsible for adding "
      "-std=c++11 when using Ceres.")
  endif()
else()
  # CMake version >= 2.8.12 supports target_compile_options().
  if (CXX11 AND COMPILER_HAS_CXX11_FLAG)
    # If Ceres is compiled using C++11, and the compiler requires -std=c++11
    # to be set, then ensure that this requirement is rolled into the exported
    # CMake target, such that client code which uses Ceres will inherit it (if
    # the CMake version supports it), iff they are NOT compiling for C.  We
    # check for not C, rather than C++ as LINKER_LANGUAGE is often NOTFOUND and
    # then uses the default (C++).
    target_compile_options(ceres PUBLIC
      $<$<NOT:$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,C>>:-std=c++11>)
  endif()
endif()

if (BUILD_SHARED_LIBS)
  # When building a shared library, mark all external libraries as
  # PRIVATE so they don't show up as a dependency.
  target_link_libraries(ceres
        LINK_PUBLIC ${CERES_LIBRARY_PUBLIC_DEPENDENCIES}
        LINK_PRIVATE ${CERES_LIBRARY_PRIVATE_DEPENDENCIES})
else (BUILD_SHARED_LIBS)
  # When building a static library, all external libraries are
  # PUBLIC(default) since the user needs to link to them.
  # They will be listed in CeresTargets.cmake.
  set(CERES_LIBRARY_DEPENDENCIES
        ${CERES_LIBRARY_PUBLIC_DEPENDENCIES}
        ${CERES_LIBRARY_PRIVATE_DEPENDENCIES})
  target_link_libraries(ceres ${CERES_LIBRARY_DEPENDENCIES})
endif (BUILD_SHARED_LIBS)

install(TARGETS ceres
        EXPORT  CeresExport
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib${LIB_SUFFIX}
        ARCHIVE DESTINATION lib${LIB_SUFFIX})
INSTALL(TARGETS ceres
        EXPORT  openMVG-targets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib${LIB_SUFFIX}
        ARCHIVE DESTINATION lib${LIB_SUFFIX})

if (BUILD_TESTING AND GFLAGS)
  add_library(gtest gmock_gtest_all.cc gmock_main.cc)
  if (BUILD_SHARED_LIBS)
    # Define gtest-specific shared library flags for compilation.
    append_target_property(gtest COMPILE_DEFINITIONS
      GTEST_CREATE_SHARED_LIBRARY)
  endif()

  add_library(test_util
              evaluator_test_utils.cc
              numeric_diff_test_utils.cc
              test_util.cc)

  if (MINIGLOG)
    # When using miniglog, it is compiled into Ceres, thus Ceres becomes
    # the library against which other libraries should link for logging.
    target_link_libraries(gtest ${GFLAGS_LIBRARIES} ceres)
    target_link_libraries(test_util ceres gtest)
  else (MINIGLOG)
    target_link_libraries(gtest ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})
    target_link_libraries(test_util ceres gtest ${GLOG_LIBRARIES})
  endif (MINIGLOG)

  macro (CERES_TEST NAME)
    add_executable(${NAME}_test ${NAME}_test.cc)
    target_link_libraries(${NAME}_test test_util ceres gtest)
    if (BUILD_SHARED_LIBS)
      # Define gtest-specific shared library flags for linking.
      append_target_property(${NAME}_test COMPILE_DEFINITIONS
        GTEST_LINKED_AS_SHARED_LIBRARY)
    endif()
    add_test(NAME ${NAME}_test
             COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}_test
             --test_srcdir
             ${CMAKE_SOURCE_DIR}/data)
  endmacro (CERES_TEST)

  ceres_test(array_utils)
  ceres_test(autodiff)
  ceres_test(autodiff_cost_function)
  ceres_test(autodiff_local_parameterization)
  ceres_test(block_jacobi_preconditioner)
  ceres_test(block_random_access_dense_matrix)
  ceres_test(block_random_access_diagonal_matrix)
  ceres_test(block_random_access_sparse_matrix)
  ceres_test(block_sparse_matrix)
  ceres_test(bundle_adjustment)
  ceres_test(c_api)
  ceres_test(canonical_views_clustering)
  ceres_test(compressed_row_sparse_matrix)
  ceres_test(conditioned_cost_function)
  ceres_test(corrector)
  ceres_test(cost_function_to_functor)
  ceres_test(covariance)
  ceres_test(cubic_interpolation)
  ceres_test(detect_structure)
  ceres_test(dense_sparse_matrix)
  ceres_test(dynamic_autodiff_cost_function)
  ceres_test(dynamic_compressed_row_sparse_matrix)
  ceres_test(dynamic_numeric_diff_cost_function)
  ceres_test(evaluator)
  ceres_test(gradient_checker)
  ceres_test(gradient_checking_cost_function)
  ceres_test(gradient_problem)
  ceres_test(gradient_problem_solver)
  ceres_test(graph)
  ceres_test(graph_algorithms)
  ceres_test(householder_vector)
  ceres_test(is_close)
  ceres_test(implicit_schur_complement)
  ceres_test(iterative_schur_complement_solver)
  ceres_test(jet)
  ceres_test(levenberg_marquardt_strategy)
  ceres_test(dogleg_strategy)
  ceres_test(line_search_preprocessor)
  ceres_test(local_parameterization)
  ceres_test(loss_function)
  ceres_test(minimizer)
  ceres_test(normal_prior)
  ceres_test(numeric_diff_cost_function)
  ceres_test(ordered_groups)
  ceres_test(parameter_block)
  ceres_test(parameter_block_ordering)
  ceres_test(partitioned_matrix_view)
  ceres_test(polynomial)
  ceres_test(problem)
  ceres_test(program)
  ceres_test(reorder_program)
  ceres_test(residual_block)
  ceres_test(residual_block_utils)
  ceres_test(rotation)
  ceres_test(schur_complement_solver)
  ceres_test(schur_eliminator)
  ceres_test(single_linkage_clustering)
  ceres_test(small_blas)
  ceres_test(solver)

  # TODO(sameeragarwal): This test should ultimately be made
  # independent of SuiteSparse.
  if (SUITESPARSE AND SUITESPARSE_FOUND)
    ceres_test(compressed_col_sparse_matrix_utils)
    target_link_libraries(compressed_col_sparse_matrix_utils_test
                          ${SUITESPARSE_LIBRARIES})
  endif (SUITESPARSE AND SUITESPARSE_FOUND)

  ceres_test(symmetric_linear_solver)
  ceres_test(triplet_sparse_matrix)
  ceres_test(trust_region_minimizer)
  ceres_test(trust_region_preprocessor)
  ceres_test(unsymmetric_linear_solver)
  ceres_test(visibility)
  ceres_test(visibility_based_preconditioner)

  # Put the large end to end test last.
  ceres_test(system)
endif (BUILD_TESTING AND GFLAGS)
