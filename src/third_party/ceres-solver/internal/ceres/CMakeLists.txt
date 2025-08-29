# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2022 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form %Cmust reproduce the above copyright notice,
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

# Build the list of dependencies for Ceres based on the current configuration.

# Avoid 'xxx.cc has no symbols' warnings from source files which are 'empty'
# when their enclosing #ifdefs are disabled.
find_package(Threads REQUIRED)
list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES Threads::Threads)
# Make dependency visible to the parent CMakeLists.txt
set(Threads_DEPENDENCY "find_dependency (Threads)" PARENT_SCOPE)

# Source files that contain public symbols and live in the ceres namespaces.
# Such symbols are expected to be marked with CERES_EXPORT and the files below
# sorted in lexicographical order.
set(CERES_EXPORTED_SRCS
    c_api.cc
    gradient_checker.cc
    gradient_problem.cc
    loss_function.cc
    manifold.cc
    normal_prior.cc
    problem.cc
    solver.cc
    types.cc
)

# Also depend on the header files so that they appear in IDEs.
file(GLOB CERES_INTERNAL_HDRS *.h)
if (MINIGLOG)
  file(GLOB MINIGLOG_HDRS miniglog/glog/*.h)
  list(APPEND CERES_INTERNAL_HDRS ${MINIGLOG_HDRS})
  if (ANDROID)
    list(APPEND CERES_LIBRARY_PUBLIC_DEPENDENCIES log)
  endif()
endif()

# Depend also on public headers so they appear in IDEs.
file(GLOB CERES_PUBLIC_HDRS ${Ceres_SOURCE_DIR}/include/ceres/*.h)
file(GLOB CERES_PUBLIC_INTERNAL_HDRS ${Ceres_SOURCE_DIR}/include/ceres/internal/*.h)
file(GLOB CERES_PUBLIC_INTERNAL_HDRS
  ${Ceres_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/ceres/internal/*.h)

# Include the specialized schur solvers.
if (SCHUR_SPECIALIZATIONS)
  file(GLOB CERES_INTERNAL_SCHUR_FILES generated/*.cc)
else (SCHUR_SPECIALIZATIONS)
  # Only the fully dynamic solver. The build is much faster this way.
  file(GLOB CERES_INTERNAL_SCHUR_FILES generated/*_d_d_d.cc)
endif (SCHUR_SPECIALIZATIONS)

# The generated specializations of the Schur eliminator include
# schur_eliminator_impl.h which defines EIGEN_CACHEFRIENDLY_PRODUCT_THRESHOLD
# to a different value than Eigen's default.  Depending on the order of files
# in the unity build this can lead to clashes.  Additionally, these files are
# already generated in a way which leads to fairly large compilation units,
# so the gains from a unity build would be marginal.
set_source_files_properties(${CERES_INTERNAL_SCHUR_FILES} PROPERTIES
  SKIP_UNITY_BUILD_INCLUSION ON)

if (NOT MINIGLOG AND GLOG_FOUND)
  list(APPEND CERES_LIBRARY_PUBLIC_DEPENDENCIES ${GLOG_LIBRARIES})
  if (gflags_FOUND)
    # If glog & gflags are both found, we assume that glog was built with
    # gflags, as it is awkward to perform a try_compile() to verify this
    # when gflags is an imported target (as it is in newer versions).
    # As glog #includes gflags/gflags.h in glog/logging.h if compiled with
    # gflags, it is thus a public dependency for Ceres in this case.
    list(APPEND CERES_LIBRARY_PUBLIC_DEPENDENCIES gflags)
  endif()
endif (NOT MINIGLOG AND GLOG_FOUND)

if (SUITESPARSE AND SuiteSparse_FOUND)
  # Define version information for use in Solver::FullReport.
  add_definitions(-DCERES_SUITESPARSE_VERSION="${SuiteSparse_VERSION}")
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES SuiteSparse::CHOLMOD
    SuiteSparse::SPQR)

  if (SuiteSparse_Partition_FOUND)
    list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES SuiteSparse::Partition)
  endif (SuiteSparse_Partition_FOUND)
endif (SUITESPARSE AND SuiteSparse_FOUND)

if (SuiteSparse_Partition_FOUND OR EIGENMETIS)
  # Define version information for use in Solver::FullReport.
  add_definitions(-DCERES_METIS_VERSION="${METIS_VERSION}")
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES METIS::METIS)
endif (SuiteSparse_Partition_FOUND OR EIGENMETIS)

if (ACCELERATESPARSE AND AccelerateSparse_FOUND)
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${AccelerateSparse_LIBRARIES})
endif()

if (USE_CUDA)
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES
    CUDA::cublas
    CUDA::cudart
    CUDA::cusolver
    CUDA::cusparse)
  set_source_files_properties(cuda_kernels_vector_ops.cu.cc PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(cuda_kernels_bsm_to_crs.cu.cc PROPERTIES LANGUAGE CUDA)
  add_library(ceres_cuda_kernels STATIC cuda_kernels_vector_ops.cu.cc cuda_kernels_bsm_to_crs.cu.cc)
  target_compile_features(ceres_cuda_kernels PRIVATE cxx_std_14)
  target_compile_definitions(ceres_cuda_kernels PRIVATE CERES_STATIC_DEFINE)
  # Enable __host__ / __device__ annotations in lambda declarations
  target_compile_options(ceres_cuda_kernels PRIVATE --extended-lambda)
  target_include_directories(ceres_cuda_kernels PRIVATE  ${Ceres_SOURCE_DIR}/include  ${Ceres_SOURCE_DIR}/internal ${Ceres_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ceres_cuda_kernels)
endif (USE_CUDA)

if (LAPACK_FOUND)
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES ${LAPACK_LIBRARIES})
endif ()

# Source files that contain private symbols and live in the ceres::internal
# namespace. The corresponding symbols (classes, functions, etc.) are expected
# to be marked with CERES_NO_EXPORT and the files below sorted in
# lexicographical order.
add_library(ceres_internal OBJECT
    ${CERES_INTERNAL_SCHUR_FILES}
    accelerate_sparse.cc
    array_utils.cc
    block_evaluate_preparer.cc
    block_jacobi_preconditioner.cc
    block_jacobian_writer.cc
    block_random_access_dense_matrix.cc
    block_random_access_diagonal_matrix.cc
    block_random_access_matrix.cc
    block_random_access_sparse_matrix.cc
    block_sparse_matrix.cc
    block_structure.cc
    callbacks.cc
    canonical_views_clustering.cc
    cgnr_solver.cc
    compressed_col_sparse_matrix_utils.cc
    compressed_row_jacobian_writer.cc
    compressed_row_sparse_matrix.cc
    conditioned_cost_function.cc
    context.cc
    context_impl.cc
    coordinate_descent_minimizer.cc
    corrector.cc
    cost_function.cc
    covariance.cc
    covariance_impl.cc
    cuda_block_sparse_crs_view.cc
    cuda_partitioned_block_sparse_crs_view.cc
    cuda_block_structure.cc
    cuda_sparse_matrix.cc
    cuda_vector.cc
    dense_cholesky.cc
    dense_normal_cholesky_solver.cc
    dense_qr.cc
    dense_qr_solver.cc
    dense_sparse_matrix.cc
    detect_structure.cc
    dogleg_strategy.cc
    dynamic_compressed_row_jacobian_writer.cc
    dynamic_compressed_row_sparse_matrix.cc
    dynamic_sparse_normal_cholesky_solver.cc
    eigensparse.cc
    evaluation_callback.cc
    evaluator.cc
    fake_bundle_adjustment_jacobian.cc
    file.cc
    first_order_function.cc
    float_suitesparse.cc
    function_sample.cc
    gradient_checking_cost_function.cc
    gradient_problem_solver.cc
    implicit_schur_complement.cc
    inner_product_computer.cc
    is_close.cc
    iteration_callback.cc
    iterative_refiner.cc
    iterative_schur_complement_solver.cc
    levenberg_marquardt_strategy.cc
    line_search.cc
    line_search_direction.cc
    line_search_minimizer.cc
    line_search_preprocessor.cc
    linear_least_squares_problems.cc
    linear_operator.cc
    linear_solver.cc
    low_rank_inverse_hessian.cc
    minimizer.cc
    parallel_invoke.cc
    parallel_utils.cc
    parallel_vector_ops.cc
    parameter_block_ordering.cc
    partitioned_matrix_view.cc
    polynomial.cc
    power_series_expansion_preconditioner.cc
    preconditioner.cc
    preprocessor.cc
    problem_impl.cc
    program.cc
    reorder_program.cc
    residual_block.cc
    residual_block_utils.cc
    schur_complement_solver.cc
    schur_eliminator.cc
    schur_jacobi_preconditioner.cc
    schur_templates.cc
    scratch_evaluate_preparer.cc
    single_linkage_clustering.cc
    solver_utils.cc
    sparse_cholesky.cc
    sparse_matrix.cc
    sparse_normal_cholesky_solver.cc
    stringprintf.cc
    subset_preconditioner.cc
    suitesparse.cc
    thread_pool.cc
    thread_token_provider.cc
    triplet_sparse_matrix.cc
    trust_region_minimizer.cc
    trust_region_preprocessor.cc
    trust_region_step_evaluator.cc
    trust_region_strategy.cc
    visibility.cc
    visibility_based_preconditioner.cc
    wall_time.cc
)

set(CERES_LIBRARY_SOURCE
    ${CERES_EXPORTED_SRCS}
    ${CERES_INTERNAL_HDRS}
    ${CERES_PUBLIC_HDRS}
    ${CERES_PUBLIC_INTERNAL_HDRS})

# Primarily for Android, but optionally for others, compile the minimal
# glog implementation into Ceres.
if (MINIGLOG)
  list(APPEND CERES_LIBRARY_SOURCE miniglog/glog/logging.cc)
endif (MINIGLOG)

# Ceres C++ compiler flags can be too strict for an external library code
# which we do not maintain.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-missing-declarations"
                        CHECK_CXX_FLAG_Wno_missing_declarations)
if (CHECK_CXX_FLAG_Wno_missing_declarations)
  set_property(SOURCE gmock_gtest_all.cc
               APPEND_STRING PROPERTY COMPILE_FLAGS "-Wno-missing-declarations")
endif()

add_library(ceres $<TARGET_OBJECTS:ceres_internal> ${CERES_LIBRARY_SOURCE})

if(BUILD_SHARED_LIBS)
  # While building shared libraries, we additionally require a static variant to
  # be able to access internal symbols which are not intended for general use.
  # Therefore, create a static library from object files and apply all the
  # compiler options from the main library to the static one.
  add_library(ceres_static STATIC $<TARGET_OBJECTS:ceres_internal> ${CERES_LIBRARY_SOURCE})
  target_include_directories(ceres_static PUBLIC $<TARGET_PROPERTY:ceres,INCLUDE_DIRECTORIES>)
  target_compile_definitions(ceres_static PUBLIC $<TARGET_PROPERTY:ceres,COMPILE_DEFINITIONS>)
  target_compile_features(ceres_static PUBLIC $<TARGET_PROPERTY:ceres,COMPILE_FEATURES>)
  target_compile_options(ceres_static PUBLIC $<TARGET_PROPERTY:ceres,COMPILE_OPTIONS>)
  target_link_libraries(ceres_static
    INTERFACE $<TARGET_PROPERTY:ceres,INTERFACE_LINK_LIBRARIES>
    PRIVATE ${CERES_LIBRARY_PRIVATE_DEPENDENCIES})
  # CERES_STATIC_DEFINE is generated by the GenerateExportHeader CMake module
  # used to autogerate export.h. The macro should not be renamed without
  # updating the corresponding generate_export_header invocation.
  target_compile_definitions(ceres_static PUBLIC CERES_STATIC_DEFINE)
else()
  # In a static library build, not additional access layer is necessary as all
  # symbols are visible.
  add_library(ceres_static ALIAS ceres)
endif()

# Create a local alias target that matches the expected installed target.
add_library(Ceres::ceres ALIAS ceres)

# Apply all compiler options from the main Ceres target. Compiler options should
# be generally defined on the main target referenced by the ceres_target CMake
# variable.
target_include_directories(ceres_internal PUBLIC $<TARGET_PROPERTY:ceres,INCLUDE_DIRECTORIES>)
target_compile_definitions(ceres_internal PUBLIC $<TARGET_PROPERTY:ceres,COMPILE_DEFINITIONS>)
target_compile_options(ceres_internal PUBLIC $<TARGET_PROPERTY:ceres,COMPILE_OPTIONS>)
target_compile_definitions(ceres_internal PRIVATE ceres_EXPORTS)
target_compile_features(ceres_internal PRIVATE $<TARGET_PROPERTY:ceres,COMPILE_FEATURES>)

# Ensure the minimum required C++ language version is fulfilled as our
# requirement by downstream clients. Consumers can choose the same or a newer
# language standard revision.
target_compile_features(ceres PUBLIC cxx_std_17)

set_target_properties(ceres PROPERTIES
  VERSION ${CERES_VERSION}
  SOVERSION 4)

target_link_libraries(ceres
  PUBLIC ${CERES_LIBRARY_PUBLIC_DEPENDENCIES}
  PRIVATE ${CERES_LIBRARY_PRIVATE_DEPENDENCIES})

# Add the Ceres headers to its target.
#
# Force the location containing the configured config.h to the front of the
# include_directories list (by default it is appended to the back) to ensure
# that if the user has an installed version of Ceres in the same location as one
# of the dependencies (e.g. /usr/local) that we find the config.h we just
# configured, not the (older) installed config.h.
target_include_directories(ceres
  BEFORE PUBLIC
  $<BUILD_INTERFACE:${Ceres_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${Ceres_SOURCE_DIR}/internal
  PUBLIC $<BUILD_INTERFACE:${Ceres_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Eigen SparseQR generates various compiler warnings related to unused and
# uninitialised local variables.  To avoid having to individually suppress these
# warnings around the #include statements for Eigen headers across all GCC/Clang
# versions, we tell CMake to treat Eigen headers as system headers.  This
# results in all compiler warnings from them being suppressed.
target_link_libraries(ceres PUBLIC Eigen3::Eigen)

# Gather the list of public & private include locations for all enabled optional
# dependencies to be added to the Ceres target.
set(CERES_LIBRARY_PRIVATE_DEPENDENCIES_INCLUDE_DIRS "")
set(CERES_LIBRARY_PUBLIC_DEPENDENCIES_INCLUDE_DIRS "")
if (MINIGLOG)
  # Force the miniglog headers to the front of the public include directories
  # to protect against the case when the user has glog installed in a standard
  # location (specifically the same as the Ceres install location) but compiled
  # Ceres with MINIGLOG anyway.  Otherwise: "glog/logging.h" in the public Ceres
  # headers used in client code would match the installed version of glog, not
  # the miniglog headers, and the client application would fail to link.
  #
  # Note that this is an imperfect fix, as we cannot control the include
  # directories in client projects, and they could easily invert this ordering
  # themselves (intentionally or otherwise) and so break their build.
  target_include_directories(ceres BEFORE PUBLIC
    $<BUILD_INTERFACE:${Ceres_SOURCE_DIR}/internal/ceres/miniglog>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ceres/internal/miniglog>)
elseif (NOT FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION)
  # Only append glog include directories if the glog found was not a CMake
  # exported target that already includes them.
  list(APPEND CERES_LIBRARY_PUBLIC_DEPENDENCIES_INCLUDE_DIRS
    ${GLOG_INCLUDE_DIRS})
endif()
if (ACCELERATESPARSE)
  list(APPEND CERES_LIBRARY_PRIVATE_DEPENDENCIES_INCLUDE_DIRS
    ${AccelerateSparse_INCLUDE_DIRS})
endif()
# Add include locations for optional dependencies to the Ceres target without
# duplication.
list(REMOVE_DUPLICATES CERES_LIBRARY_PRIVATE_DEPENDENCIES_INCLUDE_DIRS)
target_include_directories(ceres PRIVATE ${CERES_LIBRARY_PRIVATE_DEPENDENCIES_INCLUDE_DIRS})
list(REMOVE_DUPLICATES CERES_LIBRARY_PUBLIC_DEPENDENCIES_INCLUDE_DIRS)
target_include_directories(ceres PUBLIC ${CERES_LIBRARY_PUBLIC_DEPENDENCIES_INCLUDE_DIRS})

# Generate an export header for annotating symbols visibility
include(GenerateExportHeader)
generate_export_header(ceres EXPORT_FILE_NAME
  ${Ceres_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/ceres/internal/export.h)

if (USE_CUDA)
  install(TARGETS ceres_cuda_kernels
          EXPORT  CeresExport
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(USE_CUDA)

install(TARGETS ceres
        EXPORT  CeresExport
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if (BUILD_TESTING AND GFLAGS)
  add_library(gtest STATIC gmock_gtest_all.cc gmock_main.cc)

  target_include_directories(gtest PRIVATE ${Ceres_SOURCE_DIR}/internal/ceres)
  if (CMAKE_SYSTEM_NAME MATCHES "QNX")
    target_link_libraries(gtest PUBLIC regex)
  endif()
  target_link_libraries(gtest PRIVATE Ceres::ceres gflags)

  add_library(test_util STATIC
              evaluator_test_utils.cc
              numeric_diff_test_utils.cc
              test_util.cc)

  target_include_directories(test_util PUBLIC ${Ceres_SOURCE_DIR}/internal)
  target_link_libraries (test_util PUBLIC ceres_static gflags gtest)

  macro (CERES_TEST NAME)
    add_executable(${NAME}_test ${NAME}_test.cc)
    # Pull in local headers from the generated test directories when ceres_test()
    # is invoked there, as well as the private headers in this directory which
    # may be referenced without the 'ceres' path prefix and all private
    # dependencies that may be directly referenced.
    target_include_directories(${NAME}_test
      PRIVATE ${Ceres_SOURCE_DIR}/internal/ceres
              ${CERES_LIBRARY_PRIVATE_DEPENDENCIES_INCLUDE_DIRS})
    # Some tests include direct references/includes of private dependency
    # headers which are not propagated via the ceres targets, so link them
    # explicitly.
    target_link_libraries(${NAME}_test PRIVATE gtest test_util ceres_static
      ${CERES_LIBRARY_PRIVATE_DEPENDENCIES})

    # covariance_test uses SuiteSparseQR.hpp. However, since SuiteSparse import
    # targets are private (link only) dependencies not propagated to consumers,
    # we need to link against the target explicitly here.
    if (TARGET SuiteSparse::SPQR)
      target_link_libraries (${NAME}_test PRIVATE SuiteSparse::SPQR)
    endif (TARGET SuiteSparse::SPQR)

    add_test(NAME ${NAME}_test
             COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}_test
             --test_srcdir
             ${Ceres_SOURCE_DIR}/data)
  endmacro (CERES_TEST)

  ceres_test(array_utils)
  ceres_test(array_selector)
  ceres_test(autodiff)
  ceres_test(autodiff_first_order_function)
  ceres_test(autodiff_cost_function)
  ceres_test(autodiff_manifold)
  ceres_test(block_jacobi_preconditioner)
  ceres_test(block_random_access_dense_matrix)
  ceres_test(block_random_access_diagonal_matrix)
  ceres_test(block_random_access_sparse_matrix)
  ceres_test(block_sparse_matrix)
  ceres_test(c_api)
  ceres_test(canonical_views_clustering)
  ceres_test(compressed_col_sparse_matrix_utils)
  ceres_test(compressed_row_sparse_matrix)
  ceres_test(concurrent_queue)
  ceres_test(conditioned_cost_function)
  ceres_test(conjugate_gradients_solver)
  ceres_test(corrector)
  ceres_test(cost_function_to_functor)
  ceres_test(covariance)
  ceres_test(cubic_interpolation)
  ceres_test(cuda_partitioned_block_sparse_crs_view)
  ceres_test(cuda_block_sparse_crs_view)
  ceres_test(cuda_block_structure)
  ceres_test(cuda_dense_cholesky)
  ceres_test(cuda_dense_qr)
  ceres_test(cuda_kernels_vector_ops)
  ceres_test(cuda_sparse_matrix)
  ceres_test(cuda_streamed_buffer)
  ceres_test(cuda_vector)
  ceres_test(dense_linear_solver)
  ceres_test(dense_cholesky)
  ceres_test(dense_qr)
  ceres_test(dense_sparse_matrix)
  ceres_test(detect_structure)
  ceres_test(dogleg_strategy)
  ceres_test(dynamic_autodiff_cost_function)
  ceres_test(dynamic_compressed_row_sparse_matrix)
  ceres_test(dynamic_numeric_diff_cost_function)
  ceres_test(dynamic_sparse_normal_cholesky_solver)
  ceres_test(dynamic_sparsity)
  ceres_test(evaluation_callback)
  ceres_test(evaluator)
  ceres_test(fixed_array)
  ceres_test(gradient_checker)
  ceres_test(gradient_checking_cost_function)
  ceres_test(gradient_problem)
  ceres_test(gradient_problem_solver)
  ceres_test(graph)
  ceres_test(graph_algorithms)
  ceres_test(householder_vector)
  ceres_test(implicit_schur_complement)
  ceres_test(inner_product_computer)
  ceres_test(invert_psd_matrix)
  ceres_test(integer_sequence_algorithm)
  ceres_test(is_close)
  ceres_test(iterative_refiner)
  ceres_test(iterative_schur_complement_solver)
  ceres_test(jet)
  ceres_test(jet_traits)
  ceres_test(levenberg_marquardt_strategy)
  ceres_test(line_search_minimizer)
  ceres_test(line_search_preprocessor)
  ceres_test(loss_function)
  ceres_test(manifold)
  ceres_test(minimizer)
  ceres_test(normal_prior)
  ceres_test(numeric_diff_cost_function)
  ceres_test(numeric_diff_first_order_function)
  ceres_test(ordered_groups)
  ceres_test(parallel_for)
  ceres_test(parallel_utils)
  ceres_test(parameter_block)
  ceres_test(parameter_block_ordering)
  ceres_test(parameter_dims)
  ceres_test(partitioned_matrix_view)
  ceres_test(polynomial)
  ceres_test(power_series_expansion_preconditioner)
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
  ceres_test(sparse_cholesky)
  ceres_test(sparse_normal_cholesky_solver)
  ceres_test(subset_preconditioner)
  ceres_test(system)
  ceres_test(tiny_solver)
  ceres_test(tiny_solver_autodiff_function)
  ceres_test(tiny_solver_cost_function_adapter)
  ceres_test(thread_pool)
  ceres_test(triplet_sparse_matrix)
  ceres_test(trust_region_minimizer)
  ceres_test(trust_region_preprocessor)
  ceres_test(visibility)
  ceres_test(visibility_based_preconditioner)

  add_subdirectory(generated_bundle_adjustment_tests)

endif (BUILD_TESTING AND GFLAGS)

macro(add_dependencies_to_benchmark BENCHMARK_TARGET)
  target_include_directories(${BENCHMARK_TARGET}
    PRIVATE ${Ceres_SOURCE_DIR}/internal
            ${CERES_LIBRARY_PRIVATE_DEPENDENCIES_INCLUDE_DIRS})
  # Benchmarks include direct references/includes of private dependency headers
  # which are not propagated via the ceres targets, so link them explicitly.
  target_link_libraries(${BENCHMARK_TARGET}
    PRIVATE benchmark::benchmark ceres_static
            ${CERES_LIBRARY_PRIVATE_DEPENDENCIES})
endmacro()

if (BUILD_BENCHMARKS)
  add_executable(evaluation_benchmark evaluation_benchmark.cc)
  add_dependencies_to_benchmark(evaluation_benchmark)

  add_executable(small_blas_gemv_benchmark small_blas_gemv_benchmark.cc)
  add_dependencies_to_benchmark(small_blas_gemv_benchmark)

  add_executable(small_blas_gemm_benchmark small_blas_gemm_benchmark.cc)
  add_dependencies_to_benchmark(small_blas_gemm_benchmark)

  add_executable(invert_psd_matrix_benchmark invert_psd_matrix_benchmark.cc)
  add_dependencies_to_benchmark(invert_psd_matrix_benchmark)

  add_executable(schur_eliminator_benchmark schur_eliminator_benchmark.cc)
  add_dependencies_to_benchmark(schur_eliminator_benchmark)

  add_executable(jet_operator_benchmark jet_operator_benchmark.cc)
  add_dependencies_to_benchmark(jet_operator_benchmark)

  add_executable(dense_linear_solver_benchmark dense_linear_solver_benchmark.cc)
  add_dependencies_to_benchmark(dense_linear_solver_benchmark)

  add_executable(spmv_benchmark spmv_benchmark.cc)
  add_dependencies_to_benchmark(spmv_benchmark)

  add_executable(parallel_vector_operations_benchmark parallel_vector_operations_benchmark.cc)
  add_dependencies_to_benchmark(parallel_vector_operations_benchmark)

  add_executable(parallel_for_benchmark parallel_for_benchmark.cc)
  add_dependencies_to_benchmark(parallel_for_benchmark)

  add_executable(block_jacobi_preconditioner_benchmark
    block_jacobi_preconditioner_benchmark.cc)
  add_dependencies_to_benchmark(block_jacobi_preconditioner_benchmark)

  add_subdirectory(autodiff_benchmarks)
endif (BUILD_BENCHMARKS)
