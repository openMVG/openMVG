# TODO: Add support for other compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Increase the inlining threshold only for those functions marked with an
  # inline hint. This is typically far more realistic to significantly increase
  # in a large code-base than -inline-threshold as that has a larger scope.
  list(APPEND CERES_BENCHMARK_FLAGS "-mllvm" "-inlinehint-threshold=1000000")
endif()

add_executable(autodiff_benchmarks autodiff_benchmarks.cc)
add_dependencies_to_benchmark(autodiff_benchmarks)
target_compile_options(autodiff_benchmarks PRIVATE ${CERES_BENCHMARK_FLAGS})
