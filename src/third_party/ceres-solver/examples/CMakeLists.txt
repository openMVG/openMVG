# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2022 Google Inc. All rights reserved.
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

add_executable(helloworld helloworld.cc)
target_link_libraries(helloworld PRIVATE Ceres::ceres)

add_executable(helloworld_numeric_diff helloworld_numeric_diff.cc)
target_link_libraries(helloworld_numeric_diff PRIVATE Ceres::ceres)

add_executable(helloworld_analytic_diff helloworld_analytic_diff.cc)
target_link_libraries(helloworld_analytic_diff PRIVATE Ceres::ceres)

add_executable(curve_fitting curve_fitting.cc)
target_link_libraries(curve_fitting PRIVATE Ceres::ceres)

add_executable(rosenbrock rosenbrock.cc)
target_link_libraries(rosenbrock PRIVATE Ceres::ceres)

add_executable(rosenbrock_analytic_diff rosenbrock_analytic_diff.cc)
target_link_libraries(rosenbrock_analytic_diff PRIVATE Ceres::ceres)

add_executable(rosenbrock_numeric_diff rosenbrock_numeric_diff.cc)
target_link_libraries(rosenbrock_numeric_diff PRIVATE Ceres::ceres)

add_executable(curve_fitting_c curve_fitting.c)
target_link_libraries(curve_fitting_c PRIVATE Ceres::ceres)
# Force CMake to link curve_fitting_c using the C++ linker.
set_target_properties(curve_fitting_c PROPERTIES LINKER_LANGUAGE CXX)
# As this is a C file #including <math.h> we have to explicitly add the math
# library (libm). Although some compilers (dependent upon options) will accept
# the indirect link to libm via Ceres, at least GCC 4.8 on pure Debian won't.
if (HAVE_LIBM)
  target_link_libraries(curve_fitting_c PRIVATE m)
endif (HAVE_LIBM)

add_executable(ellipse_approximation ellipse_approximation.cc)
target_link_libraries(ellipse_approximation PRIVATE Ceres::ceres)

add_executable(robust_curve_fitting robust_curve_fitting.cc)
target_link_libraries(robust_curve_fitting PRIVATE Ceres::ceres)

add_executable(simple_bundle_adjuster simple_bundle_adjuster.cc)
target_link_libraries(simple_bundle_adjuster PRIVATE Ceres::ceres)

add_executable(bicubic_interpolation bicubic_interpolation.cc)
target_link_libraries(bicubic_interpolation PRIVATE Ceres::ceres)

add_executable(bicubic_interpolation_analytic bicubic_interpolation_analytic.cc)
target_link_libraries(bicubic_interpolation_analytic PRIVATE Ceres::ceres)

add_executable(iteration_callback_example iteration_callback_example.cc)
target_link_libraries(iteration_callback_example PRIVATE Ceres::ceres)

add_executable(evaluation_callback_example evaluation_callback_example.cc)
target_link_libraries(evaluation_callback_example PRIVATE Ceres::ceres)

if (GFLAGS)
  add_executable(powell powell.cc)
  target_link_libraries(powell PRIVATE Ceres::ceres gflags)

  add_executable(nist nist.cc)
  target_link_libraries(nist PRIVATE Ceres::ceres gflags)
  if (HAVE_BIGOBJ)
    target_compile_options(nist PRIVATE /bigobj)
  endif()

  add_executable(more_garbow_hillstrom more_garbow_hillstrom.cc)
  target_link_libraries(more_garbow_hillstrom PRIVATE Ceres::ceres gflags)

  add_executable(circle_fit circle_fit.cc)
  target_link_libraries(circle_fit PRIVATE Ceres::ceres gflags)

  add_executable(bundle_adjuster
                 bundle_adjuster.cc
                 bal_problem.cc)
  target_link_libraries(bundle_adjuster PRIVATE Ceres::ceres gflags)

  add_executable(libmv_bundle_adjuster
                 libmv_bundle_adjuster.cc)
  target_link_libraries(libmv_bundle_adjuster PRIVATE Ceres::ceres gflags)

  add_executable(libmv_homography
                 libmv_homography.cc)
  target_link_libraries(libmv_homography PRIVATE Ceres::ceres gflags)

  add_executable(denoising
                 denoising.cc
                 fields_of_experts.cc)
  target_link_libraries(denoising PRIVATE Ceres::ceres gflags)

  add_executable(robot_pose_mle
                 robot_pose_mle.cc)
  target_link_libraries(robot_pose_mle PRIVATE Ceres::ceres gflags)
endif (GFLAGS)

add_subdirectory(sampled_function)
add_subdirectory(slam)
