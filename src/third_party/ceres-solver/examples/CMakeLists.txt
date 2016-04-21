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

# Only Ceres itself should be compiled with CERES_BUILDING_SHARED_LIBRARY
# defined, any users of Ceres will have CERES_USING_SHARED_LIBRARY defined
# for them in Ceres' config.h if appropriate.
if (BUILD_SHARED_LIBS)
  remove_definitions(-DCERES_BUILDING_SHARED_LIBRARY)
endif()

add_executable(helloworld helloworld.cc)
target_link_libraries(helloworld ceres)

add_executable(helloworld_numeric_diff helloworld_numeric_diff.cc)
target_link_libraries(helloworld_numeric_diff ceres)

add_executable(helloworld_analytic_diff helloworld_analytic_diff.cc)
target_link_libraries(helloworld_analytic_diff ceres)

add_executable(curve_fitting curve_fitting.cc)
target_link_libraries(curve_fitting ceres)

add_executable(rosenbrock rosenbrock.cc)
target_link_libraries(rosenbrock ceres)

add_executable(curve_fitting_c curve_fitting.c)
target_link_libraries(curve_fitting_c ceres)
# Force CMake to link curve_fitting_c using the C linker, this is important
# when Ceres was compiled using C++11 to ensure that -std=c++11 is not passed
# through.
set_target_properties(curve_fitting_c PROPERTIES LINKER_LANGUAGE C)
# As this is a C file #including <math.h> we have to explicitly add the math
# library (libm). Although some compilers (dependent upon options) will accept
# the indirect link to libm via Ceres, at least GCC 4.8 on pure Debian won't.
if (NOT MSVC)
  target_link_libraries(curve_fitting_c m)
endif (NOT MSVC)

add_executable(ellipse_approximation ellipse_approximation.cc)
target_link_libraries(ellipse_approximation ceres)

add_executable(robust_curve_fitting robust_curve_fitting.cc)
target_link_libraries(robust_curve_fitting ceres)

add_executable(simple_bundle_adjuster simple_bundle_adjuster.cc)
target_link_libraries(simple_bundle_adjuster ceres)

add_executable(sampled_function sampled_function.cc)
target_link_libraries(sampled_function ceres)

if (GFLAGS)
  # The CERES_GFLAGS_NAMESPACE compile definition is NOT stored in
  # CERES_COMPILE_OPTIONS (and thus config.h) as Ceres itself does not
  # require gflags, only the tests and examples do.
  add_definitions(-DCERES_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})

  add_executable(powell powell.cc)
  target_link_libraries(powell ceres ${GFLAGS_LIBRARIES})

  add_executable(nist nist.cc)
  target_link_libraries(nist ceres ${GFLAGS_LIBRARIES})

  add_executable(more_garbow_hillstrom more_garbow_hillstrom.cc)
  target_link_libraries(more_garbow_hillstrom ceres ${GFLAGS_LIBRARIES})

  add_executable(circle_fit circle_fit.cc)
  target_link_libraries(circle_fit ceres ${GFLAGS_LIBRARIES})

  add_executable(bundle_adjuster
                 bundle_adjuster.cc
                 bal_problem.cc)
  target_link_libraries(bundle_adjuster ceres ${GFLAGS_LIBRARIES})

  add_executable(libmv_bundle_adjuster
                 libmv_bundle_adjuster.cc)
  target_link_libraries(libmv_bundle_adjuster ceres ${GFLAGS_LIBRARIES})

  add_executable(libmv_homography
                 libmv_homography.cc)
  target_link_libraries(libmv_homography ceres ${GFLAGS_LIBRARIES})

  add_executable(denoising
                 denoising.cc
                 fields_of_experts.cc)
  target_link_libraries(denoising ceres ${GFLAGS_LIBRARIES})

  add_executable(robot_pose_mle
                 robot_pose_mle.cc)
  target_link_libraries(robot_pose_mle ceres ${GFLAGS_LIBRARIES})

endif (GFLAGS)
