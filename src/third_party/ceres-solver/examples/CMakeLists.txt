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

# Only Ceres itself should be compiled with CERES_BUILDING_SHARED_LIBRARY
# defined, any users of Ceres will have CERES_USING_SHARED_LIBRARY defined
# for them in Ceres' config.h if appropriate.
IF (BUILD_SHARED_LIBS)
  REMOVE_DEFINITIONS(-DCERES_BUILDING_SHARED_LIBRARY)
ENDIF()

ADD_EXECUTABLE(helloworld helloworld.cc)
TARGET_LINK_LIBRARIES(helloworld ceres)

ADD_EXECUTABLE(helloworld_numeric_diff helloworld_numeric_diff.cc)
TARGET_LINK_LIBRARIES(helloworld_numeric_diff ceres)

ADD_EXECUTABLE(helloworld_analytic_diff helloworld_analytic_diff.cc)
TARGET_LINK_LIBRARIES(helloworld_analytic_diff ceres)

ADD_EXECUTABLE(curve_fitting curve_fitting.cc)
TARGET_LINK_LIBRARIES(curve_fitting ceres)

ADD_EXECUTABLE(rosenbrock rosenbrock.cc)
TARGET_LINK_LIBRARIES(rosenbrock ceres)

ADD_EXECUTABLE(curve_fitting_c curve_fitting.c)
TARGET_LINK_LIBRARIES(curve_fitting_c ceres)
# As this is a C file #including <math.h> we have to explicitly add the math
# library (libm). Although some compilers (dependent upon options) will accept
# the indirect link to libm via Ceres, at least GCC 4.8 on pure Debian won't.
IF (NOT MSVC)
  TARGET_LINK_LIBRARIES(curve_fitting_c m)
ENDIF (NOT MSVC)

ADD_EXECUTABLE(ellipse_approximation ellipse_approximation.cc)
TARGET_LINK_LIBRARIES(ellipse_approximation ceres)

ADD_EXECUTABLE(robust_curve_fitting robust_curve_fitting.cc)
TARGET_LINK_LIBRARIES(robust_curve_fitting ceres)

ADD_EXECUTABLE(simple_bundle_adjuster
               simple_bundle_adjuster.cc)
TARGET_LINK_LIBRARIES(simple_bundle_adjuster ceres)

IF (GFLAGS)
  # The CERES_GFLAGS_NAMESPACE compile definition is NOT stored in
  # CERES_COMPILE_OPTIONS (and thus config.h) as Ceres itself does not
  # require gflags, only the tests and examples do.
  ADD_DEFINITIONS(-DCERES_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})

  ADD_EXECUTABLE(powell powell.cc)
  TARGET_LINK_LIBRARIES(powell ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(nist nist.cc)
  TARGET_LINK_LIBRARIES(nist ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(more_garbow_hillstrom more_garbow_hillstrom.cc)
  TARGET_LINK_LIBRARIES(more_garbow_hillstrom ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(circle_fit circle_fit.cc)
  TARGET_LINK_LIBRARIES(circle_fit ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(bundle_adjuster
                 bundle_adjuster.cc
                 bal_problem.cc)
  TARGET_LINK_LIBRARIES(bundle_adjuster ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(libmv_bundle_adjuster
                 libmv_bundle_adjuster.cc)
  TARGET_LINK_LIBRARIES(libmv_bundle_adjuster ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(libmv_homography
                 libmv_homography.cc)
  TARGET_LINK_LIBRARIES(libmv_homography ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(denoising
                 denoising.cc
                 fields_of_experts.cc)
  TARGET_LINK_LIBRARIES(denoising ceres ${GFLAGS_LIBRARIES})

  ADD_EXECUTABLE(robot_pose_mle
                 robot_pose_mle.cc)
  TARGET_LINK_LIBRARIES(robot_pose_mle ceres ${GFLAGS_LIBRARIES})

ENDIF (GFLAGS)
