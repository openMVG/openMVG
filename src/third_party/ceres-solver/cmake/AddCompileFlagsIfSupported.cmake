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
# Author: sergey.vfx@gmail.com (Sergey Sharybin)

function(add_cxx_compiler_flag_if_supported
    AGGREGATED_CXX_FLAGS_VAR
    FLAG_TO_ADD_IF_SUPPORTED)
  include(CheckCXXCompilerFlag)
  # Use of whitespace or '-' in variable names (used by CheckCXXSourceCompiles
  # as #defines) will trigger errors.
  string(STRIP "${FLAG_TO_ADD_IF_SUPPORTED}" FLAG_TO_ADD_IF_SUPPORTED)
  # Build an informatively named test result variable so that it will be evident
  # which tests were performed/succeeded in the CMake output, e.g for -Wall:
  #
  # -- Performing Test CHECK_CXX_FLAG_Wall - Success
  #
  # NOTE: This variable is also used to cache test result.
  string(REPLACE "-" "_" CHECK_CXX_FLAG
    "CHECK_CXX_FLAG${FLAG_TO_ADD_IF_SUPPORTED}")
  check_cxx_compiler_flag(${FLAG_TO_ADD_IF_SUPPORTED} ${CHECK_CXX_FLAG})
  if (${CHECK_CXX_FLAG})
    set(${AGGREGATED_CXX_FLAGS_VAR}
      "${${AGGREGATED_CXX_FLAGS_VAR}} ${FLAG_TO_ADD_IF_SUPPORTED}" PARENT_SCOPE)
  endif()
endfunction()
