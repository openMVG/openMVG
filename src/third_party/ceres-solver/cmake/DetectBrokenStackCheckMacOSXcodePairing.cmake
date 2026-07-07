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
# Author: alexs.mac@gmail.com (Alex Stewart)

# As detailed in [1] the combination of macOS 10.15.x (Catalina) and
# Xcode 11.0-1 enables by default a broken version of -fstack-check which
# can break the alignment requirements for SIMD instructions resulting in
# segfaults from within Eigen. This issue was apparently fixed in Xcode 11.2
# despite not appearing in the official release notes.
#
# Although this can be worked around by compiling with -fno-stack-check, we
# instead prevent generation as the update to Xcode 11.2 is free and failing
# to include -fno-stack-check *everywhere* could still result in random
# segfaults.
#
# [1]: https://forums.developer.apple.com/thread/121887
function(detect_broken_stack_check_macos_xcode_pairing)
  if (NOT APPLE)
    return()
  endif()

  execute_process(COMMAND sw_vers -productVersion
    OUTPUT_VARIABLE MACOS_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if (MACOS_VERSION VERSION_LESS 10.15)
    # Only 10.15 (Catalina) is likely to be affected, irrespective of the Xcode
    # version. Although it is possible to recreate the issue on 10.14 (Mojave)
    # and Xcode 11.0-1 if -fstack-check is forced on, this is not the default.
    return()
  endif()

  execute_process(COMMAND xcodebuild -version
    OUTPUT_VARIABLE XCODE_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX MATCH "Xcode [0-9\\.]+" XCODE_VERSION "${XCODE_VERSION}")
  string(REGEX REPLACE "Xcode ([0-9\\.]+)" "\\1" XCODE_VERSION "${XCODE_VERSION}")

  if ((XCODE_VERSION VERSION_EQUAL 11.0) OR
      (XCODE_VERSION VERSION_EQUAL 11.1))
    message(FATAL_ERROR "Detected macOS version: ${MACOS_VERSION} and "
      "Xcode version: ${XCODE_VERSION} which combined exhibit an "
      "-fstack-check bug which can break alignment requirements for at least "
      "AVX instructions as detailed here [1]."
      "\n"
      "This bug affected Xcode 11.0 and 11.1 but only when used with 10.15 "
      "(Catalina), and was fixed in Xcode 11.2. Without the fix in place, "
      "random segfaults will occur in Eigen operations used by Ceres that use "
      "AVX instructions."
      "\n"
      "Please update to at least Xcode 11.2."
      "\n"
      "[1]: https://forums.developer.apple.com/thread/121887")
  endif()
endfunction()
