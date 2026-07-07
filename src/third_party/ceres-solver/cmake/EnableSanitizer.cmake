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

# Usage: enable_sanitizer(REQUIRED_SANITIZERS) where REQUIRED_SANITIZERS should
# contain the list of sanitizers to enable by updating CMAKE_CXX_FLAGS and
# CMAKE_EXE_LINKER_FLAGS.
#
# The specified sanitizers will be checked both for compatibility with the
# current compiler and with each other as some sanitizers are mutually
# exclusive.
macro(enable_sanitizer)
  # According to the Clang documentation [1] the following sanitizers are
  # mututally exclusive.
  # [1]: https://clang.llvm.org/docs/UsersManual.html#controlling-code-generation
  set(INCOMPATIBLE_SANITIZERS address thread memory)
  # Set the recommended additional common compile flags for any sanitizer to
  # get the best possible output, e.g [2] but make them visible in the cache
  # so that the user can edit them if required.
  # [2]: https://clang.llvm.org/docs/AddressSanitizer.html#usage
  set(COMMON_SANITIZER_COMPILE_OPTIONS
    "-g -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Common compile flags enabled for any sanitizer")

  # Check that the specified list of sanitizers to enable does not include
  # multiple entries from the incompatible list.
  set(MERGED_SANITIZERS ${ARGN} ${INCOMPATIBLE_SANITIZERS})
  list(LENGTH MERGED_SANITIZERS COMBINED_LENGTH)
  list(REMOVE_DUPLICATES MERGED_SANITIZERS)
  list(LENGTH MERGED_SANITIZERS COMBINED_LENGTH_NO_DUPLICATES)
  math(EXPR VALID_LENGTH "${COMBINED_LENGTH} - 1")
  if (COMBINED_LENGTH_NO_DUPLICATES LESS VALID_LENGTH)
    include(PrettyPrintCMakeList)
    pretty_print_cmake_list(REQUESTED_SANITIZERS ${ARGN})
    pretty_print_cmake_list(
      PRETTY_INCOMPATIBLE_SANITIZERS ${INCOMPATIBLE_SANITIZERS})
    message(FATAL_ERROR "Found incompatible sanitizers in requested set: "
      "${REQUESTED_SANITIZERS}. The following sanitizers are mutually "
      "exclusive: ${PRETTY_INCOMPATIBLE_SANITIZERS}")
  endif()

  # Until CMake 3.14 and CMAKE_REQUIRED_LINK_OPTIONS there was no equivalent to
  # CMAKE_REQUIRED_FLAGS for try_compile() for linker flags. However, in CMake
  # 3.2 CMP0056 was introduced that when enabled passes CMAKE_EXE_LINKER_FLAGS
  # to try_compile() which allows us to achieve the same effect.
  cmake_policy(SET CMP0056 NEW)
  include(CheckCXXCompilerFlag)

  unset(ADDED_SANITIZER)
  foreach(REQUESTED_SANITIZER ${ARGN})
    set(SANITIZER_FLAG -fsanitize=${REQUESTED_SANITIZER})
    # Save the current CMAKE_EXE_LINKER_FLAGS before modifying it to test for
    # the existence of the sanitizer flag so that we can revert after the test.
    set(INITIAL_CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZER_FLAG}")
    check_cxx_compiler_flag(${SANITIZER_FLAG} HAVE_SANITIZER)
    set(CMAKE_EXE_LINKER_FLAGS "${INITIAL_CMAKE_EXE_LINKER_FLAGS}")
    if (NOT HAVE_SANITIZER)
      message(FATAL_ERROR "Specified sanitizer: ${REQUESTED_SANITIZER} is not "
        "supported by the compiler.")
    endif()
    message(STATUS "Enabling sanitizer: ${REQUESTED_SANITIZER}")
    set(ADDED_SANITIZER TRUE)
    # As per the Clang documentation, the sanitizer flags must be added to both
    # the compiler and linker flags.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAG}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZER_FLAG}")
  endforeach()
  if (ADDED_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_SANITIZER_COMPILE_OPTIONS}")
  endif()
endmacro()
