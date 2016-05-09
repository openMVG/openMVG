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
# Author: sergey.vfx@gmail.com (Sergey Sharybin)
#

# FindSharedPtr.cmake - Find shared pointer header and namespace.
#
# This module defines the following variables:
#
# SHARED_PTR_FOUND: TRUE if shared_ptr found.
# SHARED_PTR_TR1_MEMORY_HEADER: True if <tr1/memory> header is to be used
# for the shared_ptr object, otherwise use <memory>.
# SHARED_PTR_TR1_NAMESPACE: TRUE if shared_ptr is defined in std::tr1 namespace,
# otherwise it's assumed to be defined in std namespace.

macro(FIND_SHARED_PTR)
  # To support CXX11 option, clear the results of all check_xxx() functions
  # s/t we always perform the checks each time, otherwise CMake fails to
  # detect that the tests should be performed again after CXX11 is toggled.
  unset(HAVE_STD_MEMORY_HEADER CACHE)
  unset(HAVE_SHARED_PTR_IN_STD_NAMESPACE CACHE)
  unset(HAVE_SHARED_PTR_IN_TR1_NAMESPACE CACHE)
  unset(HAVE_TR1_MEMORY_HEADER CACHE)
  unset(HAVE_SHARED_PTR_IN_TR1_NAMESPACE_FROM_TR1_MEMORY_HEADER CACHE)

  set(SHARED_PTR_FOUND FALSE)
  check_include_file_cxx(memory HAVE_STD_MEMORY_HEADER)
  if (HAVE_STD_MEMORY_HEADER)
    # Finding the memory header doesn't mean that shared_ptr is in std
    # namespace.
    #
    # In particular, MSVC 2008 has shared_ptr declared in std::tr1.  In
    # order to support this, we do an extra check to see which namespace
    # should be used.
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("#include <memory>
                               int main() {
                                 std::shared_ptr<int> int_ptr;
                                 return 0;
                               }"
                              HAVE_SHARED_PTR_IN_STD_NAMESPACE)

    if (HAVE_SHARED_PTR_IN_STD_NAMESPACE)
      message("-- Found shared_ptr in std namespace using <memory> header.")
      set(SHARED_PTR_FOUND TRUE)
    else (HAVE_SHARED_PTR_IN_STD_NAMESPACE)
      check_cxx_source_compiles("#include <memory>
                                 int main() {
                                   std::tr1::shared_ptr<int> int_ptr;
                                   return 0;
                                 }"
                                HAVE_SHARED_PTR_IN_TR1_NAMESPACE)
      if (HAVE_SHARED_PTR_IN_TR1_NAMESPACE)
        message("-- Found shared_ptr in std::tr1 namespace using <memory> header.")
        set(SHARED_PTR_TR1_NAMESPACE TRUE)
        set(SHARED_PTR_FOUND TRUE)
      endif (HAVE_SHARED_PTR_IN_TR1_NAMESPACE)
    endif (HAVE_SHARED_PTR_IN_STD_NAMESPACE)
  endif (HAVE_STD_MEMORY_HEADER)

  if (NOT SHARED_PTR_FOUND)
    # Further, gcc defines shared_ptr in std::tr1 namespace and
    # <tr1/memory> is to be included for this. And what makes things
    # even more tricky is that gcc does have <memory> header, so
    # all the checks above wouldn't find shared_ptr.
    check_include_file_cxx("tr1/memory" HAVE_TR1_MEMORY_HEADER)
    if (HAVE_TR1_MEMORY_HEADER)
      check_cxx_source_compiles("#include <tr1/memory>
                                 int main() {
                                   std::tr1::shared_ptr<int> int_ptr;
                                   return 0;
                                 }"
                                HAVE_SHARED_PTR_IN_TR1_NAMESPACE_FROM_TR1_MEMORY_HEADER)
      if (HAVE_SHARED_PTR_IN_TR1_NAMESPACE_FROM_TR1_MEMORY_HEADER)
        message("-- Found shared_ptr in std::tr1 namespace using <tr1/memory> header.")
          set(SHARED_PTR_TR1_MEMORY_HEADER TRUE)
          set(SHARED_PTR_TR1_NAMESPACE TRUE)
        set(SHARED_PTR_FOUND TRUE)
      endif (HAVE_SHARED_PTR_IN_TR1_NAMESPACE_FROM_TR1_MEMORY_HEADER)
    endif (HAVE_TR1_MEMORY_HEADER)
  endif (NOT SHARED_PTR_FOUND)
endmacro(FIND_SHARED_PTR)
