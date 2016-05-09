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
# Author: alexs.mac@gmail.com (Alex Stewart)
#

# FindUnorderedMap.cmake - Find unordered_map header and namespace.
#
# This module defines the following variables:
#
# UNORDERED_MAP_FOUND: TRUE if unordered_map is found.
# HAVE_UNORDERED_MAP_IN_STD_NAMESPACE: Use <unordered_map> & std.
# HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE: Use <unordered_map> & std::tr1.
# HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE: <tr1/unordered_map> & std::tr1.
macro(FIND_UNORDERED_MAP)
  # To support CXX11 option, clear the results of all check_xxx() functions
  # s/t we always perform the checks each time, otherwise CMake fails to
  # detect that the tests should be performed again after CXX11 is toggled.
  unset(HAVE_STD_UNORDERED_MAP_HEADER CACHE)
  unset(HAVE_UNORDERED_MAP_IN_STD_NAMESPACE CACHE)
  unset(HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE CACHE)
  unset(HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE CACHE)

  set(UNORDERED_MAP_FOUND FALSE)
  include(CheckIncludeFileCXX)
  check_include_file_cxx(unordered_map HAVE_STD_UNORDERED_MAP_HEADER)
  if (HAVE_STD_UNORDERED_MAP_HEADER)
    # Finding the unordered_map header doesn't mean that unordered_map
    # is in std namespace.
    #
    # In particular, MSVC 2008 has unordered_map declared in std::tr1.
    # In order to support this, we do an extra check to see which
    # namespace should be used.
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("#include <unordered_map>
                               int main() {
                                 std::unordered_map<int, int> map;
                                 return 0;
                               }"
                               HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
    if (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
      set(UNORDERED_MAP_FOUND TRUE)
      message("-- Found unordered_map/set in std namespace.")
    else (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
      check_cxx_source_compiles("#include <unordered_map>
                                 int main() {
                                   std::tr1::unordered_map<int, int> map;
                                   return 0;
                                 }"
                                 HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
      if (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
        set(UNORDERED_MAP_FOUND TRUE)
        message("-- Found unordered_map/set in std::tr1 namespace.")
      else (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
        message("-- Found <unordered_map> but cannot find either "
          "std::unordered_map or std::tr1::unordered_map.")
      endif (HAVE_UNORDERED_MAP_IN_TR1_NAMESPACE)
    endif (HAVE_UNORDERED_MAP_IN_STD_NAMESPACE)
  else (HAVE_STD_UNORDERED_MAP_HEADER)
    check_include_file_cxx("tr1/unordered_map"
      HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE)
    if (HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE)
      set(UNORDERED_MAP_FOUND TRUE)
      message("-- Found tr1/unordered_map/set in std::tr1 namespace.")
    else (HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE)
      message("-- Unable to find <unordered_map> or <tr1/unordered_map>.")
    endif (HAVE_TR1_UNORDERED_MAP_IN_TR1_NAMESPACE)
  endif (HAVE_STD_UNORDERED_MAP_HEADER)
endmacro(FIND_UNORDERED_MAP)
