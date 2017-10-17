# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2017 Google Inc. All rights reserved.
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

# Microsoft deprecated the POSIX Bessel functions: j[0,1,n]() in favour
# of _j[0,1,n](), it appears since at least MSVC 2005 [1].  This function
# checks if the underscore prefixed versions of the Bessel functions are
# defined, and sets ${HAVE_UNDERSCORE_PREFIXED_BESSEL_FUNCTIONS_VAR} to
# TRUE if they do.
#
# [1] https://msdn.microsoft.com/en-us/library/ms235384(v=vs.100).aspx
function(check_if_underscore_prefixed_bessel_functions_exist
    HAVE_UNDERSCORE_PREFIXED_BESSEL_FUNCTIONS_VAR)
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles(
    "#include <math.h>
     int main(int argc, char * argv[]) {
       double result;
       result = _j0(1.2345);
       result = _j1(1.2345);
       result = _jn(2, 1.2345);
       return 0;
     }"
     HAVE_UNDERSCORE_PREFIXED_BESSEL_FUNCTIONS)
   set(${HAVE_UNDERSCORE_PREFIXED_BESSEL_FUNCTIONS_VAR}
     ${HAVE_UNDERSCORE_PREFIXED_BESSEL_FUNCTIONS}
     PARENT_SCOPE)
endfunction()
