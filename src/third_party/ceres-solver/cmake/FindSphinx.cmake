# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2013 Google Inc. All rights reserved.
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
# Author: pablo.speciale@gmail.com (Pablo Speciale)
#

# Find the Sphinx documentation generator
#
# This modules defines
#  SPHINX_EXECUTABLE
#  SPHINX_FOUND

FIND_PROGRAM(SPHINX_EXECUTABLE
             NAMES sphinx-build
             PATHS
               /usr/bin
               /usr/local/bin
               /opt/local/bin
             DOC "Sphinx documentation generator")

IF (NOT SPHINX_EXECUTABLE)
  SET(_Python_VERSIONS 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0 1.6 1.5)

  FOREACH (_version ${_Python_VERSIONS})
    SET(_sphinx_NAMES sphinx-build-${_version})

    FIND_PROGRAM(SPHINX_EXECUTABLE
                 NAMES ${_sphinx_NAMES}
                 PATHS
                   /usr/bin
                   /usr/local/bin
                   /opt/loca/bin
                 DOC "Sphinx documentation generator")
  ENDFOREACH ()
ENDIF ()

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)

MARK_AS_ADVANCED(SPHINX_EXECUTABLE)
