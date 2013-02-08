# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2012 Google Inc. All rights reserved.
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
# Author: arnaudgelas@gmail.com (Arnaud Gelas)

FIND_PROGRAM(PYGMENTIZE NAMES pygmentize
    PATHS
        /usr/bin
        /usr/local/bin
        /opt/bin
        /opt/local/bin)

IF (${PYGMENTIZE} MATCHES PYGMENTIZE-NOTFOUND)
  MESSAGE(SEND_ERROR "Pygmentize is not installed")
ENDIF (${PYGMENTIZE} MATCHES PYGMENTIZE-NOTFOUND)

FIND_PACKAGE(LATEX REQUIRED)

FILE(DOWNLOAD http://minted.googlecode.com/files/minted.sty
     ${CMAKE_BINARY_DIR}/docs/minted.sty
     SHOW_PROGRESS)

ADD_CUSTOM_TARGET(UserGuide ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/docs ${CMAKE_BINARY_DIR}/docs
    COMMAND ${PDFLATEX_COMPILER}
        -shell-escape ${CMAKE_BINARY_DIR}/docs/ceres-solver.tex
        -output-directory ${CMAKE_BINARY_DIR}/docs
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/docs
    COMMAND ${BIBTEX_COMPILER} ceres-solver.aux
    COMMAND ${PDFLATEX_COMPILER}
        -shell-escape ${CMAKE_BINARY_DIR}/docs/ceres-solver.tex
        -output-directory ${CMAKE_BINARY_DIR}/docs
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/docs
    COMMAND ${PDFLATEX_COMPILER}
        -shell-escape ${CMAKE_BINARY_DIR}/docs/ceres-solver.tex
        -output-directory ${CMAKE_BINARY_DIR}/docs
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/docs)

INSTALL(FILES ${CMAKE_BINARY_DIR}/docs/ceres-solver.pdf
        DESTINATION share/ceres/docs)
