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
# Author: pablo.speciale@gmail.com (Pablo Speciale)
#

#[=======================================================================[.rst:
FindSphinx
==========

Module for locating Sphinx and its components.

This modules defines the following variables:

``Sphinx_FOUND``
  ``TRUE`` iff Sphinx and all of its components have been found.

``Sphinx_BUILD_EXECUTABLE``
  Path to the ``sphinx-build`` tool.
]=======================================================================]

include (FindPackageHandleStandardArgs)

find_program (Sphinx_BUILD_EXECUTABLE
  NAMES sphinx-build
  PATHS /opt/local/bin
  DOC "Sphinx documentation generator"
)

mark_as_advanced (Sphinx_BUILD_EXECUTABLE)

if (Sphinx_BUILD_EXECUTABLE)
  execute_process (
    COMMAND ${Sphinx_BUILD_EXECUTABLE} --version
    ERROR_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE _Sphinx_BUILD_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE _Sphinx_VERSION_STRING
    RESULT_VARIABLE _Sphinx_BUILD_RESULT
  )

  if (_Sphinx_BUILD_RESULT EQUAL 0)
    string (REGEX REPLACE "^sphinx-build[ \t]+([^ \t]+)$" "\\1" Sphinx_VERSION
      "${_Sphinx_VERSION_STRING}")

    if (Sphinx_VERSION MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
      set (Sphinx_VERSION_COMPONENTS 3)
      set (Sphinx_VERSION_MAJOR ${CMAKE_MATCH_1})
      set (Sphinx_VERSION_MINOR ${CMAKE_MATCH_2})
      set (Sphinx_VERSION_PATCH ${CMAKE_MATCH_3})
    endif (Sphinx_VERSION MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
  else (_Sphinx_BUILD_RESULT EQUAL 0)
    message (WARNING "Could not determine sphinx-build version: ${_Sphinx_BUILD_ERROR}")
  endif (_Sphinx_BUILD_RESULT EQUAL 0)

  unset (_Sphinx_BUILD_ERROR)
  unset (_Sphinx_BUILD_RESULT)
  unset (_Sphinx_VERSION_STRING)

  find_package (Python COMPONENTS Interpreter)
  set (_Sphinx_BUILD_RESULT FALSE)

  if (Python_Interpreter_FOUND)
    # Check for Sphinx theme dependency for documentation
    foreach (component IN LISTS Sphinx_FIND_COMPONENTS)
      string (REGEX MATCH "^(.+_theme)$" theme_component "${component}")

      if (NOT theme_component STREQUAL component)
        continue ()
      endif (NOT theme_component STREQUAL component)

      execute_process (
        COMMAND ${Python_EXECUTABLE} -c "import ${theme_component}"
        ERROR_STRIP_TRAILING_WHITESPACE
        ERROR_VARIABLE _Sphinx_BUILD_ERROR
        OUTPUT_QUIET
        RESULT_VARIABLE _Sphinx_BUILD_RESULT
      )

      if (_Sphinx_BUILD_RESULT EQUAL 0)
        set (Sphinx_${component}_FOUND TRUE)
      elseif (_Sphinx_BUILD_RESULT EQUAL 0)
        message (WARNING "Could not determine whether Sphinx component '${theme_component}' is available: ${_Sphinx_BUILD_ERROR}")
        set (Sphinx_${component}_FOUND FALSE)
      endif (_Sphinx_BUILD_RESULT EQUAL 0)

      unset (_Sphinx_BUILD_ERROR)
      unset (_Sphinx_BUILD_RESULT)
    endforeach (component)

    unset (theme_component)
  endif (Python_Interpreter_FOUND)
endif (Sphinx_BUILD_EXECUTABLE)

find_package_handle_standard_args (Sphinx
  REQUIRED_VARS Sphinx_BUILD_EXECUTABLE
  VERSION_VAR Sphinx_VERSION
  HANDLE_COMPONENTS
)
