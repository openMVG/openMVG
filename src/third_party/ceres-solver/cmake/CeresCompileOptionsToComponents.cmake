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
# Authors: alexs.mac@gmail.com (Alex Stewart)
#

# Conditionally add a value to the output list based on whether the specified
# value is found in the input list.
function(update_output_if_found INPUT_LIST_VAR OUTPUT_LIST_VAR ITEM_TO_FIND VAR_TO_COPY_IF_FOUND VAR_TO_COPY_IF_NOT_FOUND)
  list(FIND ${INPUT_LIST_VAR} "${ITEM_TO_FIND}" HAVE_ITEM)
  # list(FIND ..) returns -1 if the element was not in the list, but CMake
  # interprets if (VAR) to be true if VAR is any non-zero number, even
  # negative ones, hence we have to explicitly check for >= 0.
  if (HAVE_ITEM GREATER -1)
    list(APPEND ${OUTPUT_LIST_VAR} "${VAR_TO_COPY_IF_FOUND}")
  else()
    list(APPEND ${OUTPUT_LIST_VAR} "${VAR_TO_COPY_IF_NOT_FOUND}")
  endif()
  set(${OUTPUT_LIST_VAR} ${${OUTPUT_LIST_VAR}} PARENT_SCOPE)
endfunction()

# Helpers for update_output_if_found() to improve legibility when dealing with
# USE_XXX & NO_XXX option types in ceres_compile_options_to_components().
macro(add_to_output_if_found INPUT_LIST_VAR OUTPUT_LIST_VAR ITEM_TO_FIND VAR_TO_COPY_IF_FOUND)
  update_output_if_found(${INPUT_LIST_VAR}
    ${OUTPUT_LIST_VAR}
    "${ITEM_TO_FIND}"
    "${VAR_TO_COPY_IF_FOUND}"
    "") # Copy nothing if not found.
endmacro()

macro(add_to_output_if_not_found INPUT_LIST_VAR OUTPUT_LIST_VAR ITEM_TO_FIND VAR_TO_COPY_IF_NOT_FOUND)
  update_output_if_found(${INPUT_LIST_VAR}
    ${OUTPUT_LIST_VAR}
    "${ITEM_TO_FIND}"
    "" # Copy nothing if found
    "${VAR_TO_COPY_IF_NOT_FOUND}")
endmacro()

# Convert the Ceres compile options specified by: CURRENT_CERES_COMPILE_OPTIONS
# into the corresponding list of Ceres components (names), which may be used in:
# find_package(Ceres COMPONENTS <XXX>).
function(ceres_compile_options_to_components CURRENT_CERES_COMPILE_OPTIONS CERES_COMPONENTS_VAR)
  # To enable users to specify that they want *a* sparse linear algebra backend
  # without having to specify explicitly which one, for each sparse library we
  # add the 'meta-module': SparseLinearAlgebraLibrary in addition to their own
  # module name.
  add_to_output_if_found(CURRENT_CERES_COMPILE_OPTIONS ${CERES_COMPONENTS_VAR}
    CERES_USE_EIGEN_SPARSE "EigenSparse;SparseLinearAlgebraLibrary")
  add_to_output_if_not_found(CURRENT_CERES_COMPILE_OPTIONS ${CERES_COMPONENTS_VAR}
    CERES_NO_LAPACK "LAPACK")
  add_to_output_if_not_found(CURRENT_CERES_COMPILE_OPTIONS ${CERES_COMPONENTS_VAR}
    CERES_NO_SUITESPARSE "SuiteSparse;SparseLinearAlgebraLibrary")
  add_to_output_if_not_found(CURRENT_CERES_COMPILE_OPTIONS ${CERES_COMPONENTS_VAR}
    CERES_NO_ACCELERATE_SPARSE "AccelerateSparse;SparseLinearAlgebraLibrary")
  add_to_output_if_not_found(CURRENT_CERES_COMPILE_OPTIONS ${CERES_COMPONENTS_VAR}
    CERES_RESTRICT_SCHUR_SPECIALIZATION "SchurSpecializations")
  # Remove duplicates of SparseLinearAlgebraLibrary if multiple sparse backends
  # are present.
  list(REMOVE_DUPLICATES ${CERES_COMPONENTS_VAR})
  set(${CERES_COMPONENTS_VAR} "${${CERES_COMPONENTS_VAR}}" PARENT_SCOPE)
endfunction()
