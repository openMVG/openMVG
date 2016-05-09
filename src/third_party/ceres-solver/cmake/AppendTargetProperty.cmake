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

# Append item(s) to a property on a declared CMake target:
#
#    append_target_property(target property item_to_append1
#                                           [... item_to_appendN])
#
# The set_target_properties() CMake function will overwrite the contents of the
# specified target property.  This function instead appends to it, so can
# be called multiple times with the same target & property to iteratively
# populate it.
function(append_target_property TARGET PROPERTY)
  if (NOT TARGET ${TARGET})
    message(FATAL_ERROR "Invalid target: ${TARGET} cannot append: ${ARGN} "
      "to property: ${PROPERTY}")
  endif()
  if (NOT PROPERTY)
    message(FATAL_ERROR "Invalid property to update for target: ${TARGET}")
  endif()
  # Get the initial state of the specified property for the target s/t
  # we can append to it (not overwrite it).
  get_target_property(INITIAL_PROPERTY_STATE ${TARGET} ${PROPERTY})
  if (NOT INITIAL_PROPERTY_STATE)
    # Ensure that if the state is unset, we do not insert the XXX-NOTFOUND
    # returned by CMake into the property.
    set(INITIAL_PROPERTY_STATE "")
  endif()
  # Delistify (remove ; separators) the potentially set of items to append
  # to the specified target property.
  string(REPLACE ";" " " ITEMS_TO_APPEND "${ARGN}")
  set_target_properties(${TARGET} PROPERTIES ${PROPERTY}
    "${INITIAL_PROPERTY_STATE} ${ITEMS_TO_APPEND}")
endfunction()
