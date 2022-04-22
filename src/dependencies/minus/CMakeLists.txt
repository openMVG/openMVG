#  Copyright (c) 2022 Pierre MOULON, Ricardo Fabbri and Gabriel Andrade
#
# \author Pierre MOULON
# \author Gabriel ANDRADE Rio de Janeiro State U.
# \author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io)

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# -- MINUS: OpenMVG specific build style --------------------------------------
# optional, header-only library can be used without this library
# library is useful when the function will be used in many places
# to avoid having to include Minus many times and slowdown the compiler
#
# For OpenMVG it we do this precompile so that MINUS is always compiled with
# highly optimized compiler flags, even when OpenMVG is debug
set(minus_sources
  minus/minus/minus.hxx minus/minus/minus.h

  minus/minus/chicago14a.h minus/minus/chicago-default.h minus/minus/chicago14a-default-data.h minus/minus/chicago14a-io.h
  minus/minus/chicago14a.hxx minus/minus/chicago14a-default-data.hxx

  minus/minus/problem-defs.h

  minus/minus/Templates/minus-chicago-alltypes-allformulations+double.cxx
  minus/minus/Templates/minus_io_common+double.cxx
  minus/minus/Templates/minus_util+double.cxx
)
# this library is used since we want to precompile instantiations
# but header-only will also work, see README.md
add_library( minus ${minus_sources} )
# always build minus with specially tweaked optimizations.
# these were copied directly from minus original flags
set_target_properties(minus PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(minus PROPERTIES COMPILE_FLAGS "-O3 -march=native -ffast-math")
set_target_properties(minus PROPERTIES COMPILE_DEFINITIONS "NDEBUG")
find_package(Threads REQUIRED)
target_include_directories(minus PUBLIC ${MINUS_INCLUDE_DIRS})
target_link_libraries(minus PUBLIC Threads::Threads -pthread)
