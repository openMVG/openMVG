#!/bin/sh
#
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
# Author: keir@google.com
#         settinger@google.com
#
# Ceres build script for Android. To build the ceres.so libraray for Android,
# cd into an empty directory and run the script. Usage:
#
#   build_android.sh \
#   <path to Android NDK> \
#   <path to Eigen> \
#   <path to Ceres source>
#
#   make
#
# This script exists only to make it easier to get Ceres building on Android;
# as one can see from the code below, it is only a matter of extracting a
# standalone NDK toolchain from the NDK, and getting the right arguments to
# CMake to get it to work.
#
# Android NDK version r5 or higher is required. Jellybean does not work out of
# the box, since the android-cmake toolchain is not yet updated to for it.
#
# Note: You may wish to run 'ccmake .', the CMake curses GUI, in order to tweak
# the build parameters that are set by default. There are a few settings to
# consider changing:
#
# SCHUR_SPECIALIZATIONS:
#
# Consider if enabling the schur specializations is a big enough win for the
# problem you are solving, since compiling the schur specializations
# considerably increases the binary size. Disable them by running 'ccmake .',
# finding the SCHUR_SPECIALIZATIONS variable, hitting enter (toggling to "off"),
# pressing 'c' to generate, then 'g' to generate and exit, followed by 'make'.
#
# EXECUTABLE_OUTPUT_PATH
# LIBRARY_OUTPUT_PATH
# LIBRARY_OUTPUT_PATH_ROOT:
#
# In normal CMake builds, where you do an out of source build, the source
# directory is untouched when building. However, by default the Android CMake
# toolchain selects locations under your *source* tree for the final library
# and binary destinations. For example, if your Ceres git tree is under
# ceres-solver.git/ and the build directory you are using is
# ceres-android-bin/, the resulting binaries will live in ceres-solver.git/
# (but not the intermediate .o files!) By changing the variables
# EXECUTABLE_OUTPUT_PATH, LIBRARY_OUTPUT_PATH, and LIBRARY_OUTPUT_PATH_ROOT to
# something under e.g. ceres-android-bin/ then true out-of-ource builds work.

if [ $# -ne 3 ] ; then
  echo "usage: build_android.sh \\"
  echo "       <path to Android NDK> \\"
  echo "       <path to Eigen> \\"
  echo "       <path to Ceres source>"
  exit 1
fi

if [ -f "CMakeLists.txt" ] ; then
  echo "ERROR: Can't run from inside the source tree."
  echo "       Make a new directory that's not under"
  echo "       the main Ceres source tree."
  exit 1
fi

# For some reason, using the NDK in-place doesn't work even though the
# android-cmake toolchain appears to support it.
#
# TODO(keir): Figure out the issue with the stand-alone NDK and don't create
# the extra stand-alone toolchain. Also test with other NDK versions and add
# explicit checks to ensure a compatible version is used.
ANDROID_NDK=$1
MAKE_ANDROID_TOOLCHAIN=$ANDROID_NDK/build/tools/make-standalone-toolchain.sh
if [ ! -f $MAKE_ANDROID_TOOLCHAIN ] ; then
  echo "ERROR: First argument doesn't appear to be the NDK root; missing:"
  echo "       $MAKE_ANDROID_TOOLCHAIN"
  exit 1
fi

EIGEN_DIR=$2
if [ ! -f $EIGEN_DIR/eigen3.pc.in ] ; then
  echo "ERROR: Second argument doesn't appear to be Eigen3; missing:"
  echo "       $EIGEN_DIR/eigen3.pc.in"
  exit 1
fi

CERES_SOURCE_ROOT=$3
if [ ! -f "$CERES_SOURCE_ROOT/internal/ceres/CMakeLists.txt" ] ; then
  echo "ERROR: Third argument doesn't appear to be the Ceres source directory."
  exit 1
fi
echo "Using Ceres source directory: $CERES_SOURCE_ROOT"

# Make a standalone Android NDK toolchain if needed.
export ANDROID_STANDALONE_TOOLCHAIN="`pwd`/toolchain"
if [ ! -d "$ANDROID_STANDALONE_TOOLCHAIN" ] ; then
  echo "Extracting the Android GCC standalone toolchain to:"
  echo "    $ANDROID_STANDALONE_TOOLCHAIN..."
  $ANDROID_NDK/build/tools/make-standalone-toolchain.sh \
  --platform=android-8 \
  --install-dir=$ANDROID_STANDALONE_TOOLCHAIN
else
  echo "Found NDK standalone toolchain; skipping creation."
fi

# Get the Android CMake NDK toolchain file if needed.
if [ ! -d "android-cmake" ] ; then
  hg clone https://code.google.com/p/android-cmake/
else
  echo "Found Android-CMake toolchain; skipping download."
fi

ANDROID_CMAKE_TOOLCHAIN=android-cmake/toolchain/android.toolchain.cmake
if [ ! -f $ANDROID_CMAKE_TOOLCHAIN ] ; then
  echo "ERROR: It seems the toolchain file is missing:"
  echo "       $ANDROID_CMAKE_TOOLCHAIN"
  exit 1
fi

cmake $CERES_SOURCE_ROOT \
      -DCMAKE_TOOLCHAIN_FILE=$ANDROID_CMAKE_TOOLCHAIN \
      -DCMAKE_BUILD_TYPE=Release \
      -DEIGEN_INCLUDE=$EIGEN_DIR \
      -DBUILD_ANDROID=ON \
      -DSUITESPARSE=OFF \
      -DGFLAGS=OFF \
      -DCXSPARSE=OFF
