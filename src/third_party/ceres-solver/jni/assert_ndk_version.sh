#!/bin/bash

# Bash script to assert that the current version of the NDK is at least the
# specified version. Prints 'true' to standard out if it's the right version,
# 'false' if it's not.
#
# Typically used like this, in your jni/Android.mk:
#
#   ifneq ($(shell $(LOCAL_PATH)/assert_ndk_version.sh "r5c" "ndk-dir"), true)
#     $(error NDK version r5c or greater required)
#   endif
#
# See https://gist.github.com/2878774 for asserting SDK version.
#
# Retrieved from: https://gist.github.com/jorgenpt/1961404 on 2014-06-03.
#
# Copyright (c) 2012, Lookout, Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Author: jorgenpt@gmail.com (Jorgen Tjerno)
#         alexs.mac@gmail.com (Alex Stewart)

# Extracts 'r5c' into '5 c', also handles newer versions of the form
# 'r9d (64-bit)' and versions >= 10.
function get_major_minor() {
  # r9d (64-bit) -> '9d', also handle versions >= 10.
  local version=$(echo "$1" | sed 's/r\([0-9]\{1,2\}[a-z]\{0,1\}\).*/\1/')
  local major=$(echo "$version" | sed 's/\([0-9]\{1,2\}\).*/\1/')
  local minor=$(echo "$version" | sed 's/^[0-9]*//')
  echo "$major $minor"
}

if [[ -z "$2" ]]; then
  echo "Usage: $0 <required version> <NDK_ROOT>" >&2
  echo " For example: $0 r5c android-ndk-r9d" >&2
  exit 1
fi

# Assert that the expected version is at least 4.
declare -a expected_version
expected_version=( $(get_major_minor "$1") )
if [[ ${expected_version[0]} -le 4 ]]; then
  echo "Cannot test for versions less than r5: r4 doesn't have a version file." >&2
  echo false
  exit 1
fi

release_file="$2/RELEASE.TXT"

# NDK version r4 or earlier doesn't have a RELEASE.txt, and we just asserted
# that the person was looking for r5 or above, so that implies that this is an
# invalid version.
if [ ! -s "$release_file" ]; then
  echo false
  exit 0
fi

# Make sure the data is at least kinda sane.
version=$(grep '^r' $release_file)
declare -a actual_version
actual_version=( $(get_major_minor "$version") )
if [ -z "$version" ] || [ -z "${actual_version[0]}" ]; then
  echo "Invalid RELEASE.txt: $(cat $release_file)" >&2
  echo false
  exit 1
fi

if [[ ${actual_version[0]} -lt ${expected_version[0]} ]]; then
  echo "false"
elif [[ ${actual_version[0]} -eq ${expected_version[0]} ]]; then
  # This uses < and not -lt because they're string identifiers (a, b, c, etc)
  if [[ "${actual_version[1]}" < "${expected_version[1]}" ]]; then
    echo "false"
  else
    echo "true"
  fi
else
  echo "true"
fi
