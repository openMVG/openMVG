// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_VERSION_HPP
#define OPENMVG_VERSION_HPP

#define OPENMVG_VERSION_MAJOR 1
#define OPENMVG_VERSION_MINOR 4
#define OPENMVG_VERSION_REVISION 0

// Preprocessor to string conversion
#define OPENMVG_TO_STRING_HELPER(x) #x
#define OPENMVG_TO_STRING(x) OPENMVG_TO_STRING_HELPER(x)

// OpenMVG version as a string; for example "1.4.0".
#define OPENMVG_VERSION_STRING OPENMVG_TO_STRING(OPENMVG_VERSION_MAJOR) "." \
  OPENMVG_TO_STRING(OPENMVG_VERSION_MINOR) "." \
  OPENMVG_TO_STRING(OPENMVG_VERSION_REVISION)

#endif  // OPENMVG_VERSION_HPP
