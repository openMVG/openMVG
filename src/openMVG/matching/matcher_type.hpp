
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace openMVG{
namespace matching{

enum EMatcherType
{
  BRUTE_FORCE_L2,
  ANN_L2,
  CASCADE_HASHING_L2,
  BRUTE_FORCE_HAMMING
};

} // namespace matching
} // namespace openMVG
