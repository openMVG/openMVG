// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/sfm.hpp"

// Basic storage of OpenMVG data related to a scene
struct Document
{
  openMVG::sfm::SfM_Data sfm_data;
  std::shared_ptr<openMVG::sfm::Features_Provider> feats_provider;
  std::shared_ptr<openMVG::sfm::Matches_Provider> matches_provider;
};
