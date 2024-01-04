
// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SCENE_INITIALIZER_STELLAR_HPP
#define OPENMVG_SFM_SCENE_INITIALIZER_STELLAR_HPP

#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializer.hpp"

namespace openMVG {
namespace sfm {

// Initialize a sfm_data with a "largest" stellar reconstruction:
// - a set of N cameras that share a common node.
class SfMSceneInitializerStellar : public SfMSceneInitializer
{
public:
  SfMSceneInitializerStellar(
    SfM_Data & sfm_data,
    const Features_Provider * features_provider,
    const Matches_Provider * matches_provider);

  ~SfMSceneInitializerStellar() override = default;

  bool Process() override;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SCENE_INITIALIZER_STELLAR_HPP
