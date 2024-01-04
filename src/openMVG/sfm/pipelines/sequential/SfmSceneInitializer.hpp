
// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SCENE_INITIALIZER_HPP
#define OPENMVG_SFM_SCENE_INITIALIZER_HPP

namespace openMVG {
namespace sfm {

struct SfM_Data;
struct Features_Provider;
struct Matches_Provider;

// Interface to initialize a sfm_data with some camera poses.
// It allows to generate a 3D seed for a sequential SfM pipeline.
class SfMSceneInitializer
{
public:

  SfMSceneInitializer(
    SfM_Data & sfm_data,
    const Features_Provider * features_provider,
    const Matches_Provider * matches_provider)
    :sfm_data_(sfm_data),
     features_provider_(features_provider),
     matches_provider_(matches_provider)
  {
  };

  virtual ~SfMSceneInitializer() = default;

  virtual bool Process() {
    return true;
  };

  const SfM_Data & Get_sfm_data() const {return sfm_data_;};

protected:
  SfM_Data & sfm_data_;
  const Features_Provider * features_provider_;
  const Matches_Provider * matches_provider_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SCENE_INITIALIZER_HPP
