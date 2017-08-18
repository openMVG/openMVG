// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/tracks/tracks.hpp"

namespace openMVG {
namespace sfm {

struct HsfmSubmap;

class SubmapSfMReconstructionEngine : public SequentialSfMReconstructionEngine
{
public:
  SubmapSfMReconstructionEngine(
      const HsfmSubmap & submap,
      const tracks::STLMAPTracks & map_tracks,
      const std::string & soutDirectory,
      const std::string & sloggingFile);

  ~SubmapSfMReconstructionEngine() = default;

  bool InitLandmarkTracks();

private:
};

}
}
