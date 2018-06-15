// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_SfM.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/hypercluster.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include <iomanip>


#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {


SubmapSfMReconstructionEngine::SubmapSfMReconstructionEngine(const HsfmSubmap & submap, const tracks::STLMAPTracks & map_tracks, const std::string & soutDirectory, const std::string & sloggingFile)
  : SequentialSfMReconstructionEngine(submap.sfm_data, map_tracks, soutDirectory, sloggingFile)
{}

// we override this method with an empty function since it is
// the only thing that we don't use from SequentialSfMReconstructionEngine::Process.
// In this class, map tracks is defined at construction
bool SubmapSfMReconstructionEngine::InitLandmarkTracks()
{
  // Initialize the shared track visibility helper
  shared_track_visibility_helper_.reset(new openMVG::tracks::SharedTrackVisibilityHelper(map_tracks_));

  return true;
}

} // namespace sfm
} // namespace openMVG
