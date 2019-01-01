// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_LOCALIZATION_SEQUENTIAL2_SFM_HPP
#define OPENMVG_SFM_LOCALIZATION_SEQUENTIAL2_SFM_HPP

#include <set>
#include <string>
#include <vector>

#include "openMVG/sfm/pipelines/sfm_engine.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/tracks/tracks.hpp"

namespace htmlDocument { class htmlDocumentStream; }

namespace openMVG {
namespace sfm {

struct Features_Provider;
struct Matches_Provider;
class SfMSceneInitializer;

/// Sequential SfM Pipeline Reconstruction Engine.
/// This engine uses existing poses or starts from scratch the reconstruction
class SequentialSfMReconstructionEngine2 : public ReconstructionEngine
{
public:

  SequentialSfMReconstructionEngine2(
    SfMSceneInitializer * scene_initializer,
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~SequentialSfMReconstructionEngine2() override;

  virtual bool Process() override;

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  /// Initialize tracks
  bool InitTracksAndLandmarks();

  /// Triangulate tracks
  bool Triangulation();

  /// Adding missing view (Try to find the pose of the missing camera)
  bool AddingMissingView(const float & track_inlier_ratio);

  /// Adjust intrinsics, landmark and extrinsics according the user config.
  bool BundleAdjustment();

  /**
   * Set the default lens distortion type to use if it is declared unknown
   * in the intrinsics camera parameters by the previous steps.
   *
   * It can be declared unknown if the type cannot be deduced from the metadata.
   */
  void SetUnknownCameraType(const cameras::EINTRINSIC camType)
  {
    cam_type_ = camType;
  }

private:

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  cameras::EINTRINSIC cam_type_; // The camera type for the unknown cameras

  //-- Data provider
  Features_Provider * features_provider_;
  Matches_Provider  * matches_provider_;

  //-- Reconstruction Initialization
  SfMSceneInitializer * scene_initializer_;

  /// Putative landmark with view id visibility
  Landmarks landmarks_;
  /// Tracking (used to build landmark visibility and compute 2D-3D visibility)
  openMVG::tracks::STLMAPTracks map_tracks_;
  /// Helper to compute fast 2D-3D visibility
  std::unique_ptr<openMVG::tracks::SharedTrackVisibilityHelper> shared_track_visibility_helper_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_HPP
