// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_HPP
#define OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_HPP

#include <set>
#include <string>
#include <vector>

#include "openMVG/sfm/pipelines/sfm_engine.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/tracks/tracks.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace { template <typename T> class Histogram; }

namespace openMVG {
namespace sfm {

struct Features_Provider;
struct Matches_Provider;

/// Sequential SfM Pipeline Reconstruction Engine.
class SequentialSfMReconstructionEngine : public ReconstructionEngine
{
public:

  SequentialSfMReconstructionEngine(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~SequentialSfMReconstructionEngine() override;

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  virtual bool Process() override;

  void setInitialPair(const Pair & initialPair)
  {
    initial_pair_ = initialPair;
  }

  /// Initialize tracks
  bool InitLandmarkTracks();

  /// Select a candidate initial pair
  bool ChooseInitialPair(Pair & initialPairIndex) const;

  /// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
  bool MakeInitialPair3D(const Pair & initialPair);

  /// Automatic initial pair selection (based on a 'baseline' computation score)
  bool AutomaticInitialPairChoice(Pair & initialPair) const;

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

protected:


private:

  /// Return MSE (Mean Square Error) and a histogram of residual values.
  double ComputeResidualsHistogram(Histogram<double> * histo);

  /// List the images that the greatest number of matches to the current 3D reconstruction.
  bool FindImagesWithPossibleResection(std::vector<uint32_t> & vec_possible_indexes);

  /// Add a single Image to the scene and triangulate new possible tracks.
  bool Resection(const uint32_t imageIndex);

  /// Bundle adjustment to refine Structure; Motion and Intrinsics
  bool BundleAdjustment();

  /// Discard track with too large residual error
  bool badTrackRejector(double dPrecision, size_t count = 0);

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  Pair initial_pair_;
  cameras::EINTRINSIC cam_type_; // The camera type for the unknown cameras

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;

  // Temporary data
  openMVG::tracks::STLMAPTracks map_tracks_; // putative landmark tracks (visibility per 3D point)

  // Helper to compute if some image have some track in common
  std::unique_ptr<openMVG::tracks::SharedTrackVisibilityHelper> shared_track_visibility_helper_;

  Hash_Map<IndexT, double> map_ACThreshold_; // Per camera confidence (A contrario estimated threshold error)

  std::set<uint32_t> set_remaining_view_id_;     // Remaining camera index that can be used for resection
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_HPP
