
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP

#include "openMVG/sfm/pipelines/sfm_engine.hpp"

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"

namespace openMVG{
namespace sfm{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_RelativeMotions : public ReconstructionEngine
{
public:

  GlobalSfMReconstructionEngine_RelativeMotions(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_RelativeMotions();

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  void SetRotationAveragingMethod(ERotationAveragingMethod eRotationAveragingMethod);
  void SetTranslationAveragingMethod(ETranslationAveragingMethod eTranslation_averaging_method_);

  virtual bool Process();

protected:
  /// Compute from relative rotations the global rotations of the camera poses
  bool Compute_Global_Rotations
  (
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    Hash_Map<IndexT, Mat3> & map_globalR
  );

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );

  /// Compute the initial structure of the scene
  bool Compute_Initial_Structure
  (
    matching::PairWiseMatches & tripletWise_matches
  );

  // Adjust the scene (& remove outliers)
  bool Adjust();

private:
  /// Compute relative rotations
  void Compute_Relative_Rotations
  (
    openMVG::rotation_averaging::RelativeRotations & vec_relatives_R
  );

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  ERotationAveragingMethod eRotation_averaging_method_;
  ETranslationAveragingMethod eTranslation_averaging_method_;

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
