// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_STELLAR_ENGINE_HPP
#define OPENMVG_SFM_STELLAR_ENGINE_HPP

#include <memory>
#include <string>

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/sfm_engine.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

enum class EGraphSimplification : uint8_t {
  NONE,
  STAR_X,
  MST_X
};

/// A Stellar SfM Pipeline Reconstruction Engine:
/// - Use a local relative motion refinement,
/// - Fuse groups or relative motions.
class StellarSfMReconstructionEngine : public ReconstructionEngine
{
public:

  StellarSfMReconstructionEngine(
    const SfM_Data & sfm_data,
    const std::string & out_directory_logging,
    const std::string & logging_file = "");

  ~StellarSfMReconstructionEngine() override;

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);
  void SetGraphSimplification(const EGraphSimplification graph_simplification, const int value);

  bool Process() override;

protected:

  /// Compute relative rotations & translations
  void ComputeRelativeMotions
  (
    Hash_Map<Pair, geometry::Pose3> & relative_poses
  ) const;

  // Structure to store a stellar reconstruction
  // A series of relative motions centered in one pose Id
  struct StellarPodRelativeMotions
  {
    struct RelativeMotion
    {
      Mat3 rotation;
      Vec3 translation;
    };
    Hash_Map<Pair, RelativeMotion> relative_motions; // relative motions in the stellar pod
  };

  bool ComputeStellarReconstructions
  (
    const Hash_Map<Pair, geometry::Pose3> & relative_poses,
    Hash_Map<IndexT, StellarPodRelativeMotions> & stellar_reconstruction_per_pose
  ) const;

  /// Compute the global rotations from the relative rotations
  bool Compute_Global_Rotations
  (
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    Hash_Map<IndexT, Mat3> & map_globalR
  ) const;

  /// Compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    const Hash_Map<IndexT, StellarPodRelativeMotions> & stellar_reconstruction_per_pose
  );

  bool Compute_Initial_Structure(const int min_covisibility = 2);

  bool Adjust();

private:

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string logging_file_;

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;

  //-- Graph simplication for Stellar configurations
  EGraphSimplification graph_simplification_;
  int graph_simplification_value_; // Additional parameters (nb_tree or star satellite, ...)
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_STELLAR_ENGINE_HPP
