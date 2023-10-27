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

#include "openMVG/sfm/pipelines/sequential/sequential_SfM_base.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/multiview/solver_resection.hpp"
#include "openMVG/multiview/triangulation_method.hpp"
#include "openMVG/tracks/tracks.hpp"

namespace htmlDocument { class htmlDocumentStream; }

namespace openMVG {
namespace sfm {

struct Features_Provider;
struct Matches_Provider;

/// Sequential SfM Pipeline Reconstruction Engine.
class SequentialSfMReconstructionEngine : public SequentialSfMReconstructionEngineBase
{
public:

  SequentialSfMReconstructionEngine(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "")
    :
    SequentialSfMReconstructionEngineBase(sfm_data, soutDirectory, loggingFile)
  { }

  ~SequentialSfMReconstructionEngine() override
  { };

  virtual bool Process() override;

  // tests that in fact the distortion is not identity
  static bool isDistortionZero(const cameras::IntrinsicBase *cam) 
  {
    Vec2 v(2,3), v_ud = cam->get_ud_pixel(v); // dummy
    return (v-v_ud).norm() < 1e-8;
  }


protected:


private:

  /// Highlevel methods -------------------------------------------------------
  /// Make initial 2- or 3-view reconstruction seed (robust plus BA and initial filters)
  bool MakeInitialSeedReconstruction();

  /// Compute the initial 3D seed (First camera: {R=Id|t=0}, 
  /// 2nd and 3rd estimated {R|t} by 3-point trifocal FABBRI CVPR20)
  bool MakeInitialTriplet3D(const Triplet & current_triplet);

  /// Incremental rec once a seed rec is established
  bool ResectOneByOneTilDone();

  /// Lowerlevel methods -------------------------------------------------------
  void ReconstructAllTangents();

  /// Add a single Image to the scene and triangulate new possible tracks.
  bool Resection(const uint32_t imageIndex);
  void ResectionAddTracks(IndexT I);

  /// See if all observations are adequately filled-in
  /// Test assumptions about the code, eg links in observation feature id,
  /// and actual features
  /// To be run after a major rec
  bool ConsistencyCheck() const;

  // consistency checks for oriented datastructures that run parallel to
  // the usual landmarks/observations 
  bool ConsistencyCheckOriented() const;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_LOCALIZATION_SEQUENTIAL_SFM_HPP
