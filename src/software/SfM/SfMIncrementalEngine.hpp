
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_INCREMENTAL_ENGINE_H
#define OPENMVG_SFM_INCREMENTAL_ENGINE_H

#include "openMVG/numeric/numeric.h"

#include "software/SfM/SfMSimpleCamera.hpp"
#include "software/SfM/SfMEngine.hpp"
#include "software/SfM/SfMReconstructionData.hpp"
#include "openMVG/features/features.hpp"

#include "openMVG/tracks/tracks.hpp"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::tracks;

#include <memory>

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

namespace openMVG{

//Estimate E -> So R,t for the first pair
// Maintain a track list that explain 3D reconstructed scene
// Add images with Resection with the 3D tracks.
class IncrementalReconstructionEngine : public ReconstructionEngine
{
public:
  IncrementalReconstructionEngine(const std::string & sImagePath,
    const std::string & sMatchesPath, const std::string & sOutDirectory,
    bool bHtmlReport = false);

  ~IncrementalReconstructionEngine();

  virtual bool Process();

private:
  /// Read input data (point correspondences, K matrix)
  bool ReadInputData();

  /// Find the best initial pair
  bool InitialPairChoice( std::pair<size_t,size_t> & initialPairIndex);

  /// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
  bool MakeInitialPair3D(const std::pair<size_t,size_t> & initialPair);

  /// List the images that the greatest number of matches to the current 3D reconstruction.
  bool FindImagesWithPossibleResection(std::vector<size_t> & vec_possible_indexes);

  /// Add to the current scene the desired image indexes.
  bool Resection(std::vector<size_t> & vec_possible_indexes);

  /// Add a single Image to the scene and triangulate new possible tracks
  bool Resection(size_t imageIndex);

  /// Discard track with too large residual error
  size_t badTrackRejector(double dPrecision);

public:
  /// Give a color to all the 3D points
  void ColorizeTracks(std::vector<Vec3> & vec_color);

  const reconstructorHelper & refToReconstructorHelper() const
  {
    return _reconstructorData;
  }

  /// Bundle adjustment Refine Structure and Motion or Structure only
  void BundleAdjustment(bool bStructureAndMotion = true);

  // Return MSE (Mean Square Error) and an histogram of residual values.
  double ComputeResidualsHistogram(Histogram<double> * histo);

  const std::vector<std::string> & getFilenamesVector() const
    {return _vec_fileNames;}

  const openMVG::tracks::STLMAPTracks & getTracks() const
    {return _map_reconstructed;}

  const std::vector< std::pair<size_t, size_t> > & getImagesSize() const
    {return _vec_imageSize;}

private:

  // -----
  // Input data
  // ----

  std::vector<std::string> _vec_fileNames; // Images considered for the reconstruction
  std::map< size_t, std::vector<SIOPointFeature> > _map_feats; // feature per images

  std::vector< std::pair<size_t, size_t> > _vec_imageSize; // Size of each image
  Mat3 _K; // Intrinsic guess for the used camera.

  typedef tracks::mapPairWiseMatches STLPairWiseMatches;
  STLPairWiseMatches _map_Matches_H; // pairwise matches for Homography model
  STLPairWiseMatches _map_Matches_F; // pairwise matches for Fundamental model

  openMVG::tracks::STLMAPTracks _map_tracks; // reconstructed track (visibility per 3D point)

  // -----
  // Future reconstructed data
  // ----
  reconstructorHelper _reconstructorData; // helper to save reconstructed data (Camera and 3D points)

  openMVG::tracks::STLMAPTracks _map_reconstructed;

  std::set<size_t> _set_remainingImageId;  // Remaining image index that could be incrementally added

  std::vector<size_t> _vec_added_order; //show the added image order of the image

  // -----
  // Reporting ..
  // ----
  bool _bHtmlReport;
  std::auto_ptr<htmlDocument::htmlDocumentStream> _htmlDocStream;

};


} // namespace openMVG

#endif // OPENMVG_SFM_INCREMENTAL_ENGINE_H
