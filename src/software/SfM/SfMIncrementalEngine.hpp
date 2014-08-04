
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_INCREMENTAL_ENGINE_H
#define OPENMVG_SFM_INCREMENTAL_ENGINE_H

#include "openMVG/numeric/numeric.h"

#include "software/SfM/SfMEngine.hpp"
#include "software/SfM/SfMReconstructionData.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "openMVG/features/features.hpp"
#include <openMVG/tracks/tracks.hpp>

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
  void ColorizeTracks(std::vector<Vec3> & vec_tracksColor) const;

  const reconstructorHelper & refToReconstructorHelper() const
  {
    return _reconstructorData;
  }

  /// Bundle adjustment to refine Structure and Motion
  void BundleAdjustment();

  // Return MSE (Mean Square Error) and an histogram of residual values.
  double ComputeResidualsHistogram(Histogram<double> * histo);

  const std::vector<std::string> getFilenamesVector() const
  {
    std::vector<std::string> vec_fileNames;
    for ( std::vector<openMVG::SfMIO::CameraInfo>::const_iterator iter_camInfo = _vec_camImageNames.begin();
      iter_camInfo != _vec_camImageNames.end();
      iter_camInfo++ )
    {
      vec_fileNames.push_back( iter_camInfo->m_sImageName);
    }
    return vec_fileNames;
  }

  const openMVG::tracks::STLMAPTracks & getTracks() const
    {return _map_reconstructed;}

  const std::vector< std::pair<size_t, size_t> > getImagesSize() const
  {
    std::vector< std::pair<size_t, size_t> > vec_imageSize;
    for ( std::vector<openMVG::SfMIO::CameraInfo>::const_iterator iter_camInfo = _vec_camImageNames.begin();
      iter_camInfo != _vec_camImageNames.end();
      iter_camInfo++ )
    {
      std::vector<openMVG::SfMIO::IntrinsicCameraInfo>::const_iterator it_intrinsic = _vec_intrinsicGroups.begin();
      std::advance(it_intrinsic, iter_camInfo->m_intrinsicId);
      vec_imageSize.push_back( std::make_pair( it_intrinsic->m_w, it_intrinsic->m_h ) );
    }
    return vec_imageSize;
  }

  void setInitialPair(std::pair<size_t,size_t> initialPair)
  {
    _initialpair = initialPair;
  }

  void setIfUseBundleAdjustment(bool bUseBundleAdjustment)
  {
    _bUseBundleAdjustment = bUseBundleAdjustment;
  }
  
  void setIfRefinePrincipalPointAndRadialDisto(bool bRefinePPandDisto)
  {
    _bRefinePPandDisto = bRefinePPandDisto;
  }

  void setIfRefineFocal(bool bRefineFocal)
  {
    _bRefineFocal = bRefineFocal;
  }

private:

  // -----
  // Input data
  // ----
  std::vector<openMVG::SfMIO::CameraInfo> _vec_camImageNames;
  std::vector<openMVG::SfMIO::IntrinsicCameraInfo> _vec_intrinsicGroups;
  std::map< size_t, std::vector<SIOPointFeature> > _map_feats; // feature per images

  // Intrinsic Id per imageId
  std::map<size_t, size_t> _map_IntrinsicIdPerImageId;

  //-- Visibility information
  openMVG::matching::PairWiseMatches _map_Matches_F; // pairwise matches for Fundamental model
  openMVG::tracks::STLMAPTracks _map_tracks; // reconstructed track (visibility per 3D point)

  //-- configuration of the reconstruction
  std::pair<size_t,size_t> _initialpair;
  bool _bUseBundleAdjustment;
  bool _bRefinePPandDisto; // Boolean used to know if Principal point and Radial disto is refined
  bool _bRefineFocal;      // Boolean used to know if we refine focal or not

  // -----
  // Future reconstructed data
  // ----
  
  // helper to save reconstructed data (Camera and 3D points)
  reconstructorHelper _reconstructorData; 

  // tracks that are reconstructed during the sequential SfM process
  openMVG::tracks::STLMAPTracks _map_reconstructed;

  // Remaining image indexes that could be sequentially added
  std::set<size_t> _set_remainingImageId;

  //store in which order image have been added
  std::vector<size_t> _vec_added_order; 

  // Per camera confidence (A contrario estimated threshold error)
  std::map<size_t, double> _map_ACThreshold;

  /// List of images that belong to a common intrinsic group
  std::map<size_t, std::vector<size_t> > _map_ImagesIdPerIntrinsicGroup;
  std::map<size_t, Vec6 > _map_IntrinsicsPerGroup;

  // -----
  // Reporting ..
  // ----
  bool _bHtmlReport;
  std::auto_ptr<htmlDocument::htmlDocumentStream> _htmlDocStream;

};


} // namespace openMVG

#endif // OPENMVG_SFM_INCREMENTAL_ENGINE_H
