
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GLOBAL_SFM_ENGINE_H
#define OPENMVG_GLOBAL_SFM_ENGINE_H

#include "openMVG/numeric/numeric.h"

#include "openMVG/cameras/PinholeCamera.hpp"
#include "software/SfM/SfMEngine.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "software/SfM/SfMReconstructionData.hpp"
class SIOPointFeature;

#include "openMVG/tracks/tracks.hpp"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::tracks;

#include "openMVG/graph/triplet_finder.hpp"
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTriplets.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"

#include <memory>

namespace openMVG{

//------------------
//-- Bibliography --
//------------------
//- [1] "Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion."
//- Authors: Pierre MOULON, Pascal MONASSE and Renaud MARLET.
//- Date: December 2013.
//- Conference: ICCV.
//- [2] "Disambiguating visual relations using loop constraints."
//- Autors: Christopher ZACH, Manfred KLOPSCHITZ and Marc POLLEFEYS.
//- Date: 2010
//- Conference : CVPR.
//- [3] "Efficient and Robust Large-Scale Rotation Averaging"
//- Authors: Avishek Chatterjee and Venu Madhav Govindu
//- Date: December 2013.
//- Conference: ICCV.

// Implementation of [1]
// Some points differs from the [1] paper to ease open source port:
//-- Relative rotation inference:
//   - only the triplet rejection is performed (in [1] a Bayesian inference on cycle error is performed [2])
//-- Global rotation computation:
//    - a sparse least square,
//    - or, the L1 averaging method of [3].
//-- Linear Programming solver:
//   - in order to have the best performance it is advised to used the MOSEK LP backend.

enum ERotationAveragingMethod
{
  ROTATION_AVERAGING_L1 = 1,
  ROTATION_AVERAGING_L2 = 2
};

class GlobalReconstructionEngine : public ReconstructionEngine
{
public:
  GlobalReconstructionEngine(const std::string & sImagePath,
    const std::string & sMatchesPath,
    const std::string & sOutDirectory,
    const ERotationAveragingMethod & eRotationAveragingMethod,
    bool bHtmlReport = false);

  ~GlobalReconstructionEngine();

  virtual bool Process();

  /// Give a color to all the 3D points
  void ColorizeTracks(
    const STLMAPTracks & map_tracks, // tracks to be colorized
    std::vector<Vec3> & vec_tracksColor // output associated color
    ) const;

  //--
  // Accessors
  //--

  const reconstructorHelper & refToReconstructorHelper() const
  { return _reconstructorData;  }

  const openMVG::tracks::STLMAPTracks & getTracks() const
  { return _map_selectedTracks; }

  const std::vector<std::string> getFilenamesVector() const
  { return _vec_fileNames;  }

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

  //--
  // TYPEDEF
  //--
  typedef std::map< std::pair<size_t, size_t>, std::pair<Mat3, Vec3> > Map_RelativeRT;
  typedef std::map<size_t, PinholeCamera > Map_Camera;

private:
  /// Read input data (point correspondences, K matrix)
  bool ReadInputData();

  bool CleanGraph();

  void ComputeRelativeRt(Map_RelativeRT & vec_relatives);

  // Detect and remove the outlier relative rotations
  void rotationInference(Map_RelativeRT & map_relatives);

  // Compute the global rotations from relative rotations
  bool computeGlobalRotations(
    ERotationAveragingMethod eRotationAveragingMethod,
    const std::map<size_t, size_t> & map_cameraNodeToCameraIndex,
    const std::map<size_t, size_t> & map_cameraIndexTocameraNode,
    const Map_RelativeRT & map_relatives,
    std::map<size_t, Mat3> & map_globalR) const;

  // List the triplet of the image connection graph (_map_Matches_F)
  void tripletListing(std::vector< graphUtils::Triplet > & vec_triplets) const;

  // Relative rotations inference on relative rotations composition error along 3 length cycles (triplets).
  void tripletRotationRejection(
    std::vector< graphUtils::Triplet > & vec_triplets,
    Map_RelativeRT & map_relatives);

  // Compute relative translations over the graph of putative triplets
  void computePutativeTranslation_EdgesCoverage(
    const std::map<std::size_t, Mat3> & map_globalR,
    const std::vector< graphUtils::Triplet > & vec_triplets,
    std::vector<openMVG::lInfinityCV::relativeInfo > & vec_initialEstimates,
    matching::PairWiseMatches & newpairMatches) const;

  // Bundle adjustment : refine structure Xis and camera parameters (with optional refined parameters)
  void bundleAdjustment(
    Map_Camera & map_camera,
    std::vector<Vec3> & vec_allScenes,
    const STLMAPTracks & map_tracksSelected,
    bool bRefineRotation = true,
    bool bRefineTranslation = true,
    bool bRefineIntrinsics = false);

private:

  // -----
  // Input data
  // ----

  // Image considered for the reconstruction
  std::vector<std::string> _vec_fileNames;
  std::vector<openMVG::SfMIO::CameraInfo> _vec_camImageNames;
  std::vector<openMVG::SfMIO::IntrinsicCameraInfo> _vec_intrinsicGroups;
  std::map< size_t, std::vector<SIOPointFeature> > _map_feats; // feature per images

  matching::PairWiseMatches _map_Matches_F; // pairwise matches for Essential matrix model

  // Parameter
  ERotationAveragingMethod _eRotationAveragingMethod;

  //------
  //-- Mapping between camera node Ids and cameraIndex:
  //--------------
  // Graph node Ids mapping to camera Ids
  std::map<size_t, size_t> map_cameraNodeToCameraIndex; // graph node Id to 0->Ncam
  // camera Ids to graph node Ids
  std::map<size_t, size_t> map_cameraIndexTocameraNode; // 0->Ncam correspondance to graph node Id
  //--
  //----

  //-----
  //-- Reconstruction data
  //-----
  // Cameras (Motion)
  Map_Camera _map_camera;
  // Structure
  std::vector<Vec3> _vec_allScenes;
  // Structure visibility
  STLMAPTracks _map_selectedTracks; // reconstructed track (visibility per 3D point)
  // Scene and structure container (for disk output)
  reconstructorHelper _reconstructorData;
  //-----


  // -----
  // Reporting ..
  // ----
  bool _bHtmlReport;
  std::auto_ptr<htmlDocument::htmlDocumentStream> _htmlDocStream;

};


} // namespace openMVG

#endif // OPENMVG_GLOBAL_SFM_ENGINE_H
