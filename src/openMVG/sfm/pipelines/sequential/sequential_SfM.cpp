
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {

using namespace openMVG::geometry;
using namespace openMVG::cameras;

SequentialSfMReconstructionEngine::SequentialSfMReconstructionEngine(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory),
    _sLoggingFile(sloggingFile),
    _initialpair(Pair(0,0)),
    _camType(EINTRINSIC(PINHOLE_CAMERA_RADIAL3))
{
  if (!_sLoggingFile.empty())
  {
    // setup HTML logger
    _htmlDocStream = std::make_shared<htmlDocument::htmlDocumentStream>("SequentialReconstructionEngine SFM report.");
    _htmlDocStream->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("SequentialSfMReconstructionEngine")));
    _htmlDocStream->pushInfo("<hr>");

    _htmlDocStream->pushInfo( "Dataset info:");
    _htmlDocStream->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }
  // Init remaining image list
  for (Views::const_iterator itV = sfm_data.GetViews().begin();
    itV != sfm_data.GetViews().end(); ++itV)
  {
    _set_remainingViewId.insert(itV->second.get()->id_view);
  }
}

SequentialSfMReconstructionEngine::~SequentialSfMReconstructionEngine()
{
  if (!_sLoggingFile.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(_sLoggingFile.c_str());
    htmlFileStream << _htmlDocStream->getDoc();
  }
}

void SequentialSfMReconstructionEngine::SetFeaturesProvider(Features_Provider * provider)
{
  _features_provider = provider;
}

void SequentialSfMReconstructionEngine::SetMatchesProvider(Matches_Provider * provider)
{
  _matches_provider = provider;
}

bool SequentialSfMReconstructionEngine::Process() {

  //-------------------
  //-- Incremental reconstruction
  //-------------------

  if (!InitLandmarkTracks())
    return false;

  Pair initialPairIndex;
  if(! ChooseInitialPair(initialPairIndex))
    return false;

  // Initial pair Essential Matrix and [R|t] estimation.
  if(! MakeInitialPair3D(initialPairIndex))
    return false;

  size_t imageIndex = 0;
  size_t resectionGroupIndex = 0;
  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  std::vector<size_t> vec_possible_resection_indexes;
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    // std::cout << "Resection group start " << resectionGroupIndex << " with " << vec_possible_resection_indexes.size() << " images.\n";

    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (std::vector<size_t>::const_iterator iter = vec_possible_resection_indexes.begin();
      iter != vec_possible_resection_indexes.end(); ++iter)
    {
      const size_t currentIndex = imageIndex;
      ++imageIndex;
      const bool bResect = Resection(*iter);
      bImageAdded |= bResect;
      if (!bResect) {
        std::cerr << "\nResection of image: " << *iter << " was not possible" << std::endl;
      }
      _set_remainingViewId.erase(*iter);
    }

    if (bImageAdded)
    {
      // Scene logging as ply for visual debug
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
      Save(_sfm_data, stlplus::create_filespec(_sOutDirectory, os.str(), ".ply"), ESfM_Data(ALL));

      // std::cout << "Global Bundle start, resection group index: " << resectionGroupIndex << ".\n";
      int bundleAdjustmentIteration = 0;
      // Perform BA until all point are under the given precision
      do
      {
        // std::cout << "Resection group index: " << resectionGroupIndex << ", bundle iteration: " << bundleAdjustmentIteration << "\n";
        BundleAdjustment();
        ++bundleAdjustmentIteration;
      }
      while (badTrackRejector(4.0, 50) != 0);
    }
    ++resectionGroupIndex;
  }
  // Ensure there is no remaining outliers
  badTrackRejector(4.0, 0);

  //-- Reconstruction done.
  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << _sfm_data.GetPoses().size()
    << " from " << _sfm_data.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << _sfm_data.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  std::cout << "\nHistogram of residuals:" << h.ToString() << std::endl;

  if (!_sLoggingFile.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion process finished.";
    _htmlDocStream->pushInfo("<hr>");
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Structure from Motion (statistics):<br>"
      << "-- #Camera calibrated: " << _sfm_data.GetPoses().size()
      << " from " <<_sfm_data.GetViews().size() << " input images.<br>"
      << "-- #Tracks, #3D points: " << _sfm_data.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());

    _htmlDocStream->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

    const std::vector<double> xBin = h.GetXbinsValue();
    std::pair< std::pair<double,double>, std::pair<double,double> > range =
      autoJSXGraphViewport<double>(xBin, h.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("3DtoImageResiduals",600,300);
    jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    _htmlDocStream->pushInfo(jsxGraph.toStr());
  }
  return true;
}

/// Select a candidate initial pair
bool SequentialSfMReconstructionEngine::ChooseInitialPair(Pair & initialPairIndex) const
{
  if (_initialpair != Pair(0,0))
  {
    initialPairIndex = _initialpair;
  }
  else
  {

    // List Views that support valid intrinsic
    std::set<IndexT> valid_views;
    for (Views::const_iterator it = _sfm_data.GetViews().begin();
      it != _sfm_data.GetViews().end(); ++it)
    {
      const View * v = it->second.get();
      if( _sfm_data.GetIntrinsics().find(v->id_intrinsic) != _sfm_data.GetIntrinsics().end())
        valid_views.insert(v->id_view);
    }

    if (_sfm_data.GetIntrinsics().empty() || valid_views.empty())
    {
      std::cerr
        << "There is no defined intrinsic data in order to compute an essential matrix for the initial pair."
        << std::endl;
      return false;
    }

    std::cout << std::endl
      << "---------------------------------------------------\n"
      << "IncrementalReconstructionEngine::ChooseInitialPair\n"
      << "---------------------------------------------------\n"
      << " Pairs that have valid intrinsic and high support of points are displayed:\n"
      << " Choose one pair manually by typing the two integer indexes\n"
      << "---------------------------------------------------\n"
      << std::endl;

    // Try to list the 10 top pairs that have:
    //  - valid intrinsics,
    //  - valid estimated Fundamental matrix.
    std::vector< size_t > vec_NbMatchesPerPair;
    std::vector<openMVG::matching::PairWiseMatches::const_iterator> vec_MatchesIterator;
    const openMVG::matching::PairWiseMatches & map_Matches = _matches_provider->_pairWise_matches;
    for (openMVG::matching::PairWiseMatches::const_iterator
      iter = map_Matches.begin();
      iter != map_Matches.end(); ++iter)
    {
      const Pair current_pair = iter->first;
      if (valid_views.count(current_pair.first) &&
        valid_views.count(current_pair.second) )
      {
        vec_NbMatchesPerPair.push_back(iter->second.size());
        vec_MatchesIterator.push_back(iter);
      }
    }
    // sort the Pairs in descending order according their correspondences count
    using namespace stl::indexed_sort;
    std::vector< sort_index_packet_descend< size_t, size_t> > packet_vec(vec_NbMatchesPerPair.size());
    sort_index_helper(packet_vec, &vec_NbMatchesPerPair[0], std::min((size_t)10, vec_NbMatchesPerPair.size()));

    for (size_t i = 0; i < std::min((size_t)10, vec_NbMatchesPerPair.size()); ++i) {
      const size_t index = packet_vec[i].index;
      openMVG::matching::PairWiseMatches::const_iterator iter = vec_MatchesIterator[index];
      std::cout << "(" << iter->first.first << "," << iter->first.second <<")\t\t"
        << iter->second.size() << " matches" << std::endl;
    }

    // Ask the user to choose an initial pair (by set some view ids)
    std::cout << std::endl << " type INITIAL pair ids: X enter Y enter\n";
    int val, val2;
    if ( std::cin >> val && std::cin >> val2) {
      initialPairIndex.first = val;
      initialPairIndex.second = val2;
    }
  }

  std::cout << "\nPutative starting pair is: (" << initialPairIndex.first
      << "," << initialPairIndex.second << ")" << std::endl;

  // Check validity of the initial pair indices:
  if (_features_provider->feats_per_view.find(initialPairIndex.first) == _features_provider->feats_per_view.end() ||
      _features_provider->feats_per_view.find(initialPairIndex.second) == _features_provider->feats_per_view.end())
  {
    std::cerr << "At least one of the initial pair indices is invalid."
      << std::endl;
    return false;
  }
  return true;
}

bool SequentialSfMReconstructionEngine::InitLandmarkTracks()
{
  // Compute tracks from matches
  tracks::TracksBuilder tracksBuilder;

  {
    // List of features matches for each couple of images
    const openMVG::matching::PairWiseMatches & map_Matches = _matches_provider->_pairWise_matches;
    std::cout << std::endl << "Track building" << std::endl;

    tracksBuilder.Build(map_Matches);
    std::cout << std::endl << "Track filtering" << std::endl;
    tracksBuilder.Filter();
    std::cout << std::endl << "Track filtering : min occurence" << std::endl;
    tracksBuilder.FilterPairWiseMinimumMatches(20);
    std::cout << std::endl << "Track export to internal struct" << std::endl;
    //-- Build tracks with STL compliant type :
    tracksBuilder.ExportToSTL(_map_tracks);

    std::cout << std::endl << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats :
      //    - number of images
      //    - number of tracks
      std::set<size_t> set_imagesId;
      tracks::TracksUtilsMap::ImageIdInTracks(_map_tracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      tracks::TracksUtilsMap::TracksLength(_map_tracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)  {
        osTrack << "\t" << iter->first << "\t" << iter->second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }
  return _map_tracks.size() > 0;
}

/// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
bool SequentialSfMReconstructionEngine::MakeInitialPair3D(const Pair & current_pair)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  const size_t I = min(current_pair.first, current_pair.second);
  const size_t J = max(current_pair.first, current_pair.second);

  // a. Assert we have valid pinhole cameras
  const View * view_I = _sfm_data.GetViews().at(I).get();
  const Intrinsics::const_iterator iterIntrinsic_I = _sfm_data.GetIntrinsics().find(view_I->id_intrinsic);
  const View * view_J = _sfm_data.GetViews().at(J).get();
  const Intrinsics::const_iterator iterIntrinsic_J = _sfm_data.GetIntrinsics().find(view_J->id_intrinsic);

  if (iterIntrinsic_I == _sfm_data.GetIntrinsics().end() ||
      iterIntrinsic_J == _sfm_data.GetIntrinsics().end() )
  {
    return false;
  }

  const Pinhole_Intrinsic * cam_I = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_I->second.get());
  const Pinhole_Intrinsic * cam_J = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_J->second.get());
  if (cam_I == NULL || cam_J == NULL)
  {
    return false;
  }

  // b. Get common features between the two view
  // use the track to have a more dense match correspondence set
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  std::set<size_t> set_imageIndex;
  set_imageIndex.insert(I);
  set_imageIndex.insert(J);
  tracks::TracksUtilsMap::GetTracksInImages(set_imageIndex, _map_tracks, map_tracksCommon);

  //-- Copy point to arrays
  const size_t n = map_tracksCommon.size();
  Mat xI(2,n), xJ(2,n);
  size_t cptIndex = 0;
  for (openMVG::tracks::STLMAPTracks::const_iterator
    iterT = map_tracksCommon.begin(); iterT != map_tracksCommon.end();
    ++iterT, ++cptIndex)
  {
    tracks::submapTrack::const_iterator iter = iterT->second.begin();
    const size_t i = iter->second;
    const size_t j = (++iter)->second;

    Vec2 feat = _features_provider->feats_per_view[I][i].coords().cast<double>();
    xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
    feat = _features_provider->feats_per_view[J][j].coords().cast<double>();
    xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
  }

  // c. Robust estimation of the relative pose
  RelativePose_Info relativePose_info;

  const std::pair<size_t, size_t> imageSize_I(cam_I->w(), cam_I->h());
  const std::pair<size_t, size_t> imageSize_J(cam_J->w(), cam_J->h());

  if (!robustRelativePose(
    cam_I->K(), cam_J->K(), xI, xJ, relativePose_info, imageSize_I, imageSize_J, 4096))
  {
    std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
      << std::endl;
    return false;
  }
  std::cout << "A-Contrario initial pair residual: "
    << relativePose_info.found_residual_precision << std::endl;
  // Bound min precision at 1 pix.
  relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

  bool bRefine_using_BA = true;
  if (bRefine_using_BA)
  {
    // Refine the defined scene
    SfM_Data tiny_scene;
    tiny_scene.views.insert(*_sfm_data.GetViews().find(view_I->id_view));
    tiny_scene.views.insert(*_sfm_data.GetViews().find(view_J->id_view));
    tiny_scene.intrinsics.insert(*_sfm_data.GetIntrinsics().find(view_I->id_intrinsic));
    tiny_scene.intrinsics.insert(*_sfm_data.GetIntrinsics().find(view_J->id_intrinsic));

    // Init poses
    const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

    // Init structure
    const Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
    const Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
    Landmarks & landmarks = tiny_scene.structure;

    for (openMVG::tracks::STLMAPTracks::const_iterator
      iterT = map_tracksCommon.begin();
      iterT != map_tracksCommon.end();
      ++iterT)
    {
      // Get corresponding points
      tracks::submapTrack::const_iterator iter = iterT->second.begin();
      const size_t i = iter->second;
      const size_t j = (++iter)->second;

      const Vec2 x1_ = _features_provider->feats_per_view[I][i].coords().cast<double>();
      const Vec2 x2_ = _features_provider->feats_per_view[J][j].coords().cast<double>();

      Vec3 X;
      TriangulateDLT(P1, x1_, P2, x2_, &X);
      Observations obs;
      obs[view_I->id_view] = Observation(x1_, i);
      obs[view_J->id_view] = Observation(x2_, j);
      landmarks[iterT->first].obs = obs;
      landmarks[iterT->first].X = X;
    }
    Save(tiny_scene, stlplus::create_filespec(_sOutDirectory, "initialPair.ply"), ESfM_Data(ALL));

    // - refine only Structure and Rotations & translations (keep intrinsic constant)
    Bundle_Adjustment_Ceres::BA_options options(true, false);
    options._linear_solver_type = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.Adjust(tiny_scene, true, true, false))
    {
      return false;
    }

    // Save computed data
    const Pose3 pose_I = _sfm_data.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
    const Pose3 pose_J = _sfm_data.poses[view_J->id_pose] = tiny_scene.poses[view_J->id_pose];
    _map_ACThreshold.insert(std::make_pair(I, relativePose_info.found_residual_precision));
    _map_ACThreshold.insert(std::make_pair(J, relativePose_info.found_residual_precision));
    _set_remainingViewId.erase(view_I->id_view);
    _set_remainingViewId.erase(view_J->id_view);

    // List inliers and save them
    for (Landmarks::const_iterator iter = tiny_scene.GetLandmarks().begin();
      iter != tiny_scene.GetLandmarks().end(); ++iter)
    {
      const IndexT trackId = iter->first;
      const Landmark & landmark = iter->second;
      const Observations & obs = landmark.obs;
      Observations::const_iterator iterObs_xI = obs.begin();
      Observations::const_iterator iterObs_xJ = obs.begin();
      std::advance(iterObs_xJ, 1);

      const Observation & ob_xI = iterObs_xI->second;
      const IndexT & viewId_xI = iterObs_xI->first;

      const Observation & ob_xJ = iterObs_xJ->second;
      const IndexT & viewId_xJ = iterObs_xJ->first;

      const double angle = AngleBetweenRay(
        pose_I, cam_I, pose_J, cam_J, ob_xI.x, ob_xJ.x);
      const Vec2 residual_I = cam_I->residual(pose_I, landmark.X, ob_xI.x);
      const Vec2 residual_J = cam_J->residual(pose_J, landmark.X, ob_xJ.x);
      if ( angle > 2.0 &&
           pose_I.depth(landmark.X) > 0 &&
           pose_J.depth(landmark.X) > 0 &&
           residual_I.norm() < relativePose_info.found_residual_precision &&
           residual_J.norm() < relativePose_info.found_residual_precision)
      {
        _sfm_data.structure[trackId] = landmarks[trackId];
      }
      else  {
        // Remove this observation from the scene tracking data
        _map_tracks[trackId].erase(I);
        _map_tracks[trackId].erase(J);
        if (_map_tracks[trackId].size() < 2)  {
          _map_tracks.erase(trackId);
        }
      }
    }
    // Save outlier residual information
    Histogram<double> histoResiduals;
    std::cout << std::endl
      << "=========================\n"
      << " MSE Residual InitialPair Inlier: " << ComputeResidualsHistogram(&histoResiduals) << "\n"
      << "=========================" << std::endl;

    if (!_sLoggingFile.empty())
    {
      using namespace htmlDocument;
      _htmlDocStream->pushInfo(htmlMarkup("h1","Essential Matrix."));
      ostringstream os;
      os << std::endl
        << "-------------------------------" << "<br>"
        << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
        << view_I->s_Img_path << ","
        << view_J->s_Img_path << "<br>"
        << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
        << "-- Resection status: " << "OK" << "<br>"
        << "-- Nb points used for robust Essential matrix estimation: "
        << xI.cols() << "<br>"
        << "-- Nb points validated by robust estimation: "
        << _sfm_data.structure.size() << "<br>"
        << "-- % points validated: "
        << _sfm_data.structure.size()/static_cast<float>(xI.cols())
        << "<br>"
        << "-------------------------------" << "<br>";
      _htmlDocStream->pushInfo(os.str());

      _htmlDocStream->pushInfo(htmlMarkup("h2",
        "Residual of the robust estimation (Initial triangulation). Thresholded at: "
        + toString(relativePose_info.found_residual_precision)));

      _htmlDocStream->pushInfo(htmlMarkup("h2","Histogram of residuals"));

      std::vector<double> xBin = histoResiduals.GetXbinsValue();
      std::pair< std::pair<double,double>, std::pair<double,double> > range =
        autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

      htmlDocument::JSXGraphWrapper jsxGraph;
      jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
      jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
      jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
        relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
      jsxGraph.UnsuspendUpdate();
      jsxGraph.setViewport(range);
      jsxGraph.close();
      _htmlDocStream->pushInfo(jsxGraph.toStr());

      _htmlDocStream->pushInfo("<hr>");

      ofstream htmlFileStream( string(stlplus::folder_append_separator(_sOutDirectory) +
        "Reconstruction_Report.html").c_str());
      htmlFileStream << _htmlDocStream->getDoc();
    }
  }
  return !_sfm_data.structure.empty();
}

double SequentialSfMReconstructionEngine::ComputeResidualsHistogram(Histogram<double> * histo)
{
  // Collect residuals for each observation
  std::vector<float> vec_residuals;
  vec_residuals.reserve(_sfm_data.structure.size());
  for(Landmarks::const_iterator iterTracks = _sfm_data.GetLandmarks().begin();
      iterTracks != _sfm_data.GetLandmarks().end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for(Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = _sfm_data.GetViews().find(itObs->first)->second.get();
      const Pose3 pose = _sfm_data.GetPoseOrDie(view);
      const std::shared_ptr<IntrinsicBase> intrinsic = _sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      vec_residuals.push_back( fabs(residual(0)) );
      vec_residuals.push_back( fabs(residual(1)) );
    }
  }
  // Display statistics
  if (vec_residuals.size() > 1)
  {
    float dMin, dMax, dMean, dMedian;
    minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                            dMin, dMax, dMean, dMedian);
    if (histo)  {
      *histo = Histogram<double>(dMin, dMax, 10);
      histo->Add(vec_residuals.begin(), vec_residuals.end());
    }

    std::cout << std::endl << std::endl;
    std::cout << std::endl
      << "IncrementalReconstructionEngine::ComputeResidualsMSE." << "\n"
      << "\t-- #Tracks:\t" << _sfm_data.GetLandmarks().size() << std::endl
      << "\t-- Residual min:\t" << dMin << std::endl
      << "\t-- Residual median:\t" << dMedian << std::endl
      << "\t-- Residual max:\t "  << dMax << std::endl
      << "\t-- Residual mean:\t " << dMean << std::endl;

    return dMean;
  }
  return -1.0;
}

/// Functor to sort a vector of pair given the pair's second value
template<class T1, class T2, class Pred = std::less<T2> >
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left,
                    const std::pair<T1,T2>&right)
  {
    Pred p;
    return p(left.second, right.second);
  }
};

/**
 * @brief Estimate images on which we can compute the resectioning safely.
 *
 * @param[out] vec_possible_indexes: list of indexes we can use for resectioning.
 * @return False if there is no possible resection.
 *
 * Sort the images by the number of features already reconstructed.
 * Instead of returning the best one, we select a group of images, we can use
 * to do the resectioning in one step, which means without global Bundle Adjustment.
 * So we return all the images with at least 75% of the number of matches
 * of the best image.
 */
bool SequentialSfMReconstructionEngine::FindImagesWithPossibleResection(
  std::vector<size_t> & vec_possible_indexes)
{
  // Threshold used to select the best images
  static const double dThresholdGroup = 0.75;

  vec_possible_indexes.clear();

  if (_set_remainingViewId.empty() || _sfm_data.GetLandmarks().empty())
    return false;

  // Collect tracksIds
  std::set<size_t> reconstructed_trackId;
  std::transform(_sfm_data.GetLandmarks().begin(), _sfm_data.GetLandmarks().end(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  Pair_Vec vec_putative; // ImageId, NbPutativeCommonPoint
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (std::set<size_t>::const_iterator iter = _set_remainingViewId.begin();
        iter != _set_remainingViewId.end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      const size_t viewId = *iter;

      // Compute 2D - 3D possible content
      openMVG::tracks::STLMAPTracks map_tracksCommon;
      std::set<size_t> set_viewId;
      set_viewId.insert(viewId);
      tracks::TracksUtilsMap::GetTracksInImages(set_viewId, _map_tracks, map_tracksCommon);

      if (! map_tracksCommon.empty())
      {
        std::set<size_t> set_tracksIds;
        tracks::TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

        // Count the common possible putative point
        //  with the already 3D reconstructed trackId
        std::vector<size_t> vec_trackIdForResection;
        std::set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
          reconstructed_trackId.begin(),
          reconstructed_trackId.end(),
          std::back_inserter(vec_trackIdForResection));

#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          vec_putative.push_back( make_pair(viewId, vec_trackIdForResection.size()));
        }
      }
    }
  }

  if (vec_putative.empty())
    return false;

  // Sort by the number of matches into the 3D scene.
  std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<size_t, size_t, std::greater<size_t> >());

  // The number of matches of the best image.
  const IndexT M = vec_putative[0].second;
  vec_possible_indexes.push_back(vec_putative[0].first);

  // TEMPORARY
  // std::cout << std::endl << std::endl << "TEMPORARY return only the best image" << std::endl;
  // return true;
  // END TEMPORARY

  // Return all the images that have at least N% of the number of matches of the best image.
  const size_t threshold = static_cast<size_t>(dThresholdGroup * M);
  for (size_t i = 1; i < vec_putative.size() &&
    vec_putative[i].second > threshold; ++i)
  {
    vec_possible_indexes.push_back(vec_putative[i].first);
  }
  return true;
}

/**
 * @brief Add one image to the 3D reconstruction. To the resectioning of
 * the camera and triangulate all the new possible tracks.
 * @param[in] viewIndex: image index to add to the reconstruction.
 *
 * A. Compute 2D/3D matches
 * B. Look if intrinsic data is known or not
 * C. Do the resectioning: compute the camera pose.
 * D. Refine the pose of the found camera
 * E. Update the global scene with the new camera
 * F. Update the observations into the global scene structure
 * G. Triangulate new possible 2D tracks
 */
bool SequentialSfMReconstructionEngine::Resection(size_t viewIndex)
{
  using namespace tracks;

  std::cout << std::endl
    << "-------------------------------" << std::endl
    << "-- Resection of camera index: " << viewIndex << std::endl
    << "-------------------------------" << std::endl;

  // A. Compute 2D/3D matches
  // A1. list tracks ids used by the view
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  std::set<size_t> set_viewIndex;
  set_viewIndex.insert(viewIndex);
  TracksUtilsMap::GetTracksInImages(set_viewIndex, _map_tracks, map_tracksCommon);
  std::set<size_t> set_tracksIds;
  TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

  // A2. intersects the track list with the reconstructed
  std::set<size_t> reconstructed_trackId;
  std::transform(_sfm_data.GetLandmarks().begin(), _sfm_data.GetLandmarks().end(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  // trackIDs already reconstructed
  std::set<size_t> set_trackIdForResection;
  std::set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
    reconstructed_trackId.begin(),
    reconstructed_trackId.end(),
    std::inserter(set_trackIdForResection, set_trackIdForResection.begin()));

  if(set_trackIdForResection.empty())
  {
    // No match. The image has no connection with already reconstructed points.
    return false;
  }

  // Get back featId associated to a tracksID already reconstructed.
  // These 2D/3D associations will be used for the resection.
  std::vector<size_t> vec_featIdForResection;
  TracksUtilsMap::GetFeatIndexPerViewAndTrackId(map_tracksCommon,
    set_trackIdForResection,
    viewIndex,
    &vec_featIdForResection);

  // Create pt2D, and pt3D array
  Mat pt2D( 2, set_trackIdForResection.size());
  Mat pt3D( 3, set_trackIdForResection.size());

  // Get the view of the current camera to do the resectioning.
  const View * view_I = _sfm_data.GetViews().at(viewIndex).get();

  // B. Look if intrinsic data is known or not
  bool bKnownIntrinsic = true;
  Mat3 K = Mat3::Identity();
  const Intrinsics::const_iterator iterIntrinsic_I = _sfm_data.GetIntrinsics().find(view_I->id_intrinsic);
  Pinhole_Intrinsic * cam_I = NULL;
  if (iterIntrinsic_I == _sfm_data.GetIntrinsics().end())
  {
    bKnownIntrinsic = false;
  }
  else
  {
    cam_I = dynamic_cast<Pinhole_Intrinsic*>(iterIntrinsic_I->second.get());
    if (cam_I)
    {
      K = cam_I->K();
    }
    else
    {
      bKnownIntrinsic = false;
    }
  }

  size_t cpt = 0;
  std::set<size_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
  for (std::vector<size_t>::const_iterator iterfeatId = vec_featIdForResection.begin();
    iterfeatId != vec_featIdForResection.end();
    ++iterfeatId, ++iterTrackId, ++cpt)
  {
    pt3D.col(cpt) = _sfm_data.GetLandmarks().at(*iterTrackId).X;

    const Vec2 feat = _features_provider->feats_per_view[viewIndex][*iterfeatId].coords().cast<double>();
    if (bKnownIntrinsic)
      pt2D.col(cpt) = cam_I->get_ud_pixel(feat);
    else
      pt2D.col(cpt) = feat;
  }

  // C. Do the resectioning: compute the camera pose.
  std::vector<size_t> vec_inliers;
  Mat34 P;
  double errorMax = std::numeric_limits<double>::max();

  bool bResection = robustResection(
    std::make_pair( view_I->ui_width, view_I->ui_height ),
    pt2D, pt3D,
    &vec_inliers,
    // Use intrinsic guess if possible
    (bKnownIntrinsic) ? &K : NULL,
    &P, &errorMax);

  std::cout << std::endl
    << "-------------------------------" << std::endl
    << "-- Robust Resection of view: " << viewIndex << std::endl
    << "-- Resection status: " << bResection << std::endl
    << "-- #Points used for Resection: " << vec_featIdForResection.size() << std::endl
    << "-- #Points validated by robust Resection: " << vec_inliers.size() << std::endl
    << "-- Threshold: " << errorMax << std::endl
    << "-------------------------------" << std::endl;

  if (!_sLoggingFile.empty())
  {
    using namespace htmlDocument;
    ostringstream os;
    os << "Resection of Image index: <" << viewIndex << "> image: "
      << view_I->s_Img_path <<"<br> \n";
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << std::endl
      << "-------------------------------" << "<br>"
      << "-- Robust Resection of camera index: <" << viewIndex << "> image: "
      <<  view_I->s_Img_path <<"<br>"
      << "-- Threshold: " << errorMax << "<br>"
      << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "<br>"
      << "-- Nb points used for Resection: " << vec_featIdForResection.size() << "<br>"
      << "-- Nb points validated by robust estimation: " << vec_inliers.size() << "<br>"
      << "-- % points validated: "
      << vec_inliers.size()/static_cast<float>(vec_featIdForResection.size()) << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());
  }

  if (!bResection)
    return false;

  // D. Refine the pose of the found camera.
  // We use a local scene with only the 3D points and the nex camera.
  {
    // Decompose P matrix
    Mat3 K_, R_;
    Vec3 t_;
    KRt_From_P(P, &K_, &R_, &t_);

    // Create a SfM_DataScene with one camera and the 3D points
    SfM_Data tiny_scene;
    tiny_scene.views[view_I->id_view] = _sfm_data.GetViews().at(viewIndex);
    tiny_scene.poses[view_I->id_pose] = Pose3(R_, -R_.transpose() * t_);
    if (bKnownIntrinsic)
    {
      tiny_scene.intrinsics[view_I->id_intrinsic] = iterIntrinsic_I->second;
    }
    else
    {
      if (view_I->id_intrinsic == UndefinedIndexT)
      {
        // Update id_intrinsic to a valid value
        View * view_I = _sfm_data.GetViews().at(viewIndex).get();
        view_I->id_intrinsic = _sfm_data.GetIntrinsics().size();
      }
      // Create the new camera intrinsic group
      switch (_camType)
      {
        case PINHOLE_CAMERA:
          tiny_scene.intrinsics[view_I->id_intrinsic] =
            std::make_shared<Pinhole_Intrinsic>
            (view_I->ui_width, view_I->ui_height, K_(0,0), K_(0,2), K_(1,2));
        break;
        case PINHOLE_CAMERA_RADIAL1:
          tiny_scene.intrinsics[view_I->id_intrinsic] =
            std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (view_I->ui_width, view_I->ui_height, K_(0,0), K_(0,2), K_(1,2));
        break;
        case PINHOLE_CAMERA_RADIAL3:
          tiny_scene.intrinsics[view_I->id_intrinsic] =
            std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (view_I->ui_width, view_I->ui_height, K_(0,0), K_(0,2), K_(1,2));
        break;
        case PINHOLE_CAMERA_BROWN:
          tiny_scene.intrinsics[view_I->id_intrinsic] =
            std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (view_I->ui_width, view_I->ui_height, K_(0,0), K_(0,2), K_(1,2));
        break;
        default:
          std::cerr << "Try to create an unknown camera type." << std::endl;
          return false;
      }
    }
    // Insert the 3D points (landmarks) into the tiny_scene
    for (size_t i = 0; i < vec_inliers.size(); ++i)
    {
      const size_t idx = vec_inliers[i];
      Landmark landmark;
      landmark.X = pt3D.col(idx);
      const size_t feat_id = vec_featIdForResection[idx];
      Observation ob(_features_provider->feats_per_view[viewIndex][feat_id].coords().cast<double>(), feat_id);
      landmark.obs[view_I->id_view] = ob;
      tiny_scene.structure[i] = landmark;
    }
    Bundle_Adjustment_Ceres::BA_options options(true, false);
    options._linear_solver_type = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.Adjust(tiny_scene, true, true, false, false))
    {
      std::cerr << "Resection failed during the Bundle Adjustment." << std::endl;
      return false;
    }

    // E. Update the global scene with the new camera
    _sfm_data.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
    // - with new intrinsic data if it was unknown
    if (!bKnownIntrinsic)
    {
      _sfm_data.intrinsics[view_I->id_intrinsic] = tiny_scene.intrinsics[view_I->id_intrinsic];
    }
    _map_ACThreshold.insert(std::make_pair(viewIndex, errorMax));
  }

  // F. Update the observations into the global scene structure
  // F1. Add the new 2D observations to the reconstructed tracks
  // F2. Remove outliers from the tracks
  cpt = 0;
  std::vector<size_t>::iterator iterfeatId = vec_featIdForResection.begin();
  for (std::set<size_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
    iterTrackId != set_trackIdForResection.end(); ++iterTrackId, ++iterfeatId, ++cpt)
  {
    const IntrinsicBase * intrinsic = _sfm_data.GetIntrinsics().at(view_I->id_intrinsic).get();
    const Pose3 pose = _sfm_data.GetPoseOrDie(view_I);
    const Vec3 X = pt3D.col(cpt);
    const Vec2 x = _features_provider->feats_per_view[viewIndex][*iterfeatId].coords().cast<double>();
    const Vec2 residual = intrinsic->residual(pose, X, x);
    if (residual.norm() < errorMax &&
        pose.depth(X) > 0)
    {
      // Inlier, add the point to the reconstructed track
      _sfm_data.structure[*iterTrackId].obs[viewIndex] = Observation(x, *iterfeatId);
    }
    else {
      // Remove this observation from the scene tracking data
      _map_tracks[*iterTrackId].erase(viewIndex);
      // Remove the track itself, if there is only one view left.
      if (_map_tracks[*iterTrackId].size() < 2) {
          _map_tracks.erase(*iterTrackId);
      }
    }
  }

  // G. Triangulate new possible 2D tracks
  // We have tracks with only 2D observations but with 2 cameras reconstructed,
  // so we can triangulate these new points.
  {
    // For all reconstructed images look for common content in the tracks.
    const std::set<IndexT> valid_views = Get_Valid_Views(_sfm_data);
    for (std::set<IndexT>::const_iterator iterI = valid_views.begin();
      iterI != valid_views.end(); ++iterI)
    {
      const size_t indexI = *iterI;
      // Ignore the current view
      if (indexI == viewIndex) {  continue; }

      const size_t I = std::min(viewIndex, indexI);
      const size_t J = std::max(viewIndex, indexI);

      // Find matches between I and J
      // map_tracksCommon: All common tracks between I and J
      map_tracksCommon.clear();
      set_viewIndex.clear();
      set_viewIndex.insert(I); set_viewIndex.insert(J);
      TracksUtilsMap::GetTracksInImages(set_viewIndex, _map_tracks, map_tracksCommon);

      if (map_tracksCommon.empty()) { continue; } // no common tracks

      // All common tracks between I and J
      set_tracksIds.clear();
      TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

      // All tracks seen in I and J views but not already reconstructed in 3D.
      std::vector<IndexT> vec_tracksToAdd;
      std::set_difference(set_tracksIds.begin(), set_tracksIds.end(),
        reconstructed_trackId.begin(), reconstructed_trackId.end(),
        back_inserter(vec_tracksToAdd));

      // Do we have new tracks to add to the scene?
      if (vec_tracksToAdd.empty())
        continue;

      const View * view_1 = _sfm_data.GetViews().at(I).get();
      const View * view_2 = _sfm_data.GetViews().at(J).get();
      const IntrinsicBase * cam_1 = _sfm_data.GetIntrinsics().at(view_1->id_intrinsic).get();
      const IntrinsicBase * cam_2 = _sfm_data.GetIntrinsics().at(view_2->id_intrinsic).get();
      const Pose3 pose_1 = _sfm_data.GetPoseOrDie(view_1);
      const Pose3 pose_2 = _sfm_data.GetPoseOrDie(view_2);
      const Mat34 P1 = cam_1->get_projective_equivalent(pose_1);
      const Mat34 P2 = cam_2->get_projective_equivalent(pose_2);

      // All 2D/2D matches to triangulate
      IndMatches vec_index;
      TracksUtilsMap::TracksToIndexedMatches(map_tracksCommon, vec_tracksToAdd, &vec_index);
      assert(vec_index.size() == vec_tracksToAdd.size());

      // Triangulate 2D tracks
      IndexT new_track_count = 0;
      for (size_t i=0; i < vec_index.size(); ++i)
      {
        const size_t trackId = vec_tracksToAdd[i];

        // Get corresponding points and triangulate it
        const Vec2 x1 = _features_provider->feats_per_view[I][vec_index[i]._i].coords().cast<double>();
        const Vec2 x2 = _features_provider->feats_per_view[J][vec_index[i]._j].coords().cast<double>();

        Vec3 X_euclidean = Vec3::Zero();
        assert(reconstructed_trackId.count(trackId) == 0);
        const Vec2 x1_ud = cam_1->get_ud_pixel(x1);
        const Vec2 x2_ud = cam_2->get_ud_pixel(x2);
        TriangulateDLT(P1, x1_ud, P2, x2_ud, &X_euclidean);

        // Check triangulation results
        //  - Check angle (small angle leads imprecise triangulation)
        //  - Check positive depth
        //  - Check residual values
        const double angle = AngleBetweenRay(pose_1, cam_1, pose_2, cam_2, x1, x2);
        const Vec2 residual_1 = cam_1->residual(pose_1, X_euclidean, x1);
        const Vec2 residual_2 = cam_2->residual(pose_2, X_euclidean, x2);
        if ( angle > 2.0 &&
             pose_1.depth(X_euclidean) > 0 &&
             pose_2.depth(X_euclidean) > 0 &&
             residual_1.norm() < std::max(4.0, _map_ACThreshold[I]) &&
             residual_2.norm() < std::max(4.0, _map_ACThreshold[J]))
        {
          // Add a new track
          ++new_track_count;
          _sfm_data.structure[trackId].X = X_euclidean;
          _sfm_data.structure[trackId].obs[I] = Observation(x1, vec_index[i]._i);
          _sfm_data.structure[trackId].obs[J] = Observation(x2, vec_index[i]._j);
          reconstructed_trackId.insert(trackId);
        }
      }
      std::cout << "--Triangulated 3D points [" << I << "-" << J <<"]: "
        << "\t #Validated/#Possible: " << new_track_count << "/" << vec_index.size() << std::endl
        <<" #3DPoint for the entire scene: " << _sfm_data.GetLandmarks().size() << std::endl;
    }
  }
  return true;
}

/// Bundle adjustment to refine Structure; Motion and Intrinsics
void SequentialSfMReconstructionEngine::BundleAdjustment()
{
  Bundle_Adjustment_Ceres::BA_options options;
  if (_sfm_data.GetPoses().size() > 100)
  {
    options._preconditioner_type = ceres::JACOBI;
    options._linear_solver_type = ceres::SPARSE_SCHUR;
  }
  else
  {
    options._linear_solver_type = ceres::DENSE_SCHUR;
  }
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  bool bOk = bundle_adjustment_obj.Adjust(_sfm_data, true, true, !_bFixedIntrinsics);
}

/**
 * @brief Discard tracks with too large residual error
 *
 * Remove observation/tracks that have:
 *  - too large residual error
 *  - too small angular value
 *
 * @return True if more than 'count' outliers have been removed.
 */
size_t SequentialSfMReconstructionEngine::badTrackRejector(double dPrecision, size_t count)
{
  const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError( _sfm_data, dPrecision, 2);
  const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(_sfm_data, 2.0);

  return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
}

} // namespace sfm
} // namespace openMVG

