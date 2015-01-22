
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/SfM/SfMIncrementalEngine.hpp"
#include "software/SfM/SfMRobust.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
// Bundle Adjustment includes
#include "openMVG/bundle_adjustment/pinhole_brown_Rt_ceres_functor.hpp"
#include "openMVG/bundle_adjustment/problem_data_container.hpp"

//-- Feature tracking
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::tracks;

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"
#include "third_party/stlAddition/stlMap.hpp"
#include "openMVG/matching/indexed_sort.hpp"

#include <numeric>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <sstream>


namespace openMVG{

typedef SIOPointFeature FeatureT;
typedef std::vector<FeatureT> featsT;

IncrementalReconstructionEngine::IncrementalReconstructionEngine(const std::string & sImagePath,
  const std::string & sMatchesPath, const std::string & sOutDirectory, bool bHtmlReport)
  : ReconstructionEngine(sImagePath, sMatchesPath, sOutDirectory),
  _initialpair(std::pair<size_t,size_t>(0,0)),
  _bRefinePPandDisto(true),
  _bUseBundleAdjustment(true)
{
  _bHtmlReport = bHtmlReport;
  if (!stlplus::folder_exists(sOutDirectory)) {
    stlplus::folder_create(sOutDirectory);
  }
  if (_bHtmlReport)
  {
    _htmlDocStream = auto_ptr<htmlDocument::htmlDocumentStream>(
      new htmlDocument::htmlDocumentStream("openMVG A Contrario Incremental SFM report."));

    _htmlDocStream->pushInfo(
      htmlDocument::htmlMarkup("h1","openMVG A Contrario Incremental SFM report."));

    _htmlDocStream->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("Current directory: ") +
      sImagePath));
    _htmlDocStream->pushInfo("<hr>");
  }
}

IncrementalReconstructionEngine::~IncrementalReconstructionEngine()
{
  if (_bHtmlReport)
  {
    ofstream htmlFileStream( string(stlplus::folder_append_separator(_sOutDirectory) +
      "Reconstruction_Report.html").c_str());
    htmlFileStream << _htmlDocStream->getDoc();
  }
}

void pauseProcess()
{
  unsigned char i;
  std::cout << "\nPause : type key and press enter: ";
  cin>>i;
}

bool IncrementalReconstructionEngine::Process()
{
  //-------------------
  // Load data
  //-------------------

  if(!ReadInputData())
    return false;

  //-------------------
  //-- Incremental reconstruction
  //-------------------
  std::pair<size_t,size_t> initialPairIndex;
  if(! InitialPairChoice(initialPairIndex))
    return false;

  // Initial pair Essential Matrix and [R|t] estimation.
  if(! MakeInitialPair3D(initialPairIndex))
    return false;

  BundleAdjustment(); // Adjust 3D point and camera parameters.

  size_t round = 0;
  // Compute robust Resection of remaining image
  std::vector<size_t> vec_possible_resection_indexes;
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    bool bImageAdded = Resection(vec_possible_resection_indexes);
    if (bImageAdded)
    {
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << round << "_Resection";
      _reconstructorData.exportToPly( stlplus::create_filespec(_sOutDirectory, os.str(), ".ply"));
    }
    ++round;
    if (bImageAdded && _bUseBundleAdjustment)
    {
      // Perform BA until all point are under the given precision
      do
      {
        BundleAdjustment();
      }
      while (badTrackRejector(4.0) != 0);
    }
  }

  //-- Reconstruction done.
  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << _reconstructorData.map_Camera.size()
    << " from " << _vec_camImageNames.size() << " input images.\n"
    << "-- #Tracks, #3D points: " << _reconstructorData.map_3DPoints.size() << "\n"
    << "-------------------------------" << "\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  std::cout << "\nHistogram of residuals:" << h.ToString() << std::endl;

  if (_bHtmlReport)
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion process finished.";
    _htmlDocStream->pushInfo("<hr>");
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Structure from Motion (statistics):<br>"
      << "-- #Camera calibrated: " << _reconstructorData.map_Camera.size()
      << " from " <<_vec_camImageNames.size() << " input images.<br>"
      << "-- #Tracks, #3D points: " << _reconstructorData.map_3DPoints.size() << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());

    _htmlDocStream->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

    std::vector<double> xBin = h.GetXbinsValue();
    std::pair< std::pair<double,double>, std::pair<double,double> > range;
    range = autoJSXGraphViewport<double>(xBin, h.GetHist());

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

bool IncrementalReconstructionEngine::ReadInputData()
{
  if (!stlplus::is_folder(_sImagePath) ||
    !stlplus::is_folder(_sMatchesPath) ||
    !stlplus::is_folder(_sOutDirectory))
  {
    std::cerr << std::endl
      << "One of the required directory is not a valid directory" << std::endl;
    return false;
  }

  const std::string sListsFile = stlplus::create_filespec(_sMatchesPath,"lists","txt");
  const std::string sComputedMatchesFile_F = stlplus::create_filespec(_sMatchesPath,"matches.f","txt");
  if (!stlplus::is_file(sListsFile)||
    !stlplus::is_file(sComputedMatchesFile_F) )
  {
    std::cerr << std::endl
      << "One of the input required file is not a present (lists.txt, matches.f.txt)" << std::endl;
    return false;
  }

  // a. Read images names
  {
    if (!openMVG::SfMIO::loadImageList( _vec_camImageNames,
                                        _vec_intrinsicGroups,
                                        sListsFile) )
    {
      std::cerr << "\nEmpty image list." << std::endl;
      return false;
    }
    else
    {
      // Find to which intrinsic groups each image belong
      for (std::vector<openMVG::SfMIO::CameraInfo>::const_iterator iter = _vec_camImageNames.begin();
        iter != _vec_camImageNames.end(); ++iter)
      {
        const openMVG::SfMIO::CameraInfo & camInfo = *iter;
        // Find the index of the camera
        size_t idx = std::distance((std::vector<openMVG::SfMIO::CameraInfo>::const_iterator)_vec_camImageNames.begin(), iter);
        _set_remainingImageId.insert(idx);
        _map_IntrinsicIdPerImageId[idx] = camInfo.m_intrinsicId;
      }
    }
  }

  // b. Read matches (Fundamental)
  if (!matching::PairedIndMatchImport(sComputedMatchesFile_F, _map_Matches_F)) {
    std::cerr<< "Unable to read the Fundamental matrix matches" << std::endl;
    return false;
  }

  // c. Compute tracks from matches
  TracksBuilder tracksBuilder;

  {
    std::cout << std::endl << "Track building" << std::endl;
    tracksBuilder.Build(_map_Matches_F);
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
      TracksUtilsMap::ImageIdInTracks(_map_tracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(_map_tracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)  {
        osTrack << "\t" << iter->first << "\t" << iter->second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Read features:
  for (size_t i = 0; i < _vec_camImageNames.size(); ++i)  {
    const size_t camIndex = i;
    if (!loadFeatsFromFile(
      stlplus::create_filespec(_sMatchesPath, stlplus::basename_part(_vec_camImageNames[camIndex].m_sImageName), ".feat"),
      _map_feats[camIndex])) {
      std::cerr << "Bad reading of feature files" << std::endl;
      return false;
    }
  }

  if (_bHtmlReport)
  {
    _htmlDocStream->pushInfo( "Dataset info:");
    _htmlDocStream->pushInfo( "Number of images: " +
      htmlDocument::toString(_vec_camImageNames.size()) + "<br>");

    std::ostringstream osTrack;
    {
      //-- Display stats :
      //    - number of images
      //    - number of tracks
      std::set<size_t> set_imagesId;
      TracksUtilsMap::ImageIdInTracks(_map_tracks, set_imagesId);
      osTrack << "------------------" << "<br>"
        << "-- Tracks Stats --" << "<br>"
        << " Tracks number: " << tracksBuilder.NbTracks() << "<br>"
        << " Images Id: " << "<br>";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "<br>------------------" << "<br>";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(_map_tracks, map_Occurence_TrackLength);
      osTrack << "<table> <tr> <td>TrackLength, </td> <td>Occurrence<td> </tr>";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)
      {
        osTrack << "<tr><td>" << iter->first << "</td><td>"
          << iter->second << "</td></tr>";
      }
      osTrack << "</table>";
    }

    _htmlDocStream->pushInfo(osTrack.str());
    _htmlDocStream->pushInfo("<br>");

    _htmlDocStream->pushInfo("<hr>");
  }
  return true;
}

/// Find the best initial pair
bool IncrementalReconstructionEngine::InitialPairChoice( std::pair<size_t, size_t> & initialPairIndex)
{
  if (_initialpair != std::pair<size_t,size_t>(0,0))
  {
    initialPairIndex = _initialpair;
  }
  else
  {
    std::cout << std::endl
      << "---------------------------------------------------\n"
      << "IncrementalReconstructionEngine::InitialPairChoice\n"
      << "---------------------------------------------------\n"
      << " The best F matrix validated pairs are displayed\n"
      << " Choose one pair manually by typing the two integer indexes\n"
      << "---------------------------------------------------\n"
      << std::endl;

    // Display to the user the 10 top Fundamental matches pair
    std::vector< size_t > vec_NbMatchesPerPair;
    for (openMVG::matching::PairWiseMatches::const_iterator
      iter = _map_Matches_F.begin();
      iter != _map_Matches_F.end(); ++iter)
    {
      vec_NbMatchesPerPair.push_back(iter->second.size());
    }
    // sort in descending order
    using namespace indexed_sort;
    std::vector< sort_index_packet_descend< size_t, size_t> > packet_vec(vec_NbMatchesPerPair.size());
    sort_index_helper(packet_vec, &vec_NbMatchesPerPair[0], std::min((size_t)10, _map_Matches_F.size()));

    for (size_t i = 0; i < std::min((size_t)10, _map_Matches_F.size()); ++i) {
      size_t index = packet_vec[i].index;
      openMVG::matching::PairWiseMatches::const_iterator iter = _map_Matches_F.begin();
      std::advance(iter, index);
      std::cout << "(" << iter->first.first << "," << iter->first.second <<")\t\t"
        << iter->second.size() << " matches" << std::endl;
    }

    // Manual choice of the initial pair
    std::cout << std::endl << " type INITIAL pair Indexes: X enter Y enter\n";
    int val, val2;
    if ( std::cin >> val && std::cin >> val2) {
      initialPairIndex.first = val;
      initialPairIndex.second = val2;
    }
  }

  std::cout << "\nPutative starting pair is: (" << initialPairIndex.first
      << "," << initialPairIndex.second << ")" << std::endl;

  // Check validity of the initial pair indices:
  if (_map_feats.find(initialPairIndex.first) == _map_feats.end() ||
      _map_feats.find(initialPairIndex.second) == _map_feats.end())
  {
    std::cerr << "At least one of the initial pair indices is invalid."
      << std::endl;
    return false;
  }
  return true;
}

bool IncrementalReconstructionEngine::MakeInitialPair3D(const std::pair<size_t,size_t> & initialPair)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  const size_t I = min(initialPair.first,initialPair.second);
  const size_t J = max(initialPair.first,initialPair.second);

  // Check validity of the initial pair indices:
  if (_map_feats.find(I) == _map_feats.end() ||
       _map_feats.find(J) == _map_feats.end())
  {
    std::cerr << "At least one of the initial pair indices is invalid."
      << std::endl;
    return false;
  }

  // a.coords Get common tracks between the two images
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  std::set<size_t> set_imageIndex;
  set_imageIndex.insert(I);
  set_imageIndex.insert(J);
  TracksUtilsMap::GetTracksInImages(set_imageIndex, _map_tracks, map_tracksCommon);

  // b. Get corresponding features
  std::vector<SIOPointFeature> & vec_featI = _map_feats[I];
  std::vector<SIOPointFeature> & vec_featJ = _map_feats[J];

  //-- Copy point to array in order to estimate essential matrix:
  const size_t n = map_tracksCommon.size();
  Mat x1(2,n), x2(2,n);
  size_t cptIndex = 0;
  for (openMVG::tracks::STLMAPTracks::const_iterator
    iterT = map_tracksCommon.begin();
    iterT != map_tracksCommon.end();
    ++iterT)
  {
    //Get corresponding point
    tracks::submapTrack::const_iterator iter = iterT->second.begin();
    const size_t i = iter->second;
    const size_t j = (++iter)->second;
    const SIOPointFeature & imaL = vec_featI[i];
    const SIOPointFeature & imaR = vec_featJ[j];

    x1.col(cptIndex) = imaL.coords().cast<double>();
    x2.col(cptIndex) = imaR.coords().cast<double>();
    ++cptIndex;
  }

  Mat3 E = Mat3::Identity();
  std::vector<size_t> vec_inliers;
  double errorMax = std::numeric_limits<double>::max();

  const openMVG::SfMIO::IntrinsicCameraInfo & intrinsicCamI = _vec_intrinsicGroups[_vec_camImageNames[I].m_intrinsicId];
  const openMVG::SfMIO::IntrinsicCameraInfo & intrinsicCamJ = _vec_intrinsicGroups[_vec_camImageNames[J].m_intrinsicId];

  if (!intrinsicCamI.m_bKnownIntrinsic ||
      !intrinsicCamJ.m_bKnownIntrinsic)
  {
    std::cerr << "The selected image pair doesn't have intrinsic parameters, we cannot estimate an initial reconstruction from the two provided index: [I,J]=[" << I <<"," << J << "]" << std::endl;
    return false;
  }

  if (!SfMRobust::robustEssential(intrinsicCamI.m_K, intrinsicCamJ.m_K,
    x1, x2,
    &E, &vec_inliers,
    std::make_pair(intrinsicCamI.m_w, intrinsicCamI.m_h),
    std::make_pair(intrinsicCamJ.m_w, intrinsicCamJ.m_h),
    &errorMax))
  {
    std::cerr << " /!\\ Robust estimation failed to compute E for the initial pair" << std::endl;
    return false;
  }

  std::cout << std::endl
    << "-- Robust Essential Matrix estimation " << std::endl
    << "--  #tentative/#inliers: " << x1.cols() << " / " << vec_inliers.size() << std::endl
    << "--  Threshold: " << errorMax << std::endl;

  //--> Estimate the best possible Rotation/Translation from E
  Mat3 RJ;
  Vec3 tJ;
  if (!SfMRobust::estimate_Rt_fromE(intrinsicCamI.m_K, intrinsicCamJ.m_K,
        x1, x2, E, vec_inliers,
        &RJ, &tJ))
  {
    std::cout << " /!\\ Failed to compute initial R|t for the initial pair" << std::endl;
    return false;
  }
  std::cout << std::endl
    << "-- Rotation|Translation matrices: --" << std::endl
    << RJ << std::endl << std::endl << tJ << std::endl;

  //-> Triangulate the common tracks
  //--> Triangulate the point

  BrownPinholeCamera camI(intrinsicCamI.m_focal, intrinsicCamI.m_K(0,2), intrinsicCamI.m_K(1,2), Mat3::Identity(), Vec3::Zero());
  BrownPinholeCamera camJ(intrinsicCamJ.m_focal, intrinsicCamJ.m_K(0,2), intrinsicCamJ.m_K(1,2), RJ, tJ);

  std::vector<IndMatch> vec_index;
  for (openMVG::tracks::STLMAPTracks::const_iterator
    iterT = map_tracksCommon.begin();
    iterT != map_tracksCommon.end();
    ++iterT)
  {
    tracks::submapTrack::const_iterator iter = iterT->second.begin();
    tracks::submapTrack::const_iterator iter2 = iterT->second.begin();
    std::advance(iter2,1);
    vec_index.push_back(IndMatch(iter->second, iter2->second));
  }
  std::vector<Vec3> vec_3dPoint;
  std::vector<double> vec_triangulationResidual;

  SfMRobust::triangulate2View_Vector(camI._P, camJ._P,
    vec_featI, vec_featJ,
    vec_index, &vec_3dPoint, &vec_triangulationResidual);

  //- Add reconstructed point to the reconstruction data
  //- Write corresponding that the track have a corresponding 3D point
  cptIndex = 0;
  for (openMVG::tracks::STLMAPTracks::const_iterator
    iterT = map_tracksCommon.begin();
    iterT != map_tracksCommon.end();
    ++iterT, cptIndex++)
  {
    size_t trackId = iterT->first;
    const Vec2 xI = vec_featI[vec_index[cptIndex]._i].coords().cast<double>();
    const Vec2 xJ = vec_featJ[vec_index[cptIndex]._j].coords().cast<double>();

    if (find(vec_inliers.begin(), vec_inliers.end(), cptIndex) != vec_inliers.end())
    {
      //-- Depth
      //-- Residuals
      //-- Angle between rays ?
      if (camI.Depth(vec_3dPoint[cptIndex]) > 0
        && camJ.Depth(vec_3dPoint[cptIndex]) > 0)  {
        double angle = BrownPinholeCamera::AngleBetweenRay(camI, camJ, xI, xJ);
        if (angle > 2.)  {
          _reconstructorData.map_3DPoints[trackId] = vec_3dPoint[cptIndex];
          _reconstructorData.set_trackId.insert(trackId);
          _map_reconstructed[trackId].insert(make_pair(I,vec_index[cptIndex]._i));
          _map_reconstructed[trackId].insert(make_pair(J,vec_index[cptIndex]._j));
        }
      }
    }
    else  {
      //Remove this track entry with ImageIndexes
      _map_tracks[trackId].erase(I);
      _map_tracks[trackId].erase(J);
      if (_map_tracks[trackId].size() < 2)  {
        _map_tracks[trackId].clear();
        _map_tracks.erase(trackId);
      }
    }
  }

  std::cout << "--#Triangulated 3D points: " << vec_inliers.size() << "\n";
  std::cout << "--#Triangulated 3D points under threshold: " << _reconstructorData.map_3DPoints.size()  << "\n";
  std::cout << "--#Putative correspondences: " << x1.cols()  << "\n";

  _set_remainingImageId.erase(I);
  _set_remainingImageId.erase(J);

  if (!_reconstructorData.map_3DPoints.empty())
  {
    // Add information related to the View (I,J) to the reconstruction data
    _reconstructorData.set_imagedId.insert(I);
    _reconstructorData.set_imagedId.insert(J);

    _reconstructorData.map_Camera.insert(std::make_pair(I, camI));
    _reconstructorData.map_Camera.insert(std::make_pair(J, camJ));

    // Add the images to the reconstructed intrinsic group
    _map_ImagesIdPerIntrinsicGroup[_map_IntrinsicIdPerImageId[I]].push_back(I);
    _map_ImagesIdPerIntrinsicGroup[_map_IntrinsicIdPerImageId[J]].push_back(J);
    // Initialize the first intrinsics groups:
    Vec6 & intrinsicI = _map_IntrinsicsPerGroup[_map_IntrinsicIdPerImageId[I]];
    intrinsicI<< camI._f, camI._ppx, camI._ppy, camI._k1, camI._k2, camI._k3;
    Vec6 & intrinsicJ = _map_IntrinsicsPerGroup[_map_IntrinsicIdPerImageId[J]];
    intrinsicJ<< camJ._f, camJ._ppx, camJ._ppy, camJ._k1, camJ._k2, camJ._k3;

    _map_ACThreshold.insert(std::make_pair(I, errorMax));
    _map_ACThreshold.insert(std::make_pair(J, errorMax));

    _vec_added_order.push_back(I);
    _vec_added_order.push_back(J);
  }


  _reconstructorData.exportToPly(stlplus::create_filespec(_sOutDirectory,"sceneStart","ply"));

  Histogram<double> histoResiduals;
  std::cout << std::endl
    << "=========================\n"
    << " MSE Residual InitialPair Inlier: " << ComputeResidualsHistogram(&histoResiduals) << "\n"
    << "=========================" << std::endl;

  if (_bHtmlReport)
  {
    using namespace htmlDocument;
    _htmlDocStream->pushInfo(htmlMarkup("h1","Essential Matrix."));
    ostringstream os;
    os << std::endl
      << "-------------------------------" << "<br>"
      << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
      << stlplus::basename_part(_vec_camImageNames[I].m_sImageName) << ","
      << stlplus::basename_part(_vec_camImageNames[J].m_sImageName) << "<br>"
      << "-- Threshold: " << errorMax << "<br>"
      << "-- Resection status: " << "OK" << "<br>"
      << "-- Nb points used for robust Essential matrix estimation: "
        << map_tracksCommon.size() << "<br>"
      << "-- Nb points validated by robust estimation: "
        << vec_inliers.size() << "<br>"
      << "-- % points validated: "
      << vec_inliers.size()/static_cast<float>(map_tracksCommon.size())
        << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());

    _htmlDocStream->pushInfo(htmlMarkup("h2","Residual of the robust estimation (Initial triangulation). Thresholded at: " + toString(errorMax)));

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("InitialPairTriangulationResiduals",600,300);
    jsxGraph.addYChart(vec_triangulationResidual, "point");
    jsxGraph.addLine(0,errorMax, vec_triangulationResidual.size(), errorMax);
    jsxGraph.UnsuspendUpdate();
    std::pair< std::pair<double,double>, std::pair<double,double> > range = autoJSXGraphViewport<double>(vec_triangulationResidual);
    jsxGraph.setViewport(range);
    jsxGraph.close();
    _htmlDocStream->pushInfo(jsxGraph.toStr());


    _htmlDocStream->pushInfo(htmlMarkup("h2","Histogram of residuals"));

    std::vector<double> xBin = histoResiduals.GetXbinsValue();
    range = autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

    jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
    jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
    jsxGraph.addLine(errorMax,0, errorMax, histoResiduals.GetHist().front());
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    _htmlDocStream->pushInfo(jsxGraph.toStr());

    _htmlDocStream->pushInfo("<hr>");

    ofstream htmlFileStream( string(stlplus::folder_append_separator(_sOutDirectory) +
      "Reconstruction_Report.html").c_str());
    htmlFileStream << _htmlDocStream->getDoc();
  }
  return !_reconstructorData.map_3DPoints.empty();
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

/// List the images that the greatest number of matches to the current 3D reconstruction.
bool IncrementalReconstructionEngine::FindImagesWithPossibleResection(std::vector<size_t> & vec_possible_indexes)
{
  vec_possible_indexes.clear();

  if (_set_remainingImageId.empty() || _map_reconstructed.empty())  {
    return false;
  }

  // Estimate the image on which we could compute a resection safely
  // -> We first find the camera with the greatest number of matches
  //     with the current 3D existing 3D point => M
  // -> Then add any camera with at least 0.75M matches.
  // Keep only the best one.

  const double dPourcent = 0.75;

  std::vector< std::pair<size_t, size_t> > vec_putative; // ImageId, NbPutativeCommonPoint
  for (std::set<size_t>::const_iterator iter = _set_remainingImageId.begin();
        iter != _set_remainingImageId.end(); ++iter)
  {
    const size_t imageIndex = *iter;

    // Compute 2D - 3D possible content
    openMVG::tracks::STLMAPTracks map_tracksCommon;
    std::set<size_t> set_imageIndex;
    set_imageIndex.insert(imageIndex);
    TracksUtilsMap::GetTracksInImages(set_imageIndex, _map_tracks, map_tracksCommon);

    if (! map_tracksCommon.empty())
    {
      std::set<size_t> set_tracksIds;
      TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

      // Count the common possible putative point
      //  with the already 3D reconstructed trackId
      std::vector<size_t> vec_trackIdForResection;
      set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
        _reconstructorData.set_trackId.begin(),
        _reconstructorData.set_trackId.end(),
        std::back_inserter(vec_trackIdForResection));

      vec_putative.push_back( make_pair(imageIndex, vec_trackIdForResection.size()));
    }
  }

  if (vec_putative.empty()) {
    return false;
  }
  else  {
    sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<size_t, size_t, std::greater<size_t> >());

    size_t M = vec_putative[0].second;
    vec_possible_indexes.push_back(vec_putative[0].first);

    // TEMPORARY
    // std::cout << std::endl << std::endl << "TEMPORARY return only the best image" << std::endl;
    // return true;
    // END TEMPORARY

    const size_t threshold = static_cast<size_t>(dPourcent * M);
    for (size_t i = 1; i < vec_putative.size() &&
      vec_putative[i].second > threshold; ++i)
    {
      vec_possible_indexes.push_back(vec_putative[i].first);
    }
    return true;
  }
}

/// Add to the current scene the desired image indexes.
bool IncrementalReconstructionEngine::Resection(std::vector<size_t> & vec_possible_indexes)
{
  bool bOk = false;
  for (std::vector<size_t>::const_iterator iter = vec_possible_indexes.begin();
    iter != vec_possible_indexes.end();
    ++iter)
  {
    _vec_added_order.push_back(*iter);
    bool bResect = Resection(*iter);
    bOk |= bResect;
    if (!bResect) {
    // Resection was not possible (we remove the image from the remaining list)
      std::cerr << std::endl
        << "Resection of image: " << *iter << " was not possible" << std::endl;
    }
    _set_remainingImageId.erase(*iter);
  }
  return bOk;
}

/// Add a single Image to the scene and triangulate new possible tracks
bool IncrementalReconstructionEngine::Resection(size_t imageIndex)
{
  std::cout << std::endl
    << "-------------------------------" << std::endl
    << "-- Resection of camera index: " << imageIndex << std::endl
    << "-------------------------------" << std::endl;

  // Compute 2D - 3D possible content
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  std::set<size_t> set_imageIndex;
  set_imageIndex.insert(imageIndex);
  TracksUtilsMap::GetTracksInImages(set_imageIndex, _map_tracks, map_tracksCommon);

  std::set<size_t> set_tracksIds;
  TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

  // Intersect 3D reconstructed trackId with the one that contain the Image Id of interest
  std::set<size_t> set_trackIdForResection;
  set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
    _reconstructorData.set_trackId.begin(),
    _reconstructorData.set_trackId.end(),
    std::inserter(set_trackIdForResection, set_trackIdForResection.begin()));

  // Load feature corresponding to imageIndex
  const std::vector<SIOPointFeature> & vec_featsImageIndex = _map_feats[imageIndex];

  // Get back featId and tracksID that will be used for the resection
  std::vector<size_t> vec_featIdForResection;
  TracksUtilsMap::GetFeatIndexPerViewAndTrackId(map_tracksCommon,
    set_trackIdForResection,
    imageIndex,
    &vec_featIdForResection);

  std::cout << std::endl << std::endl
    << " Tracks in: " << imageIndex << std::endl
    << " \t" << map_tracksCommon.size() << std::endl
    << " #Reconstructed tracks:" << std::endl
    << " \t" << _reconstructorData.set_trackId.size() << std::endl
    << " #Tracks Valid for resection:" << std::endl
    << " \t" << set_trackIdForResection.size() << std::endl;

  // Normally it must not crash even if it have 0 matches
  if(set_trackIdForResection.empty())
  {
    // Too few matches (even 0 before....) images with empty connection
    _set_remainingImageId.erase(imageIndex);
    return false;
  }

  // Create pt2D, and pt3D array
  Mat pt2D( 2, set_trackIdForResection.size());
  Mat pt3D( 3, set_trackIdForResection.size());

  size_t cpt = 0;
  std::set<size_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
  for (std::vector<size_t>::const_iterator iterfeatId = vec_featIdForResection.begin();
    iterfeatId != vec_featIdForResection.end();
    ++iterfeatId, ++iterTrackId, ++cpt)
  {
    pt3D.col(cpt) = _reconstructorData.map_3DPoints[*iterTrackId];
    pt2D.col(cpt) = vec_featsImageIndex[*iterfeatId].coords().cast<double>();
  }

  //-------------
  std::vector<size_t> vec_inliers;
  Mat34 P;
  double errorMax = std::numeric_limits<double>::max();

  const openMVG::SfMIO::IntrinsicCameraInfo & intrinsicCam = _vec_intrinsicGroups[_vec_camImageNames[imageIndex].m_intrinsicId];

  bool bResection = SfMRobust::robustResection(
    std::make_pair( intrinsicCam.m_w, intrinsicCam.m_h ),
    pt2D, pt3D,
    &vec_inliers,
    // If intrinsics guess exist use it, else use a standard 6 points pose resection
    (intrinsicCam.m_bKnownIntrinsic == true) ? & intrinsicCam.m_K : NULL,
    &P, &errorMax, _bRefineFocal );

  std::cout << std::endl
    << "-------------------------------" << std::endl
    << "-- Robust Resection of camera index: " << imageIndex << std::endl
    << "-- Resection status: " << bResection << std::endl
    << "-- #Points used for Resection: " << vec_featIdForResection.size() << std::endl
    << "-- #Points validated by robust Resection: " << vec_inliers.size() << std::endl
    << "-- Threshold: " << errorMax << std::endl
    << "-------------------------------" << std::endl;

  if (_bHtmlReport)
  {
    using namespace htmlDocument;
    ostringstream os;
    os << "Resection of Image index: <" << imageIndex << "> image: "
      << stlplus::basename_part(_vec_camImageNames[imageIndex].m_sImageName) <<"<br> \n";
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << std::endl
      << "-------------------------------" << "<br>"
      << "-- Robust Resection of camera index: <" << imageIndex << "> image: "
      << stlplus::basename_part(_vec_camImageNames[imageIndex].m_sImageName) <<"<br>"
      << "-- Threshold: " << errorMax << "<br>"
      << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "<br>"
      << "-- Nb points used for Resection: " << vec_featIdForResection.size() << "<br>"
      << "-- Nb points validated by robust estimation: " << vec_inliers.size() << "<br>"
      << "-- % points validated: "
      << vec_inliers.size()/static_cast<float>(vec_featIdForResection.size()) << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());
  }

  if (!bResection) {
    return false;
  }

  //-- Add the camera to the _reconstruction data.
  _reconstructorData.set_imagedId.insert(imageIndex);
  {
    Mat3 K, R;
    Vec3 t;
    KRt_From_P(P, &K, &R, &t);
    std::cout << "\n Resection Calibration Matrix \n" << K << std::endl << std::endl;

    BrownPinholeCamera cam(K(0,0), K(0,2), K(1,2), R, t);
    _reconstructorData.map_Camera.insert(std::make_pair(imageIndex, cam));

    // Add the image to the reconstructed intrinsic group
    _map_ImagesIdPerIntrinsicGroup[_map_IntrinsicIdPerImageId[imageIndex]].push_back(imageIndex);

    // if intrinsic groups is empty fill a new:
    if (_map_IntrinsicsPerGroup.find(_map_IntrinsicIdPerImageId[imageIndex])  == _map_IntrinsicsPerGroup.end())
    {
      Vec6 & intrinsic = _map_IntrinsicsPerGroup[_map_IntrinsicIdPerImageId[imageIndex]];
      intrinsic << cam._f, cam._ppx, cam._ppy, cam._k1, cam._k2, cam._k3;
    }
  }
  _map_ACThreshold.insert(std::make_pair(imageIndex, errorMax));
  _set_remainingImageId.erase(imageIndex);

  // Evaluate residuals:
  std::vector<double> vec_ResectionResidual;
  for (size_t i = 0; i < pt3D.cols(); ++i)
  {
    double dResidual = PinholeCamera::Residual(P, pt3D.col(i), pt2D.col(i));
    vec_ResectionResidual.push_back(dResidual);
  }

  if (_bHtmlReport)
  {
    using namespace htmlDocument;
    // Export graphical residual statistics
    _htmlDocStream->pushInfo(htmlMarkup("h2","Residual of the robust estimation (Resection). Thresholded at: " + toString(errorMax)));
    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init(std::string("ResectionResidual_residualPlot_" + toString(imageIndex)),600,300);
    jsxGraph.addYChart(vec_ResectionResidual, "point");
    jsxGraph.addLine(0,errorMax, vec_ResectionResidual.size(), errorMax);
    jsxGraph.UnsuspendUpdate();
    std::pair< std::pair<double,double>, std::pair<double,double> > range = autoJSXGraphViewport<double>(vec_ResectionResidual);
    jsxGraph.setViewport(range);
    jsxGraph.close();
    _htmlDocStream->pushInfo(jsxGraph.toStr());

    if(bResection)  {

      Histogram<double> histo(0, 2*errorMax, 10);
      histo.Add(vec_ResectionResidual.begin(), vec_ResectionResidual.end());
      std::vector<double> xBin = histo.GetXbinsValue();
      std::pair< std::pair<double,double>, std::pair<double,double> > range = autoJSXGraphViewport<double>(xBin, histo.GetHist());

      jsxGraph.init(std::string("ResectionResidual_residualHisto_" + toString(imageIndex)),600,300);
      jsxGraph.addXYChart(xBin, histo.GetHist(), "line,point");
      jsxGraph.addLine(errorMax,0, errorMax, histo.GetHist().front());
      jsxGraph.UnsuspendUpdate();
      jsxGraph.setViewport(range);
      jsxGraph.close();
      _htmlDocStream->pushInfo(jsxGraph.toStr());
    }
    _htmlDocStream->pushInfo("<hr>");
  }

  // Add new entry to reconstructed track and
  //  remove outlier from the tracks
  cpt = 0;
  std::vector<size_t>::iterator iterfeatId = vec_featIdForResection.begin();
  for (std::set<size_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
    iterTrackId != set_trackIdForResection.end(); ++iterTrackId, ++cpt, ++iterfeatId)
  {
    if (vec_ResectionResidual[cpt] < errorMax) {
      // Inlier, add the point to the reconstructed track
      _map_reconstructed[*iterTrackId].insert(make_pair(imageIndex,*iterfeatId));
    }
    else {
      // Outlier remove this entry from the tracks
      _map_tracks[*iterTrackId].erase(imageIndex);
    }
  }

  // Add new possible tracks (triangulation)
  // Triangulate new possible tracks:
  // For all Union [ CurrentId, [PreviousReconstructedIds] ]
  //   -- If trackId not yet registered:
  //      -- Triangulate and add tracks id.

  {
    // For all reconstructed image look if common content in the track
    for (std::set<size_t>::const_iterator iterI = _reconstructorData.set_imagedId.begin();
      iterI != _reconstructorData.set_imagedId.end(); ++iterI)
    {
      const size_t & indexI = *iterI;
      if (indexI == imageIndex) {  continue; }
      size_t I = std::min(imageIndex, indexI);
      size_t J = std::max(imageIndex, indexI);

      // Compute possible content (match between indexI, indexJ)
      map_tracksCommon.clear(); set_imageIndex.clear();
      set_imageIndex.insert(I); set_imageIndex.insert(J);
      TracksUtilsMap::GetTracksInImages(set_imageIndex, _map_tracks, map_tracksCommon);

      if (map_tracksCommon.empty()) { continue; } // no common content

      set_tracksIds.clear();
      TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

      //-- Compute if we have something to add to the scene ?
      std::vector<size_t> vec_tracksToAdd;
      //-- Do we have new Track to add ?
      set_difference(set_tracksIds.begin(), set_tracksIds.end(),
        _reconstructorData.set_trackId.begin(), _reconstructorData.set_trackId.end(),
        back_inserter(vec_tracksToAdd));

      if (!vec_tracksToAdd.empty())
      {
        const Mat34 & P1 = _reconstructorData.map_Camera.find(I)->second._P;
        const Mat34 & P2 = _reconstructorData.map_Camera.find(J)->second._P;

        const std::vector<SIOPointFeature> & vec_featI = _map_feats[I];
        const std::vector<SIOPointFeature> & vec_featJ = _map_feats[J];

        std::vector<IndMatch> vec_index;
        TracksUtilsMap::TracksToIndexedMatches(map_tracksCommon, vec_tracksToAdd, &vec_index);

        std::vector<Vec3> vec_3dPoint;
        std::vector<double> vec_triangulationResidual;
        vec_triangulationResidual.reserve(vec_index.size());
        vec_3dPoint.reserve(vec_index.size());
        SfMRobust::triangulate2View_Vector(P1, P2,
          vec_featI, vec_featJ,
          vec_index, &vec_3dPoint, &vec_triangulationResidual);

        bool bVisual = false; // Turn to true to see possible new triangulated points
        if (bVisual)
        {
          ostringstream os;
          os << "scene_" << I << "-" << J;
          plyHelper::exportToPly(vec_3dPoint,
            stlplus::create_filespec(_sOutDirectory, os.str(), "ply"));
        }

        // Analyze 3D reconstructed point
        //  - Check positive depth
        //  - Check angle (small angle leads imprecise triangulation)
        const BrownPinholeCamera & cam1 = _reconstructorData.map_Camera.find(I)->second;
        const BrownPinholeCamera & cam2 = _reconstructorData.map_Camera.find(J)->second;

        // Threshold associated to each camera (min value of 4 pixels, to let BA refine if necessary)
        const double maxThI = std::max(4.0, _map_ACThreshold[I]);
        const double maxThJ = std::max(4.0, _map_ACThreshold[J]);

        //- Add reconstructed point to the reconstruction data
        size_t cardPointsBefore = _reconstructorData.map_3DPoints.size();
        for (size_t i = 0; i < vec_tracksToAdd.size(); ++i)
        {
          const size_t trackId = vec_tracksToAdd[i];
          const Vec3 & cur3DPt = vec_3dPoint[i];

          if ( _reconstructorData.set_trackId.find(trackId) == _reconstructorData.set_trackId.end())
          {
            const Vec2 x1 = vec_featI[vec_index[i]._i].coords().cast<double>();
            const Vec2 x2 = vec_featJ[vec_index[i]._j].coords().cast<double>();

            bool bReproj =
              cam1.Residual(cur3DPt, x1) < maxThI &&
              cam2.Residual(cur3DPt, x2) < maxThJ;

            if ( bReproj
                && cam1.Depth(cur3DPt) > 0
                && cam2.Depth(cur3DPt) > 0)
            {
              double angle = BrownPinholeCamera::AngleBetweenRay(cam1, cam2, x1, x2);
              if(angle>2) {
                _reconstructorData.map_3DPoints[trackId] = vec_3dPoint[i];
                _reconstructorData.set_trackId.insert(trackId);
                _map_reconstructed[trackId].insert(make_pair(I, vec_index[i]._i));
                _map_reconstructed[trackId].insert(make_pair(J, vec_index[i]._j));
              }
            }
          }
        }

        std::cout << "--Triangulated 3D points [" << I << "-" << J <<"]: "
          << "\t #Validated/#Possible: " << _reconstructorData.map_3DPoints.size() - cardPointsBefore
          << "/" << vec_3dPoint.size() << std::endl
          <<" #3DPoint for the entire scene: " << _reconstructorData.set_trackId.size() << std::endl;

        if(bVisual) {
          ostringstream sFileName;
          sFileName << "incremental_" << indexI << "-" << imageIndex;
          _reconstructorData.exportToPly(stlplus::create_filespec(_sOutDirectory, sFileName.str(), "ply"));
        }
      }
    }
  }
  return true;
}

size_t IncrementalReconstructionEngine::badTrackRejector(double dPrecision)
{
  // Go through the track and look for too large residual

  std::set<size_t> set_camIndex;
  std::transform(_reconstructorData.map_Camera.begin(),
    _reconstructorData.map_Camera.end(),
    std::inserter(set_camIndex,set_camIndex.begin()),
    RetrieveKey());

  std::map<size_t, std::set<size_t> > map_trackToErase; // trackid, imageIndexes
  std::set<size_t> set_trackToErase;

  for (std::map<size_t,Vec3>::const_iterator iter = _reconstructorData.map_3DPoints.begin();
    iter != _reconstructorData.map_3DPoints.end(); ++iter)
  {
    const size_t trackId = iter->first;
    const Vec3 & pt3D = iter->second;

    double maxAngle = 0.0;
    Vec3 originRay;
    // Look through the track and add point position
    const tracks::submapTrack & track = _map_reconstructed[trackId];
    for( tracks::submapTrack::const_iterator iterTrack = track.begin();
      iterTrack != track.end(); ++iterTrack)
    {
      const size_t imageId = iterTrack->first;
      const size_t featId = iterTrack->second;
      const BrownPinholeCamera & cam = _reconstructorData.map_Camera.find(imageId)->second;

      if ( set_camIndex.find(imageId) != set_camIndex.end())  {
        const std::vector<SIOPointFeature> & vec_feats = _map_feats[imageId];
        const SIOPointFeature & ptFeat = vec_feats[featId];

        double dResidual2D = cam.Residual(pt3D, ptFeat.coords().cast<double>());

        const Vec3 camPos = cam._C;
        const Vec3 dir = (pt3D - camPos).normalized();
        if (iterTrack == track.begin())
        {
          originRay = dir;
        }
        else
        {
          double dot = originRay.dot(dir);
          double angle = R2D(acos(clamp(dot, -1.0 + 1.e-8, 1.0 - 1.e-8)));
          maxAngle = max(angle, maxAngle);
        }

        // If residual too large, remove the measurement
        if (dResidual2D > dPrecision) {
          map_trackToErase[trackId].insert(imageId);
        }
      }
    }
    if (maxAngle < 3)
    {
      for( tracks::submapTrack::const_iterator iterTrack = track.begin();
        iterTrack != track.end(); ++iterTrack)  {
          const size_t imageId = iterTrack->first;
          map_trackToErase[trackId].insert(imageId);
      }
    }
  }

  size_t rejectedTrack = 0, rejectedMeasurement = 0;

  for (std::map<size_t, std::set<size_t> >::const_iterator iterT = map_trackToErase.begin();
    iterT != map_trackToErase.end(); ++iterT)
  {
    const size_t trackId = iterT->first;

    const std::set<size_t> setI = iterT->second;
    // Erase the image index reference
    for (std::set<size_t>::const_iterator iterTT = setI.begin();
      iterTT != setI.end(); ++iterTT, ++rejectedMeasurement) {
      _map_reconstructed[trackId].erase(*iterTT);
    }

    // If remaining tracks is too small, remove it
    if (_map_reconstructed[trackId].size() < 2) {
      _map_reconstructed[trackId].clear();
      _map_reconstructed.erase(trackId);
      _reconstructorData.set_trackId.erase(trackId);
      _reconstructorData.map_3DPoints.erase(trackId);
      ++rejectedTrack;
    }
  }

  std::cout << "\n#rejected track: " << set_trackToErase.size() << std::endl
    << "#rejected Entire track: " << rejectedTrack << std::endl
    << "#rejected Measurement: " << rejectedMeasurement << std::endl;
  return rejectedTrack + rejectedMeasurement;
}

void IncrementalReconstructionEngine::ColorizeTracks(std::vector<Vec3> & vec_tracksColor) const
{
  // Colorize each track
  //  Start with the most representative image
  //    and iterate to provide a color to each 3D point
  {
    C_Progress_display my_progress_bar(_map_reconstructed.size(),
                                       std::cout,
                                       "\nCompute scene structure color\n");

    vec_tracksColor.resize(_map_reconstructed.size());

    //Build a list of contiguous index for the trackIds
    std::map<size_t, size_t> trackIds_to_contiguousIndexes;
    size_t cpt = 0;
    for (openMVG::tracks::STLMAPTracks::const_iterator it = _map_reconstructed.begin();
      it != _map_reconstructed.end(); ++it, ++cpt)
    {
      trackIds_to_contiguousIndexes[it->first] = cpt;
    }

    // The track list that will be colored (point removed during the process)
    openMVG::tracks::STLMAPTracks mapTrackToColor(_map_reconstructed);
    while( !mapTrackToColor.empty() )
    {
      // Find the most representative image
      //  a. Count the number of visible point for each image
      //  b. Sort to find the most representative image

      std::map<size_t, size_t> map_IndexCardinal; // ImageIndex, Cardinal
      for (openMVG::tracks::STLMAPTracks::const_iterator
        iterT = mapTrackToColor.begin();
       iterT != mapTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = iterT->first;
        const tracks::submapTrack & track = mapTrackToColor[trackId];
        for( tracks::submapTrack::const_iterator iterTrack = track.begin();
          iterTrack != track.end(); ++iterTrack)
        {
          const size_t imageId = iterTrack->first;
          if (map_IndexCardinal.find(imageId) == map_IndexCardinal.end())
            map_IndexCardinal[imageId] = 1;
          else
            ++map_IndexCardinal[imageId];
        }
      }

      // Find the image that is the most represented
      std::vector<size_t> vec_cardinal;
      std::transform(map_IndexCardinal.begin(),
        map_IndexCardinal.end(),
        std::back_inserter(vec_cardinal),
        RetrieveValue());
      using namespace indexed_sort;
      std::vector< sort_index_packet_descend< size_t, size_t> > packet_vec(vec_cardinal.size());
      sort_index_helper(packet_vec, &vec_cardinal[0]);

      //First index is the image with the most of matches
      std::map<size_t, size_t>::const_iterator iterTT = map_IndexCardinal.begin();
      std::advance(iterTT, packet_vec[0].index);
      const size_t indexImage = iterTT->first;
      Image<RGBColor> image;
      ReadImage(
        stlplus::create_filespec(
          _sImagePath,
          stlplus::basename_part(_vec_camImageNames[indexImage].m_sImageName),
          stlplus::extension_part(_vec_camImageNames[indexImage].m_sImageName) ).c_str(), &image);

      // Iterate through the track
      std::set<size_t> set_toRemove;
      for (openMVG::tracks::STLMAPTracks::const_iterator
        iterT = mapTrackToColor.begin();
        iterT != mapTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = iterT->first;
        const tracks::submapTrack & track = mapTrackToColor[trackId];
        tracks::submapTrack::const_iterator iterF = track.find(indexImage);

        if (iterF != track.end())
        {
          // Color the track
          const size_t featId = iterF->second;
          const SIOPointFeature & feat = _map_feats.find(indexImage)->second[featId];
          RGBColor color = image(feat.y(), feat.x());

          vec_tracksColor[ trackIds_to_contiguousIndexes[trackId] ] = Vec3(color.r(), color.g(), color.b());
          set_toRemove.insert(trackId);
          ++my_progress_bar;
        }
      }
      // Remove colored track
      for (std::set<size_t>::const_iterator iter = set_toRemove.begin();
        iter != set_toRemove.end(); ++iter)
      {
        mapTrackToColor.erase(*iter);
      }
    }
  }
}

void IncrementalReconstructionEngine::BundleAdjustment()
{
  std::cout << std::endl
    << "---------------------------------\n"
    << "--      BUNDLE ADJUSTMENT      --\n"
    << "---------------------------------" << std::endl;

  //-- All the data that I must fill:
  using namespace std;

  const size_t nbCams = _reconstructorData.map_Camera.size();
  const size_t nbIntrinsics = _map_ImagesIdPerIntrinsicGroup.size();
  const size_t nbPoints3D = _reconstructorData.map_3DPoints.size();

  // Count the number of measurement (sum of the reconstructed track length)
  size_t nbmeasurements = 0;
  for (std::map<size_t,Vec3>::const_iterator iter = _reconstructorData.map_3DPoints.begin();
    iter != _reconstructorData.map_3DPoints.end();
    ++iter)
  {
    const size_t trackId = iter->first;
    // Look through the track and add point position
    const tracks::submapTrack & track = _map_reconstructed[trackId];
    nbmeasurements += track.size();
  }

  std::cout << "#Cams: " << nbCams << std::endl
    << "#Intrinsics: " << nbIntrinsics << std::endl
    << "#Points3D: " << nbPoints3D << std::endl
    << "#measurements: " << nbmeasurements << std::endl;

  // Setup a BA problem
  using namespace openMVG::bundle_adjustment;
  BA_Problem_data_camMotionAndIntrinsic<6,6> ba_problem;
  // Will refine extrinsics[R,t] per camera and grouped intrinsics [focal,ppx,ppy,k1,k2,k3].

    // Configure the size of the problem
  ba_problem.num_cameras_ = nbCams;
  ba_problem.num_intrinsics_ = nbIntrinsics;
  ba_problem.num_points_ = nbPoints3D;
  ba_problem.num_observations_ = nbmeasurements;

  ba_problem.point_index_.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_extrinsic.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_intrinsic.reserve(ba_problem.num_observations_);
  ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

  ba_problem.num_parameters_ =
    6 * ba_problem.num_cameras_ // #[Rotation|translation] = [3x1]|[3x1]
    + 6 * ba_problem.num_intrinsics_ // #[f,ppx,ppy,k1,k2,k3] = [6x1]
    + 3 * ba_problem.num_points_; // #[X] = [3x1]
  ba_problem.parameters_.reserve(ba_problem.num_parameters_);

  // Setup extrinsic parameters
  std::set<size_t> set_camIndex;
  std::map<size_t,size_t> map_camIndexToNumber_extrinsic, map_camIndexToNumber_intrinsic;
  size_t cpt = 0;
  for (reconstructorHelper::Map_BrownPinholeCamera::const_iterator iter = _reconstructorData.map_Camera.begin();
    iter != _reconstructorData.map_Camera.end();  ++iter, ++cpt)
  {
    // in order to map camera index to contiguous number
    set_camIndex.insert(iter->first);
    map_camIndexToNumber_extrinsic.insert(std::make_pair(iter->first, cpt));

    const Mat3 R = iter->second._R;
    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // translation
    const Vec3 t = iter->second._t;
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
  }

  // Setup intrinsic parameters groups
  cpt = 0;
  for (std::map<size_t, Vec6 >::const_iterator iterIntrinsicGroup = _map_IntrinsicsPerGroup.begin();
    iterIntrinsicGroup != _map_IntrinsicsPerGroup.end();
    ++iterIntrinsicGroup, ++cpt)
  {
    const Vec6 & intrinsic = iterIntrinsicGroup->second;
    ba_problem.parameters_.push_back(intrinsic(0)); // FOCAL_LENGTH
    ba_problem.parameters_.push_back(intrinsic(1)); // PRINCIPAL_POINT_X
    ba_problem.parameters_.push_back(intrinsic(2)); // PRINCIPAL_POINT_Y
    ba_problem.parameters_.push_back(intrinsic(3)); // K1
    ba_problem.parameters_.push_back(intrinsic(4)); // K2
    ba_problem.parameters_.push_back(intrinsic(5)); // K3
  }

  cpt = 0;
  // Link each camera to it's intrinsic group
  for (std::map<size_t, std::vector<size_t> >::const_iterator iterIntrinsicGroup = _map_ImagesIdPerIntrinsicGroup.begin();
    iterIntrinsicGroup != _map_ImagesIdPerIntrinsicGroup.end(); ++ iterIntrinsicGroup, ++cpt)
  {
    const std::vector<size_t> vec_imagesId = iterIntrinsicGroup->second;
    for (std::vector<size_t>::const_iterator iter_vec = vec_imagesId.begin();
      iter_vec != vec_imagesId.end(); ++iter_vec)
    {
      const size_t camIndex = *iter_vec;
      map_camIndexToNumber_intrinsic.insert(std::make_pair(camIndex, cpt));
    }
  }

  // Setup 3D points
  for (std::map<size_t,Vec3>::const_iterator iter = _reconstructorData.map_3DPoints.begin();
    iter != _reconstructorData.map_3DPoints.end();
    ++iter)
  {
    const Vec3 & pt3D = iter->second;
    ba_problem.parameters_.push_back(pt3D[0]);
    ba_problem.parameters_.push_back(pt3D[1]);
    ba_problem.parameters_.push_back(pt3D[2]);
  }

  // Fill measurements
  cpt = 0;
  for (std::map<size_t,Vec3>::const_iterator iter = _reconstructorData.map_3DPoints.begin();
    iter != _reconstructorData.map_3DPoints.end();
    ++iter)
  {
    const size_t trackId = iter->first;
    // Look through the track and add point position
    const tracks::submapTrack & track = _map_reconstructed[trackId];

    for( tracks::submapTrack::const_iterator iterTrack = track.begin();
      iterTrack != track.end();
      ++iterTrack)
    {
      const size_t imageId = iterTrack->first;
      const size_t featId = iterTrack->second;

      // If imageId reconstructed:
      //  - Add measurements (the feature position)
      //  - Add camidx (map the image number to the camera index)
      //  - Add ptidx (the 3D corresponding point index) (must be increasing)

      if ( set_camIndex.find(imageId) != set_camIndex.end())
      {
        const std::vector<SIOPointFeature> & vec_feats = _map_feats[imageId];
        const SIOPointFeature & ptFeat = vec_feats[featId];

        ba_problem.observations_.push_back( ptFeat.x() );
        ba_problem.observations_.push_back( ptFeat.y() );

        ba_problem.point_index_.push_back(cpt);
        ba_problem.camera_index_extrinsic.push_back(map_camIndexToNumber_extrinsic[imageId]);
        ba_problem.camera_index_intrinsic.push_back(map_camIndexToNumber_intrinsic[imageId]);
      }
    }
    ++cpt;
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(4.0);
  for (size_t i = 0; i < ba_problem.num_observations(); ++i) {
    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<pinhole_brown_reprojectionError::ErrorFunc_Refine_Camera_3DPoints, 2, 6, 6, 3>(
            new pinhole_brown_reprojectionError::ErrorFunc_Refine_Camera_3DPoints(
                & ba_problem.observations()[2 * i + 0]));

    problem.AddResidualBlock(cost_function,
                             p_LossFunction, // replaced by NULL if you don't want a LossFunction
                             ba_problem.mutable_camera_intrinsic_for_observation(i),
                             ba_problem.mutable_camera_extrinsic_for_observation(i),
                             ba_problem.mutable_point_for_observation(i));
  }

  //-- Lock the first camera to better deal with scene orientation ambiguity
  if (_vec_added_order.size()>0 && map_camIndexToNumber_extrinsic.size()>0)
  {
    // First camera is the first one that have been used
    problem.SetParameterBlockConstant(
      ba_problem.mutable_camera_extrinsic_for_observation(
        map_camIndexToNumber_extrinsic[_vec_added_order[0]]));
  }

  // Parametrization used to restrict camera intrinsics
  {
    //-- Optional:
    //  - bRefinePPandDisto
    //   -> true: Pinhole camera model with Brown distortion model may vary
    //   -> false: only the focal length may vary
    //
    //  - _bRefineFocal
    //   -> true: refine focal length, Principal point and radial distortion
    //   -> false: fixed focal length (refine Principal point and radial distortion)
    ceres::SubsetParameterization *constant_transform_parameterization = NULL;
    std::vector<int> vec_constant_intrinsic;
    if ( !_bRefinePPandDisto ) {
      // Parametrization used to set as constant principal point and radial distortion
      // Last five elements are ppx,ppy and radial disto factors.
      vec_constant_intrinsic.push_back(1); // PRINCIPAL_POINT_X FIXED
      vec_constant_intrinsic.push_back(2); // PRINCIPAL_POINT_Y FIXED
      vec_constant_intrinsic.push_back(3); // K1 FIXED
      vec_constant_intrinsic.push_back(4); // K2 FIXED
      vec_constant_intrinsic.push_back(5); // K3 FIXED
    }

    if ( !_bRefineFocal ) {
      // Parametrization used to set as constant the camera focal length (fixed focal length)
      vec_constant_intrinsic.push_back(0);
    }

    if ( !vec_constant_intrinsic.empty() &&
         vec_constant_intrinsic.size() != ba_problem.NINTRINSICPARAM
         // if all parameters are set as constant, better to use the SetParameterBlockConstant
       )
      constant_transform_parameterization =
        new ceres::SubsetParameterization(6, vec_constant_intrinsic);

    // Loop over intrinsics (to configure varying and fix parameters)
    for (size_t iIntrinsicGroupId = 0; iIntrinsicGroupId < ba_problem.num_intrinsics(); ++iIntrinsicGroupId)
    {
      if ( !_bRefinePPandDisto && !_bRefineFocal ){
        problem.SetParameterBlockConstant(ba_problem.mutable_cameras_intrinsic(iIntrinsicGroupId));
      }
      else{
        if ( !vec_constant_intrinsic.empty() ) {
          problem.SetParameterization(ba_problem.mutable_cameras_intrinsic(iIntrinsicGroupId),
            constant_transform_parameterization);
        }
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options options;
  options.preconditioner_type = ceres::JACOBI;
  options.linear_solver_type = ceres::SPARSE_SCHUR;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
  else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse back end for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;
#ifdef USE_OPENMP
  options.num_threads = omp_get_max_threads();
  options.num_linear_solver_threads = omp_get_max_threads();
#endif // USE_OPENMP

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (summary.IsSolutionUsable())
  {
    // Display statistics about the minimization
    std::cout << std::endl
      << "Bundle Adjustment statistics:\n"
      << " Initial RMSE: " << std::sqrt( summary.initial_cost / (ba_problem.num_observations_*2.)) << "\n"
      << " Final RMSE: " << std::sqrt( summary.final_cost / (ba_problem.num_observations_*2.)) << "\n"
      << std::endl;

    // Get back 3D points
    cpt = 0;
    for (std::map<size_t,Vec3>::iterator iter = _reconstructorData.map_3DPoints.begin();
      iter != _reconstructorData.map_3DPoints.end(); ++iter, ++cpt)
    {
      const double * pt = ba_problem.mutable_points() + cpt*3;
      Vec3 & pt3D = iter->second;
      pt3D = Vec3(pt[0], pt[1], pt[2]);
    }

    // Get back camera external and intrinsic parameters
    for (reconstructorHelper::Map_BrownPinholeCamera::iterator iter = _reconstructorData.map_Camera.begin();
      iter != _reconstructorData.map_Camera.end(); ++iter)
    {
      const size_t imageId = iter->first;
      const size_t extrinsicId = map_camIndexToNumber_extrinsic[imageId];
      // Get back extrinsic pointer
      const double * camE = ba_problem.mutable_cameras_extrinsic(extrinsicId);
      Mat3 R;
      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(camE, R.data());
      Vec3 t(camE[3], camE[4], camE[5]);

      // Get back the intrinsic group of the camera
      const size_t intrinsicId = map_camIndexToNumber_intrinsic[imageId];
      const double * camIntrinsics = ba_problem.mutable_cameras_intrinsic(intrinsicId);
      // Update the camera with update intrinsic and extrinsic parameters
      using namespace pinhole_brown_reprojectionError;
      BrownPinholeCamera & sCam = iter->second;
      sCam = BrownPinholeCamera(
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_FOCAL_LENGTH],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_PRINCIPAL_POINT_X],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_PRINCIPAL_POINT_Y],
        R,
        t,
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K1],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K2],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K3]);
    }

    //-- Update each intrinsic parameters group
    cpt = 0;
    for (std::map<size_t, Vec6 >::iterator iterIntrinsicGroup = _map_IntrinsicsPerGroup.begin();
      iterIntrinsicGroup != _map_IntrinsicsPerGroup.end();
      ++iterIntrinsicGroup, ++cpt)
    {
      const double * camIntrinsics = ba_problem.mutable_cameras_intrinsic() + cpt * 6;
      Vec6 & intrinsic = iterIntrinsicGroup->second;
      using namespace pinhole_brown_reprojectionError;
      intrinsic << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_FOCAL_LENGTH],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_PRINCIPAL_POINT_X],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_PRINCIPAL_POINT_Y],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K1],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K2],
        camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K3];

       std::cout << " for camera Idx=[" << cpt << "]: " << std::endl
        << "\t focal: " << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_FOCAL_LENGTH] << std::endl
        << "\t ppx: " << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_PRINCIPAL_POINT_X] << std::endl
        << "\t ppy: " << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_PRINCIPAL_POINT_Y] << std::endl
        << "\t k1: " << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K1] << std::endl
        << "\t k2: " << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K2] << std::endl
        << "\t k3: " << camIntrinsics[pinhole_brown_reprojectionError::OFFSET_K3] << std::endl;
    }
  }
}

double IncrementalReconstructionEngine::ComputeResidualsHistogram(Histogram<double> * histo)
{
  std::set<size_t> set_camIndex;
  for (reconstructorHelper::Map_BrownPinholeCamera::const_iterator iter = _reconstructorData.map_Camera.begin();
    iter != _reconstructorData.map_Camera.end();
    ++iter)
  {
    set_camIndex.insert(iter->first);
  }

  // For each 3D point sum their reprojection error

  std::vector<float> vec_residuals;
  vec_residuals.reserve(_reconstructorData.map_3DPoints.size());

  for (std::map<size_t,Vec3>::const_iterator iter = _reconstructorData.map_3DPoints.begin();
    iter != _reconstructorData.map_3DPoints.end();
    ++iter)
  {
    const size_t trackId = iter->first;
    const Vec3 & pt3D = iter->second;
    // Look through the track and add point position
    const tracks::submapTrack & track = _map_reconstructed[trackId];

    for( tracks::submapTrack::const_iterator iterTrack = track.begin();
      iterTrack != track.end();
      ++iterTrack)
    {
      const size_t imageId = iterTrack->first;
      const size_t featId = iterTrack->second;

      if ( set_camIndex.find(imageId) != set_camIndex.end())
      {
        const std::vector<SIOPointFeature> & vec_feats = _map_feats[imageId];
        const SIOPointFeature & ptFeat = vec_feats[featId];
        const BrownPinholeCamera & cam = _reconstructorData.map_Camera.find(imageId)->second;

        double dResidual = cam.Residual(pt3D, ptFeat.coords().cast<double>());
        vec_residuals.push_back(dResidual);
      }
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
      << "\t-- #Tracks:\t" << _map_reconstructed.size() << std::endl
      << "\t-- Residual min:\t" << dMin << std::endl
      << "\t-- Residual median:\t" << dMedian << std::endl
      << "\t-- Residual max:\t "  << dMax << std::endl
      << "\t-- Residual mean:\t " << dMean << std::endl;

    return dMean;
  }
  return -1.0;
}

} // namespace openMVG
