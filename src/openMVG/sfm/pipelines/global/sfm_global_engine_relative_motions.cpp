
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/multiview/essential.hpp"

#include "third_party/progress/progress.hpp"

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

GlobalSfMReconstructionEngine_RelativeMotions::GlobalSfMReconstructionEngine_RelativeMotions(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory), _sLoggingFile(sloggingFile), _normalized_features_provider(NULL) {

  if (!_sLoggingFile.empty())
  {
    // setup HTML logger
    _htmlDocStream = std::make_shared<htmlDocument::htmlDocumentStream>("GlobalReconstructionEngine SFM report.");
    _htmlDocStream->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("GlobalSfMReconstructionEngine_RelativeMotions")));
    _htmlDocStream->pushInfo("<hr>");

    _htmlDocStream->pushInfo( "Dataset info:");
    _htmlDocStream->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }

  // Set default motion Averaging methods
  _eRotationAveragingMethod = ROTATION_AVERAGING_L2;
  _eTranslationAveragingMethod = TRANSLATION_AVERAGING_L1;
}

GlobalSfMReconstructionEngine_RelativeMotions::~GlobalSfMReconstructionEngine_RelativeMotions()
{
  if (!_sLoggingFile.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(_sLoggingFile.c_str());
    htmlFileStream << _htmlDocStream->getDoc();
  }
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetFeaturesProvider(Features_Provider * provider)
{
  _features_provider = provider;

  // Copy features and save a normalized version
  _normalized_features_provider = std::make_shared<Features_Provider>(*provider);
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (Hash_Map<IndexT, PointFeatures>::iterator iter = _normalized_features_provider->feats_per_view.begin();
    iter != _normalized_features_provider->feats_per_view.end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif
    {
      // get the related view & camera intrinsic and compute the corresponding bearing vectors
      const View * view = _sfm_data.GetViews().at(iter->first).get();
      const std::shared_ptr<IntrinsicBase> cam = _sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
      for (PointFeatures::iterator iterPt = iter->second.begin();
        iterPt != iter->second.end(); ++iterPt)
      {
        const Vec3 bearingVector = (*cam)(cam->get_ud_pixel(iterPt->coords().cast<double>()));
        iterPt->coords() << (bearingVector.head(2) / bearingVector(2)).cast<float>();
      }
    }
  }
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetMatchesProvider(Matches_Provider * provider)
{
  _matches_provider = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetRotationAveragingMethod
(
  ERotationAveragingMethod eRotationAveragingMethod
)
{
  _eRotationAveragingMethod = eRotationAveragingMethod;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetTranslationAveragingMethod
(
  ETranslationAveragingMethod eTranslationAveragingMethod
)
{
  _eTranslationAveragingMethod = eTranslationAveragingMethod;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process() {

  //-------------------
  // Only keep the largest biedge connected subgraph
  //-------------------
  {
    const Pair_Set pairs = _matches_provider->getPairs();
    const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs, _sOutDirectory);
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
      return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, _matches_provider->_pairWise_matches);
  }

  Compute_Relative_Rotations(_relatives_Rt);

  if (!Compute_Global_Rotations())
  {
    std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Global_Translations())
  {
    std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
    return false;
  }
  if (!Compute_Initial_Structure())
  {
    std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
    return false;
  }
  if (!Adjust())
  {
    std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
    return false;
  }

  //-- Export statistics about the SfM process
  if (!_sLoggingFile.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion statistics.";
    _htmlDocStream->pushInfo("<hr>");
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- View count: " << _sfm_data.GetViews().size() << "<br>"
      << "-- Intrinsic count: " << _sfm_data.GetIntrinsics().size() << "<br>"
      << "-- Pose count: " << _sfm_data.GetPoses().size() << "<br>"
      << "-- Track count: "  << _sfm_data.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());
  }

  return true;
}

/// Compute from relative rotations the global rotations of the camera poses
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Rotations()
{
  if(_relatives_Rt.empty())
    return false;

  // Convert RelativeInfo_Map to appropriate input for solving the global rotations
  // - store the relative rotations and set a weight
  using namespace openMVG::rotation_averaging;
  RelativeRotations vec_relativeRotEstimate;
  {
    const PairWiseMatches & map_matches = _matches_provider->_pairWise_matches;
    std::set<IndexT> set_indeximage;
    for (PairWiseMatches::const_iterator iter = map_matches.begin(); iter != map_matches.end();  ++iter)
    {
      const IndexT I = iter->first.first;
      const IndexT J = iter->first.second;
      set_indeximage.insert(I);
      set_indeximage.insert(J);
    }

    std::cout << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << _relatives_Rt.size() << "\n"
      << "  #global rotations: " << set_indeximage.size() << std::endl;

    // Setup input relative rotation data for global rotation computation
    vec_relativeRotEstimate.reserve(_relatives_Rt.size());
    for(RelativeInfo_Map::const_iterator iter = _relatives_Rt.begin();
      iter != _relatives_Rt.end(); ++iter)
    {
      const openMVG::relativeInfo & rel = *iter;
      PairWiseMatches::const_iterator iterMatches = map_matches.find(rel.first); // If the pair support some matches
      if (iterMatches != map_matches.end())
      {
        vec_relativeRotEstimate.push_back(RelativeRotation(
          rel.first.first, rel.first.second,
          rel.second.first, 1.0f));
      }
    }

    //- weight computation: for each pair w = min(1, (#PairMatches/Median(#PairsMatches)))
    {
      //-- Compute the median number of matches
      std::vector<double> vec_count;
      for(RelativeInfo_Map::const_iterator iter = _relatives_Rt.begin();
        iter != _relatives_Rt.end(); ++iter)
      {
        const openMVG::relativeInfo & rel = *iter;
        // Find the number of support point for this pair
        PairWiseMatches::const_iterator iterMatches = map_matches.find(rel.first);
        if (iterMatches != map_matches.end())
        {
          vec_count.push_back(iterMatches->second.size());
        }
      }
      std::partial_sort(vec_count.begin(), vec_count.begin() + vec_count.size() / 2.0, vec_count.end());
      const float thTrustPair = vec_count[vec_count.size() / 2.0];

      for (RelativeRotations::iterator iterRelRot = vec_relativeRotEstimate.begin();
        iterRelRot != vec_relativeRotEstimate.end(); ++iterRelRot)
      {
        RelativeRotation & relRot = *iterRelRot;
        const Pair current_pair(relRot.i, relRot.j);
        float weight = 1.f; // the relative rotation weight
        PairWiseMatches::const_iterator iterMatches = map_matches.find(current_pair);
        if (iterMatches != map_matches.end())
        {
          weight = std::min((float)iterMatches->second.size()/thTrustPair, 1.f);
        }
        relRot.weight = weight;
      }
    }
  }

  // Global Rotation solver:
  ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod = TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;

  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  const bool bRotationAveraging = rotation_averaging_solver.Run(
    _eRotationAveragingMethod, eRelativeRotationInferenceMethod,
    vec_relativeRotEstimate, _map_globalR);

  if (bRotationAveraging)
  {
    // Log input graph to the HTML report
    if (!_sLoggingFile.empty() && !_sOutDirectory.empty())
    {
      // List the plausible remaining edges
      std::set<IndexT> set_ViewIds;
        std::transform(_sfm_data.GetViews().begin(), _sfm_data.GetViews().end(),
          std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      const std::string sGraph_name = "global_rotation_graph";
      graph::indexedGraph putativeGraph(set_ViewIds, rotation_averaging_solver.GetUsedPairs());
      graph::exportToGraphvizData(
        stlplus::create_filespec(_sOutDirectory, sGraph_name),
        putativeGraph.g);

      using namespace htmlDocument;
      std::ostringstream os;

      os << "<br>" << sGraph_name << "<br>"
         << "<img src=\""
         << stlplus::create_filespec(_sOutDirectory, sGraph_name, "svg")
         << "\" height=\"600\">\n";

      _htmlDocStream->pushInfo(os.str());
    }
  }
  return bRotationAveraging;
}

/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Translations()
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run(
    _eTranslationAveragingMethod,
    _sfm_data,
    _normalized_features_provider.get(),
    _matches_provider,
    _map_globalR,
    _tripletWise_matches);

  if (!_sLoggingFile.empty())
  {
    Save(_sfm_data,
      stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

/// Compute the initial structure of the scene
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Initial_Structure()
{
  // Build tracks from selected triplets (Union of all the validated triplet tracks (_tripletWise_matches))
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
    tracksBuilder.Build(_tripletWise_matches);
    tracksBuilder.Filter(3);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = _sfm_data.structure;
    IndexT idx(0);
    for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
      itTracks != map_selectedTracks.end();
      ++itTracks, ++idx)
    {
      const submapTrack & track = itTracks->second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
      {
        const size_t imaIndex = it->first;
        const size_t featIndex = it->second;
        const PointFeature & pt = _features_provider->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

    std::cout << std::endl << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats:
      //    - number of images
      //    - number of tracks
      std::set<size_t> set_imagesId;
      TracksUtilsMap::ImageIdInTracks(map_selectedTracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(map_selectedTracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)  {
        osTrack << "\t" << iter->first << "\t" << iter->second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Compute 3D position of the landmarks (structure) by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = _sfm_data.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(_sfm_data);

    std::cout << "\n#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(_sfm_data.GetLandmarks().size()) << std::endl;
    std::cout << std::endl << "  Triangulation took (s): " << timer.elapsed() << std::endl;

    // Export initial structure
    if (!_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !_sfm_data.structure.empty();
}

// Adjust the scene (& remove outliers)
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, false, true, false);
  if (b_BA_Status)
  {
    if (!_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, false);
    if (b_BA_Status && !_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && !_bFixedIntrinsics) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, true);
    if (b_BA_Status && !_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = _sfm_data.structure.size();
  RemoveOutliers_PixelResidualError(_sfm_data, 4.0);
  const size_t pointcount_pixelresidual_filter = _sfm_data.structure.size();
  RemoveOutliers_AngleError(_sfm_data, 2.0);
  const size_t pointcount_angular_filter = _sfm_data.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!_sLoggingFile.empty())
  {
    Save(_sfm_data,
      stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(_sfm_data, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = _sfm_data.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";

    b_BA_Status = bundle_adjustment_obj.Adjust(_sfm_data, true, true, !_bFixedIntrinsics);
    if (b_BA_Status && !_sLoggingFile.empty())
    {
      Save(_sfm_data,
        stlplus::create_filespec(stlplus::folder_part(_sLoggingFile), "structure_04_unstable_removed", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  return b_BA_Status;
}

void GlobalSfMReconstructionEngine_RelativeMotions::Compute_Relative_Rotations(RelativeInfo_Map & vec_relatives)
{
  // For each pair of view of the matching graph, compute the rotation from pairwise point matches:
  const Pair_Set & pair_set = _matches_provider->getPairs();
  // copy to a vector for use with threading
  const Pair_Vec pair_vec(pair_set.begin(), pair_set.end());
  C_Progress_display my_progress_bar( pair_vec.size(), std::cout, "\nCompute_Relative_Rotations\n" );
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < pair_vec.size(); ++i)
  {
    ++my_progress_bar;

    const Pair current_pair = pair_vec[i];
    const size_t I = current_pair.first;
    const size_t J = current_pair.second;

    const View * view_I = _sfm_data.views[I].get();
    const View * view_J = _sfm_data.views[J].get();

    // Check that valid camera are existing for the pair index
    if (_sfm_data.GetIntrinsics().find(view_I->id_intrinsic) == _sfm_data.GetIntrinsics().end() ||
      _sfm_data.GetIntrinsics().find(view_J->id_intrinsic) == _sfm_data.GetIntrinsics().end())
      continue;

    const IndMatches & vec_matchesInd = _matches_provider->_pairWise_matches.find(current_pair)->second;

    Mat x1(2, vec_matchesInd.size()), x2(2, vec_matchesInd.size());
    for (size_t k = 0; k < vec_matchesInd.size(); ++k)
    {
      x1.col(k) = _normalized_features_provider->feats_per_view[I][vec_matchesInd[k]._i].coords().cast<double>();
      x2.col(k) = _normalized_features_provider->feats_per_view[J][vec_matchesInd[k]._j].coords().cast<double>();
    }

    const IntrinsicBase * cam_I = _sfm_data.GetIntrinsics().at(view_I->id_intrinsic).get();
    const IntrinsicBase * cam_J = _sfm_data.GetIntrinsics().at(view_J->id_intrinsic).get();
    if ( !isValid(cam_I->getType()) || !isValid(cam_J->getType()))
    {
      continue;
    }

    // Compute relative poses from the view graph (thanks to a robust essential matrix estimation):

    RelativePose_Info relativePose_info;
    // Compute max authorized error as geometric mean of camera plane tolerated residual error
    relativePose_info.initial_residual_tolerance = std::pow(
      cam_I->imagePlane_toCameraPlaneError(2.5) *
      cam_J->imagePlane_toCameraPlaneError(2.5),
      1./2.);
    // Since we use normalized features:
    const std::pair<size_t, size_t> imageSize_I(1., 1.), imageSize_J(1.,1.);
    const Mat3 K  = Mat3::Identity();

    if (!robustRelativePose(K, K, x1, x2, relativePose_info, imageSize_I, imageSize_J, 256))
    {
      std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
        << std::endl;
      continue;
    }
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
      for (size_t k = 0; k < x1.cols(); ++k) {
        const Vec2 x1_ = _features_provider->feats_per_view[I][vec_matchesInd[k]._i].coords().cast<double>();
        const Vec2 x2_ = _features_provider->feats_per_view[J][vec_matchesInd[k]._j].coords().cast<double>();
        Vec3 X;
        TriangulateDLT(P1, x1_, P2, x2_, &X);
        Observations obs;
        obs[view_I->id_view] = Observation(x1_, vec_matchesInd[k]._i);
        obs[view_J->id_view] = Observation(x2_, vec_matchesInd[k]._j);
        landmarks[k].obs = obs;
        landmarks[k].X = X;
      }
      // - refine only Structure and Rotations & translations (keep intrinsic constant)
      Bundle_Adjustment_Ceres::BA_options options(false, false);
      options._linear_solver_type = ceres::DENSE_SCHUR;
      Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
      if (bundle_adjustment_obj.Adjust(tiny_scene, true, true, false))
      {
        // --> to debug: save relative pair geometry on disk
        // std::ostringstream os;
        // os << current_pair.first << "_" << current_pair.second << ".ply";
        // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
        //
        const Mat3 R1 = tiny_scene.poses[view_I->id_pose].rotation();
        const Mat3 R2 = tiny_scene.poses[view_J->id_pose].rotation();
        const Vec3 t1 = tiny_scene.poses[view_I->id_pose].translation();
        const Vec3 t2 = tiny_scene.poses[view_J->id_pose].translation();
        // Compute relative motion and save it
        Mat3 Rrel;
        Vec3 trel;
        RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
        // Update found relative pose
        relativePose_info.relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
      }
    }
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
    {
      vec_relatives[current_pair] = std::make_pair(
        relativePose_info.relativePose.rotation(),
        relativePose_info.relativePose.translation());

      // TODO: Extend later to a pose graph and not a view graph
      // convert relative view index to relative pose indexes
      //const Pair relativePoseIndexes(view_I->id_pose, view_J->id_pose);
      //vec_relatives[relativePoseIndexes] = std::make_pair(
      //  relativePose_info.relativePose.rotation(),
      //  relativePose_info.relativePose.translation());
    }
  }
  // Log input graph to the HTML report
  if (!_sLoggingFile.empty() && !_sOutDirectory.empty())
  {
    std::set<IndexT> set_ViewIds;
      std::transform(_sfm_data.GetViews().begin(), _sfm_data.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
    graph::indexedGraph putativeGraph(set_ViewIds, getPairs(_matches_provider->_pairWise_matches));
    graph::exportToGraphvizData(
      stlplus::create_filespec(_sOutDirectory, "input_largest_cc_relative_motions_graph"),
      putativeGraph.g);

    using namespace htmlDocument;
    std::ostringstream os;

    os << "<br>" << "input_largest_cc_relative_motions_graph" << "<br>"
       << "<img src=\""
       << stlplus::create_filespec(_sOutDirectory, "input_largest_cc_relative_motions_graph", "svg")
       << "\" height=\"600\">\n";

    _htmlDocStream->pushInfo(os.str());
  }
}

} // namespace sfm
} // namespace openMVG

