// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/stellar/sfm_stellar_engine.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/triplet_finder.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTij.hpp"

#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/relative_pose_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/pipelines/stellar/stellar_definitions.hpp"
#include "openMVG/sfm/pipelines/stellar/stellar_solver.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include <ceres/types.h>

#include "lemon/kruskal.h"

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

StellarSfMReconstructionEngine::StellarSfMReconstructionEngine(
  const SfM_Data & sfm_data,
  const std::string & out_directory_logging,
  const std::string & logging_file)
  : ReconstructionEngine(sfm_data, out_directory_logging),
    html_doc_stream_(nullptr),
    logging_file_(logging_file),
    features_provider_(nullptr),
    matches_provider_(nullptr),
    graph_simplification_(EGraphSimplification::MST_X),
    graph_simplification_value_(5)
{
  if (!logging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("Stellar SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("StellarSfMReconstructionEngine")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }
}

StellarSfMReconstructionEngine::~StellarSfMReconstructionEngine()
{
  if (!logging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(logging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void StellarSfMReconstructionEngine::SetFeaturesProvider(
  Features_Provider * provider)
{
  features_provider_ = provider;
}

void StellarSfMReconstructionEngine::SetMatchesProvider(
  Matches_Provider * provider)
{
  matches_provider_ = provider;
}

void StellarSfMReconstructionEngine::SetGraphSimplification
(
  const EGraphSimplification graph_simplification,
  const int value
)
{
  graph_simplification_ = graph_simplification;
  graph_simplification_value_ = value;
}

bool StellarSfMReconstructionEngine::Process()
{
  // Compute a relative pose for each pair
  Hash_Map<Pair, geometry::Pose3> relative_poses;
  ComputeRelativeMotions(relative_poses);

  // From the found relative motion, compute the stellar/nuplet reconstructions
  Hash_Map<IndexT, StellarPodRelativeMotions> stellar_reconstruction_per_pose;
  if (!ComputeStellarReconstructions(relative_poses, stellar_reconstruction_per_pose))
    return false;

  OPENMVG_LOG_INFO << "#Stellar reconstruction pod: " << stellar_reconstruction_per_pose.size();

  // Perform the rotation averaging and compute the global rotations
  //a- extract largest CC first
  openMVG::rotation_averaging::RelativeRotations relatives_R;
  Pair_Set pairs;
  std::set<IndexT> base_nodes;
  for (const auto & pod_it : stellar_reconstruction_per_pose)
  {
    for (const auto & rel_it : pod_it.second.relative_motions)
    {
      pairs.insert({rel_it.first.first, rel_it.first.second});
      base_nodes.insert(rel_it.first.first);
      base_nodes.insert(rel_it.first.second);
    }
  }

  const auto largest_cc_nodes = graph::KeepLargestCC_Nodes<Pair_Set, IndexT>(pairs);
  OPENMVG_LOG_INFO << "Largest CC => Keeping " << largest_cc_nodes.size() << " nodes from " << base_nodes.size() << " base nodes";
  
  for (const auto & pod_it : stellar_reconstruction_per_pose)
  {
    for (const auto & rel_it : pod_it.second.relative_motions)
    {
      if(largest_cc_nodes.count(rel_it.first.first) && largest_cc_nodes.count(rel_it.first.second)) {
        relatives_R.emplace_back(
          rel_it.first.first, rel_it.first.second, // {i,j}
          rel_it.second.rotation, // R_ij, the relative rotation
          1.0); // Weight
      }
    }
  }
  
  //b- the rotation averaging step 
  Hash_Map<IndexT, Mat3> global_rotations;
  Pair_Set selected_pairs;
  if (!Compute_Global_Rotations(relatives_R, global_rotations, selected_pairs))
  {
    OPENMVG_LOG_ERROR << "GlobalSfM:: Rotation Averaging failure!";
    return false;
  }

  
  // Perform the translation averaging and compute the global positions
  //a- be sure rotation averaging has not disconnected the graph
  const auto largest_cc_nodes_rot = graph::KeepLargestCC_Nodes<Pair_Set, IndexT>(selected_pairs);
  assert(largest_cc_nodes_rot.size() == global_rotations.size());
  
  //b- translation averaging sterp
  if (!Compute_Global_Translations(global_rotations, selected_pairs, stellar_reconstruction_per_pose))
  {
    OPENMVG_LOG_ERROR << "GlobalSfM:: Translation Averaging failure!";
    return false;
  }

  {
    const int min_covisibility = 2;
    if (!Compute_Initial_Structure(min_covisibility))
    {
      OPENMVG_LOG_ERROR << "GlobalSfM:: Cannot initialize an initial structure!";
      return false;
    }

    if (!Adjust())
    {
      OPENMVG_LOG_ERROR << "GlobalSfM:: Non-linear adjustment failure!";
      return false;
    }
  }
  return true;
}

Pair_Set selectMST(const openMVG::matching::PairWiseMatches & matches,
               const int nbTree)
{
  Pair_Set selected_pairs;

  // Build a graph, affect edge weight depending of supporting match
  // Do N time : select an MST (set weight to 0)
  // return the supporting edges

  graph::indexedGraph graph_builder(getPairs(matches));
  using map_EdgeMap = graph::indexedGraph::GraphT::EdgeMap<double>;
  map_EdgeMap edge_map(graph_builder.g);

  using EdgeIterator = graph::indexedGraph::GraphT::EdgeIt;
  for ( EdgeIterator it( graph_builder.g ); it != lemon::INVALID; ++it )
  {
    const auto & found = matches.at(std::make_pair(
      (*graph_builder.node_map_id)[graph_builder.g.u(it)],
      (*graph_builder.node_map_id)[graph_builder.g.v(it)]));
    edge_map[it] = - static_cast<double>(found.size());
  }

  int nb_select = nbTree;
  // Select Many MST
  while (nb_select>0)
  {
    std::vector<graph::indexedGraph::GraphT::Edge> tree_edge_vec;
    lemon::kruskal(graph_builder.g, edge_map, std::back_inserter(tree_edge_vec));

    //E-- Export compute MST
    for (const auto & edge_it : tree_edge_vec )
    {
      selected_pairs.insert(
       {(*graph_builder.node_map_id)[graph_builder.g.u(edge_it)],
       (*graph_builder.node_map_id)[graph_builder.g.v(edge_it)]});
       // Update the weight to avoid resampling the same tree next time
      edge_map[edge_it] = 0;
    }
    --nb_select;
  }
  OPENMVG_LOG_INFO << "Graph simplification:\n"
    <<"\tUsing #trees: " << nbTree << " selected #pairs: " << selected_pairs.size() << "\n"
    <<"\tOriginal graph #pairs: " << matches.size() << "\n"
    <<"\tKeep: " << selected_pairs.size() / (float)matches.size() << " % of original graph.";
  return selected_pairs;
}

void StellarSfMReconstructionEngine::ComputeRelativeMotions(
  Hash_Map<Pair, Pose3> & relative_poses
) const
{
  /// If there is a JSON cache use it:
  /*if (stlplus::is_file(sOut_directory_ + "/relative_motion_cache"))
  {
    std::ifstream stream((sOut_directory_ + "/relative_motion_cache").c_str());

    cereal::JSONInputArchive archive(stream);
    archive(relative_poses);

    std::cout
      << "\n::Compute_Relative_Motions\n"
      << "Loaded #RelativePoses: " << relative_poses.size() << std::endl;
    return;
  }*/

  Pair_Set selected_pairs;
  switch (graph_simplification_)
  {
  case EGraphSimplification::NONE:
  { // No simplification
    selected_pairs = matching::getPairs(matches_provider_->pairWise_matches_);
  }
  break;
  case EGraphSimplification::STAR_X: // Keep the X edges per n-Uplet
  {
    // Limit the complexity of each sub-star
    int max_stellar_pod_size = graph_simplification_value_; // 10

    const Pair_Set pairs = matching::getPairs(matches_provider_->pairWise_matches_);

    // Enhance the selected pairs in order to ensure that each NUplets have at maximum max_stellar_pod_size edges
    {
      // List all stellar configurations
      using StellarPods = Hash_Map<IndexT, Pair_Set>;
      StellarPods stellar_pods;
      for (const auto & it : pairs)
      {
        const Pair & pairIt = it;
        stellar_pods[pairIt.first].insert(pairIt);
        stellar_pods[pairIt.second].insert(pairIt);
      }
      // Keep only the Nth best pairs for each node
      for (auto & stellar_pod_it : stellar_pods)
      {
        Pair_Set & pairs = stellar_pod_it.second;
        unsigned int count = 0;
        std::vector<std::pair<unsigned int, unsigned int> > matches_per_edges;
        for (const Pair & pair_it : pairs)
        {
          matches_per_edges.emplace_back(matches_provider_->pairWise_matches_.at(pair_it).size(), count++);
        }
        std::sort(matches_per_edges.begin(), matches_per_edges.end());
        Pair_Set cpy;
        for ( std::vector<std::pair<unsigned int, unsigned int> >::const_reverse_iterator iter = matches_per_edges.rbegin();
          iter != matches_per_edges.rend() && cpy.size() < max_stellar_pod_size; ++ iter)
        {
          const unsigned int offset = iter->second;
          Pair_Set::const_iterator iterP = pairs.begin();
          std::advance(iterP, offset);
          cpy.insert(*iterP);
          selected_pairs.insert(*iterP);
        }
      }
    }
    // In order to guarantee the same connectivity as the initial graph, we add a MST to the pairs
    const int nb_trees = 1;
    const auto mst_pairs = selectMST(matches_provider_->pairWise_matches_, nb_trees);
    selected_pairs.insert(mst_pairs.begin(), mst_pairs.end());
  }
  break;
  case EGraphSimplification::MST_X:
  {
      // Select X random MST
      const int nb_trees = graph_simplification_value_;
      //Note: nb_trees could be dynamic depending of the node degree (median)
      selected_pairs = selectMST(matches_provider_->pairWise_matches_, nb_trees);
  }
  break;
  
  default:
    OPENMVG_LOG_ERROR << "Unknown graph simplification method";
    break;
  }

  OPENMVG_LOG_INFO << "Will use: " << selected_pairs.size() << " relative pose pairs.";
  OPENMVG_LOG_INFO << "Initial graph got: " << matches_provider_->pairWise_matches_.size() << " view pairs";

  // Compute a relative pose for each edge of the pose pair graph
  relative_poses = [&]
  {
    Relative_Pose_Engine relative_pose_engine;
    if (!relative_pose_engine.Process(
        selected_pairs,
        sfm_data_,
        matches_provider_,
        features_provider_))
      return Relative_Pose_Engine::Relative_Pair_Poses();
    else
      return relative_pose_engine.Get_Relative_Poses();
  }();

  /// Export to a JSON cache (for the relative motions):
  /*
  {
    std::ofstream stream((sOut_directory_ + "/relative_motion_cache").c_str());

    cereal::JSONOutputArchive archive(stream);
    archive(relative_poses);
  }
  */
}

bool StellarSfMReconstructionEngine::ComputeStellarReconstructions
(
  const Hash_Map<Pair, geometry::Pose3> & relative_poses,
  Hash_Map<IndexT, StellarPodRelativeMotions > & stellar_reconstruction_per_pose
) const
{
  OPENMVG_LOG_INFO << "::Compute_Stellar_Reconstructions";

  // List all stellar configurations
  using StellarPods = Hash_Map<IndexT, Pair_Set>;
  StellarPods stellar_pods;
  for (const auto & it : relative_poses)
  {
    const Pair & pairIt = it.first;
    stellar_pods[pairIt.first].insert(pairIt);
    stellar_pods[pairIt.second].insert(pairIt);
  }

  // Display some debug information about the stellar pods
  {
    OPENMVG_LOG_INFO << "Stellar debug: \n"
      << "#Poses: " << stellar_pods.size();

    for (const auto & stellar_pod_it : stellar_pods)
    {
      const IndexT node_id = stellar_pod_it.first;
      const Pair_Set & pairs = stellar_pod_it.second;
      OPENMVG_LOG_INFO << node_id << " => #pairs: " << pairs.size();
    }
  }

  // Compute the local reconstruction of the pods and save its relative pose motions

  system::LoggerProgress my_progress_bar( stellar_pods.size(), "- Stellar pod reconstruction -" );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int stellar_it = 0; stellar_it < static_cast<int>(stellar_pods.size()); ++stellar_it)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
    ++my_progress_bar;

    auto stellar_pod_it = stellar_pods.cbegin();
    std::advance(stellar_pod_it, stellar_it);
    const IndexT node_id = stellar_pod_it->first;
    const Pair_Set & pairs = stellar_pod_it->second;

    const bool use_all_matches = false;
    const bool use_threading = false;
    Stellar_Solver stellar_pod_solver(
      pairs,
      relative_poses,
      sfm_data_,
      matches_provider_,
      features_provider_,
      use_all_matches,
      use_threading);

    Poses poses;
    if (!stellar_pod_solver.Solve(poses))
    {
      OPENMVG_LOG_WARNING << "Cannot solve this Stellar pod";
      continue;
    }

#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
    {
      // Update the relative motions for this pod
      for (const auto pair_it : pairs)
      {
        const int I = pair_it.first;
        const int J = pair_it.second;
        Mat3 R_refined;
        Vec3 t_refined;
        RelativeCameraMotion(
          poses[I].rotation(), poses[I].translation(),
          poses[J].rotation(), poses[J].translation(),
          &R_refined, &t_refined);

        stellar_reconstruction_per_pose[node_id].relative_motions.insert(
          {pair_it, {R_refined, t_refined}});
      }
    }
  }

  return true;
}

bool StellarSfMReconstructionEngine::Compute_Global_Rotations
(
  const rotation_averaging::RelativeRotations & relatives_R,
  Hash_Map<IndexT, Mat3> & global_rotations,
  Pair_Set & used_pairs
) const
{
  if (relatives_R.empty())
    return false;
  // Log statistics about the relative rotation graph
  {
    std::set<IndexT> set_pose_ids;
    for (const auto & relative_R : relatives_R)
    {
      set_pose_ids.insert(relative_R.i);
      set_pose_ids.insert(relative_R.j);
    }

    OPENMVG_LOG_INFO << "-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << relatives_R.size() << "\n"
      << "  #global rotations: " << set_pose_ids.size();
  }

  // Global Rotation solver:
  const ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod =
    //TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    TRIPLET_ROTATION_INFERENCE_NONE;

  system::Timer timer;
  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  const bool b_rotation_averaging =
    rotation_averaging_solver.Run(
      ROTATION_AVERAGING_L2, eRelativeRotationInferenceMethod,
      relatives_R, global_rotations);

  OPENMVG_LOG_INFO << "Found #global_rotations: " << global_rotations.size() << "\n"
    << "in " << timer.elapsed() << " seconds.";

  if (b_rotation_averaging)
  {
    used_pairs = rotation_averaging_solver.GetUsedPairs();
    // Compute fitting error
    std::vector<float> vec_rotation_fitting_error;
    vec_rotation_fitting_error.reserve(relatives_R.size());
    for (const auto & relative_R : relatives_R)
    {
      const Mat3 & Rij = relative_R.Rij;
      const IndexT i = relative_R.i;
      const IndexT j = relative_R.j;
      if (global_rotations.count(i)==0 || global_rotations.count(j)==0)
        continue;
      const Mat3 & Ri = global_rotations[i];
      const Mat3 & Rj = global_rotations[j];
      const Mat3 eRij(Rj.transpose()*Rij*Ri);
      const double angularErrorDegree = R2D(getRotationMagnitude(eRij));
      vec_rotation_fitting_error.push_back(angularErrorDegree);
    }

    // TODO list inliers/outlier with a boolean array
    // Export those pairs to see the problem
    // Check if a pair if marked many time (since it can be listed many time)

    // Display some statistics about the relative to global rotation fitting estimation
    if (!vec_rotation_fitting_error.empty())
    {
      Histogram<float> histo(0.0f, *max_element(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend()), 20);
      histo.Add(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend());
      OPENMVG_LOG_INFO
        << "Relative rotations fitting error to global rotations:"
        << histo.ToString();
      {
        Histogram<float> histo(0.0f, 5.0f, 20);
        histo.Add(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend());
        OPENMVG_LOG_INFO
          << "Relative rotations fitting error to global rotations {0,5}:"
          << histo.ToString();
      }
      OPENMVG_LOG_INFO << "Statistics about global rotation fitting:";

      minMaxMeanMedian<float>(vec_rotation_fitting_error.cbegin(), vec_rotation_fitting_error.cend(), std::cout);
    }

    // Log input graph to the HTML report
    if (!logging_file_.empty() && !sOut_directory_.empty())
    {
      // Log a relative pose graph
      {
        std::set<IndexT> set_pose_ids;
        Pair_Set relative_pose_pairs;
        for (const auto & view : sfm_data_.GetViews())
        {
          const IndexT pose_id = view.second->id_pose;
          set_pose_ids.insert(pose_id);
        }
        const std::string sGraph_name = "global_relative_rotation_pose_graph_final";
        graph::indexedGraph putativeGraph(set_pose_ids, rotation_averaging_solver.GetUsedPairs());
        graph::exportToGraphvizData(
          stlplus::create_filespec(sOut_directory_, sGraph_name),
          putativeGraph);

        using namespace htmlDocument;
        std::ostringstream os;

        os << "<br>" << sGraph_name << "<br>"
           << "<img src=\""
           << stlplus::create_filespec(sOut_directory_, sGraph_name, "svg")
           << "\" height=\"600\">\n";

        html_doc_stream_->pushInfo(os.str());
      }
    }
  }

  return b_rotation_averaging;
}

/// Compute global translations
bool StellarSfMReconstructionEngine::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  const Pair_Set & used_pairs,
  const Hash_Map<IndexT, StellarPodRelativeMotions> & stellar_reconstruction_per_pose
)
{
  // Collect the relative motions between valid image pairs:
  //  pairs belonging to a stellar reconstruction and for which global rotations has been estimated
  Pair_Set valid_pairs;
  std::vector<RelativeInfo_Vec> stellar_relative_motions;
  std::set<IndexT> poses_ids;
  IndexT relative_motion_count = 0;
  {
    for (const auto & stellar_pod_it : stellar_reconstruction_per_pose)
    {
      const IndexT node_id = stellar_pod_it.first;
      if (global_rotations.count(node_id) == 0)
        continue;
      RelativeInfo_Vec relative_motion_infos;
      for (const auto & relative_info_it : stellar_pod_it.second.relative_motions)
      {
        const IndexT i = relative_info_it.first.first;
        const IndexT j = relative_info_it.first.second;

        if (global_rotations.count(i)==0 || global_rotations.count(j)==0 || used_pairs.count(relative_info_it.first) == 0) 
          continue;
        
        valid_pairs.insert({i,j});
        poses_ids.insert(i);
        poses_ids.insert(j);

        // update the relative rotation to use global rotations
        const Mat3 & Ri = global_rotations.at(i);
        const Mat3 & Rj = global_rotations.at(j);
        const Mat3 relative_rotation = Rj * Ri.transpose();

        relative_motion_infos.emplace_back(
          relative_info_it.first,
          std::make_pair(relative_rotation, relative_info_it.second.translation));
      }
      if (!relative_motion_infos.empty())
      {
        relative_motion_count += relative_motion_infos.size();
        stellar_relative_motions.push_back(relative_motion_infos);
      }
    }
  }

  const size_t iNview = poses_ids.size();
  OPENMVG_LOG_INFO << "-------------------------------" << "\n"
    << " Global translations computation: " << "\n"
    << "   - Ready to compute " << iNview << " global translations." << "\n"
    << "     from #relative translations: " << relative_motion_count;

  openMVG::system::Timer timer_translation;

  // Re-index the pair pose indices for solving the translations
  // Build a mapping from the native index to a contiguous sequence
  Hash_Map<IndexT, IndexT> reindex_forward, reindex_backward;
  reindex(valid_pairs, reindex_forward, reindex_backward);
  for (auto & stellar_group_it : stellar_relative_motions)
  {
    for (auto & relative_motion_it : stellar_group_it)
    {
      const auto & i = relative_motion_it.first.first;
      const auto & j = relative_motion_it.first.second;
      // reindex the pose ids
      relative_motion_it.first = {reindex_forward[i], reindex_forward[j]};
    }
  }

  std::vector<Vec3> vec_translations;
  if (!solve_translations_problem_softl1(stellar_relative_motions, vec_translations))
  {
    OPENMVG_LOG_ERROR << "Compute global translations: failed";
    return false;
  }

  //-- Export statistics:
  {

    std::ostringstream os;
    os << "Translation averaging statistics.";
    os.str("");
    os << "-------------------------------" << "\n"
      << "-- #relative estimates: " << relative_motion_count << "\n"
      << "-- #global translations: " <<  vec_translations.size() << "\n"
      << " timing (s): " << timer_translation << ".\n"
      << "-------------------------------" << "\n";
    OPENMVG_LOG_INFO << os.str();
  }

  // A valid solution was found:
  // - Update the view poses according the found camera translations
  for (size_t i = 0; i < vec_translations.size(); ++i)
  {
    const Vec3 & t = vec_translations[i];
    const IndexT pose_id = reindex_backward.at(i);
    const Mat3 & Ri = global_rotations.at(pose_id);
    sfm_data_.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
  }

  if (!logging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(logging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return true;
}

/// Compute the initial structure of the scene
bool StellarSfMReconstructionEngine::Compute_Initial_Structure(const int min_covisibility)
{
  // TODO => Remove the pairwise matches that does not belong to valid poses IDS!!

  // Build tracks
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
    tracksBuilder.Build(matches_provider_->pairWise_matches_);
    tracksBuilder.Filter(min_covisibility);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = sfm_data_.structure;
    IndexT idx(0);
    for (const auto & tracks_it : map_selectedTracks)
    {
      const submapTrack & track = tracks_it.second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (const auto & track_it : track)
      {
        const size_t ima_index = track_it.first;
        const size_t feat_index = track_it.second;
        const PointFeature & pt = features_provider_->feats_per_view.at(ima_index)[feat_index];
        obs[ima_index] = std::move(Observation(pt.coords().cast<double>(), feat_index));
      }
      ++idx;
    }
  }

  // Ensure tracks are visible by some poses
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 2; // 2 min
  eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength);

  // Compute 3D position of the landmark of the structure by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = sfm_data_.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(sfm_data_);
    //
    // Remove unstable depth points
    DepthCleaning(sfm_data_);

    OPENMVG_LOG_INFO << "#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(sfm_data_.GetLandmarks().size());
    OPENMVG_LOG_INFO << "  Triangulation took (s): " << timer.elapsed();

    // Export initial structure
    if (!logging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(logging_file_), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !sfm_data_.structure.empty();
}

bool StellarSfMReconstructionEngine::Adjust()
{
  if (ReconstructionEngine::intrinsic_refinement_options_ != cameras::Intrinsic_Parameter_Type::NONE)
  {
    Bundle_Adjustment_Ceres::BA_Ceres_options options;
    if ( sfm_data_.GetPoses().size() > 100 &&
        (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
         ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
         ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
        )
    // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
    {
      options.preconditioner_type_ = ceres::JACOBI;
      options.linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
    else
    {
      options.linear_solver_type_ = ceres::DENSE_SCHUR;
    }
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    const Optimize_Options ba_refine_options
      ( ReconstructionEngine::intrinsic_refinement_options_,
        Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
        Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
        Control_Point_Parameter(),
        this->b_use_motion_prior_
      );
    bool b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
    if (b_BA_Status)
    {
      DepthCleaning(sfm_data_);
    }
    if (b_BA_Status && !logging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(logging_file_), "structure_ba", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

      // Remove outliers (max_angle, residual error)
      // Remove unstable triangulations and camera poses
      RemoveOutliers_AngleError(sfm_data_, 2.0);
      RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
      eraseUnstablePosesAndObservations(sfm_data_);

      b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);

      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(logging_file_), "structure_ba_outlier_removed", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));

      return b_BA_Status;
    }
  }
  return true;
}

} // namespace sfm
} // namespace openMVG
