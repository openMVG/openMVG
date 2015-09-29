
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/global/mutexSet.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/sfm/pipelines/global/triplet_t_ACRansac_kernelAdaptator.hpp"
#include "third_party/histogram/histogram.hpp"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;

/// Use features in normalized camera frames
bool GlobalSfM_Translation_AveragingSolver::Run(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  SfM_Data & sfm_data,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Compute the relative translations and save them to vec_initialRijTijEstimates:
  Compute_translations(
    sfm_data,
    normalized_features_provider,
    matches_provider,
    map_globalR,
    tripletWise_matches);

  const bool b_translation = Translation_averaging(
    eTranslationAveragingMethod,
    sfm_data,
    map_globalR);

  // Filter matches to keep only them link to view that have valid poses
  // (necessary since multiple components exists before translation averaging)
  std::set<IndexT> valid_view_ids;
  for (const auto & view : sfm_data.GetViews())
  {
    if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
      valid_view_ids.insert(view.first);
  }
  KeepOnlyReferencedElement(valid_view_ids, tripletWise_matches);

  return b_translation;
}

bool GlobalSfM_Translation_AveragingSolver::Translation_averaging(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR)
{
  //-------------------
  //-- GLOBAL TRANSLATIONS ESTIMATION from initial triplets t_ij guess
  //-------------------

  // Keep the largest Biedge connected component graph of relative translations
  Pair_Set pairs;
  std::transform(_vec_initialRijTijEstimates.begin(), _vec_initialRijTijEstimates.end(),
    std::inserter(pairs, pairs.begin()), stl::RetrieveKey());
  const std::set<IndexT> set_remainingIds =
    openMVG::graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs, "./");
  KeepOnlyReferencedElement(set_remainingIds, _vec_initialRijTijEstimates);

  const std::string _sOutDirectory("./");
  {
    const std::set<IndexT> index = getIndexT(_vec_initialRijTijEstimates);

    const size_t iNview = index.size();
    std::cout << "\n-------------------------------" << "\n"
      << " Global translations computation: " << "\n"
      << "   - Ready to compute " << iNview << " global translations." << "\n"
      << "     from #relative translations: " << _vec_initialRijTijEstimates.size() << std::endl;

    if (iNview < 3)
    {
      // Too tiny image set to perform motion averaging
      return false;
    }
    //-- Update initial estimates from [minId,maxId] to range [0->Ncam]
    RelativeInfo_Vec vec_initialRijTijEstimates_cpy = _vec_initialRijTijEstimates;
    const Pair_Set pairs = getPairs(vec_initialRijTijEstimates_cpy);
    Hash_Map<IndexT,IndexT> _reindexForward, _reindexBackward;
    reindex(pairs, _reindexForward, _reindexBackward);
    for(size_t i = 0; i < vec_initialRijTijEstimates_cpy.size(); ++i)
    {
      openMVG::relativeInfo & rel = vec_initialRijTijEstimates_cpy[i];
      rel.first = Pair(_reindexForward[rel.first.first], _reindexForward[rel.first.second]);
    }

    openMVG::system::Timer timerLP_translation;

    switch(eTranslationAveragingMethod)
    {
      case TRANSLATION_AVERAGING_L1:
      {
        double gamma = -1.0;
        std::vector<double> vec_solution;
        {
          vec_solution.resize(iNview*3 + vec_initialRijTijEstimates_cpy.size()/3 + 1);
          using namespace openMVG::linearProgramming;
          #ifdef OPENMVG_HAVE_MOSEK
            MOSEK_SolveWrapper solverLP(vec_solution.size());
          #else
            OSI_CLP_SolverWrapper solverLP(vec_solution.size());
          #endif

          lInfinityCV::Tifromtij_ConstraintBuilder_OneLambdaPerTrif cstBuilder(vec_initialRijTijEstimates_cpy);

          LP_Constraints_Sparse constraint;
          //-- Setup constraint and solver
          cstBuilder.Build(constraint);
          solverLP.setup(constraint);
          //--
          // Solving
          const bool bFeasible = solverLP.solve();
          std::cout << " \n Feasibility " << bFeasible << std::endl;
          //--
          if (bFeasible)  {
            solverLP.getSolution(vec_solution);
            gamma = vec_solution[vec_solution.size()-1];
          }
          else  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
          }
        }

        const double timeLP_translation = timerLP_translation.elapsed();
        //-- Export triplet statistics:
        {

          std::ostringstream os;
          os << "Translation fusion statistics.";
          os.str("");
          os << "-------------------------------" << "\n"
            << "-- #relative estimates: " << vec_initialRijTijEstimates_cpy.size()
            << " converge with gamma: " << gamma << ".\n"
            << " timing (s): " << timeLP_translation << ".\n"
            << "-------------------------------" << "\n";
          std::cout << os.str() << std::endl;
        }

        std::cout << "Found solution:\n";
        std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

        std::vector<double> vec_camTranslation(iNview*3,0);
        std::copy(&vec_solution[0], &vec_solution[iNview*3], &vec_camTranslation[0]);

        std::vector<double> vec_camRelLambdas(&vec_solution[iNview*3], &vec_solution[iNview*3 + vec_initialRijTijEstimates_cpy.size()/3]);
        std::cout << "\ncam position: " << std::endl;
        std::copy(vec_camTranslation.begin(), vec_camTranslation.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\ncam Lambdas: " << std::endl;
        std::copy(vec_camRelLambdas.begin(), vec_camRelLambdas.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;

        // Update the view poses according the found camera centers
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
          const IndexT pose_id = _reindexBackward[i];
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_SOFTL1:
      {
        std::vector<Vec3> vec_translations;
        if (!solve_translations_problem_softl1(
          vec_initialRijTijEstimates_cpy, true, iNview, vec_translations))
        {
          std::cerr << "Compute global translations: failed" << std::endl;
          return false;
        }

        // A valid solution was found:
        // - Update the view poses according the found camera translations
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 & t = vec_translations[i];
          const IndexT pose_id = _reindexBackward[i];
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_L2_DISTANCE_CHORDAL:
      {
        std::vector<int> vec_edges;
        vec_edges.reserve(vec_initialRijTijEstimates_cpy.size() * 2);
        std::vector<double> vec_poses;
        vec_poses.reserve(vec_initialRijTijEstimates_cpy.size() * 3);
        std::vector<double> vec_weights;
        vec_weights.reserve(vec_initialRijTijEstimates_cpy.size());

        for(int i=0; i < vec_initialRijTijEstimates_cpy.size(); ++i)
        {
          const openMVG::relativeInfo & rel = vec_initialRijTijEstimates_cpy[i];
          vec_edges.push_back(rel.first.first);
          vec_edges.push_back(rel.first.second);
          // Since index have been remapped
          // (use the backward indexing to retrieve the second global rotation)
          const IndexT secondId = _reindexBackward[rel.first.second];
          const View * view = sfm_data.views.at(secondId).get();
          const Mat3 & Ri = map_globalR.at(view->id_pose);
          const Vec3 direction = -(Ri.transpose() * rel.second.second.normalized());

          vec_poses.push_back(direction(0));
          vec_poses.push_back(direction(1));
          vec_poses.push_back(direction(2));

          vec_weights.push_back(1.0);
        }

        const double function_tolerance = 1e-7, parameter_tolerance = 1e-8;
        const int max_iterations = 500;

        const double loss_width = 0.0; // No loss in order to compare with TRANSLATION_AVERAGING_L1

        std::vector<double> X(iNview*3, 0.0);
        if(!solve_translations_problem_l2_chordal(
          &vec_edges[0],
          &vec_poses[0],
          &vec_weights[0],
          vec_initialRijTijEstimates_cpy.size(),
          loss_width,
          &X[0],
          function_tolerance,
          parameter_tolerance,
          max_iterations))  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
        }

        // Update camera center for each view
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 C(X[i*3], X[i*3+1], X[i*3+2]);
          const IndexT pose_id = _reindexBackward[i]; // undo the reindexing
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, C);
        }
      }
      break;
      default:
      {
        std::cerr << "Unknown translation averaging method" << std::endl;
        return false;
      }
    }
  }
  return true;
}

void GlobalSfM_Translation_AveragingSolver::Compute_translations(
  const SfM_Data & sfm_data,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches &tripletWise_matches)
{
  std::cout << "\n-------------------------------" << "\n"
    << " Relative translations computation: " << "\n"
    << "-------------------------------" << std::endl;

  // Compute relative translations over the graph of global rotations
  //  thanks to an edge coverage algorithm
  ComputePutativeTranslation_EdgesCoverage(
    sfm_data,
    map_globalR,
    normalized_features_provider,
    matches_provider,
    _vec_initialRijTijEstimates,
    tripletWise_matches);
}

//-- Perform a trifocal estimation of the graph contain in vec_triplets with an
// edge coverage algorithm. Its complexity is sub-linear in term of edges count.
void GlobalSfM_Translation_AveragingSolver::ComputePutativeTranslation_EdgesCoverage(
  const SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  RelativeInfo_Vec & vec_initialEstimates,
  matching::PairWiseMatches & newpairMatches)
{
  openMVG::system::Timer timerLP_triplet;

  //--
  // Compute the relative translations using triplets of rotations over the rotation graph.
  //--
  //
  // 1. List plausible triplets over the global rotation pose graph Ids.
  //   - list all edges that have support in the rotation pose graph
  //
  Pair_Set rotation_pose_id_graph;
  std::set<IndexT> set_pose_ids;
  std::transform(map_globalR.begin(), map_globalR.end(),
    std::inserter(set_pose_ids, set_pose_ids.begin()), stl::RetrieveKey());
  // List shared correspondences (pairs) between poses
  for (const auto & match_iterator : matches_provider->_pairWise_matches)
  {
    const Pair pair = match_iterator.first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();
    // Consider iff the pair is supported by the rotation graph
    if (v1->id_pose != v2->id_pose &&
        set_pose_ids.count(v1->id_pose) && set_pose_ids.count(v2->id_pose))
    {
      rotation_pose_id_graph.insert(
        std::make_pair(v1->id_pose, v2->id_pose));
    }
  }
  // List putative triplets (from global rotations Ids)
  const std::vector< graph::Triplet > vec_triplets =
    graph::tripletListing(rotation_pose_id_graph);
  std::cout << "#Triplets: " << vec_triplets.size() << std::endl;

  {
    // Compute triplets of translations
    // Avoid to cover each edge of the graph by using an edge coverage algorithm
    // An estimated triplets of translation mark three edges as estimated.

    //-- precompute the number of track per triplet:
    Hash_Map<IndexT, IndexT> map_tracksPerTriplets;
    #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < (int)vec_triplets.size(); ++i)
    {
      // List matches that belong to the triplet of poses
      const graph::Triplet & triplet = vec_triplets[i];
      PairWiseMatches map_triplet_matches;
      const std::set<IndexT> set_pose_ids = {triplet.i, triplet.j, triplet.k};
      // List shared correspondences (pairs) between the triplet poses
      for (const auto & match_iterator : matches_provider->_pairWise_matches)
      {
        const Pair pair = match_iterator.first;
        const View * v1 = sfm_data.GetViews().at(pair.first).get();
        const View * v2 = sfm_data.GetViews().at(pair.second).get();
        // Consider iff the pair is supported by the triplet & rotation graph
        const bool b_different_pose_id = v1->id_pose != v2->id_pose;
        const int covered_pose =
          set_pose_ids.count(v1->id_pose) +
          set_pose_ids.count(v2->id_pose);
        // Different pose Id and the current edge cover the triplet edge
        if (b_different_pose_id && covered_pose == 2 )
        {
          map_triplet_matches.insert( match_iterator );
        }
      }
      // Compute tracks:
      {
        openMVG::tracks::TracksBuilder tracksBuilder;
        tracksBuilder.Build(map_triplet_matches);
        tracksBuilder.Filter(3);
        #ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
        #endif
        map_tracksPerTriplets[i] = tracksBuilder.NbTracks(); //count the # of matches in the UF tree
      }
    }

    typedef Pair myEdge;

    //-- List all edges
    std::set<myEdge > set_edges;

    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graph::Triplet & triplet = vec_triplets[i];
      const IndexT I = triplet.i, J = triplet.j , K = triplet.k;
      // Add three edges
      set_edges.insert(std::make_pair(I,J));
      set_edges.insert(std::make_pair(I,K));
      set_edges.insert(std::make_pair(J,K));
    }
    // Move set to a vector
    std::vector<myEdge > vec_edges(std::begin(set_edges), std::end(set_edges));
    std::set<myEdge >().swap(set_edges); // release memory

    openMVG::sfm::MutexSet<myEdge> m_mutexSet;

    C_Progress_display my_progress_bar(
      vec_edges.size(),
      std::cout,
      "\nRelative translations computation (edge coverage algorithm)\n");

    bool bVerbose = false;

    #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
    #endif
    for (int k = 0; k < vec_edges.size(); ++k)
    {
      #ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
      #endif
      {
        ++my_progress_bar;
      }

      const myEdge & edge = vec_edges[k];
      //-- If current edge already computed continue
      if (m_mutexSet.count(edge) || m_mutexSet.size() == vec_edges.size())
      {
        if (bVerbose)
          std::cout << "EDGES WAS PREVIOUSLY COMPUTED" << std::endl;
        continue;
      }

      std::vector<size_t> vec_possibleTriplets;
      // Find the triplets that contain the given edge
      for (size_t i = 0; i < vec_triplets.size(); ++i)
      {
        const graph::Triplet & triplet = vec_triplets[i];
        if (triplet.contain(edge))
        {
          vec_possibleTriplets.push_back(i);
        }
      }

      //-- Sort the triplet according the number of matches they have on their edges
      std::vector<size_t> vec_commonTracksPerTriplets;
      for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
      {
        vec_commonTracksPerTriplets.push_back(map_tracksPerTriplets[vec_possibleTriplets[i]]);
      }
      //-- If current edge is already computed continue
      if (m_mutexSet.count(edge))
        continue;

      using namespace stl::indexed_sort;
      std::vector< sort_index_packet_descend < size_t, size_t> > packet_vec(vec_commonTracksPerTriplets.size());
      sort_index_helper(packet_vec, &vec_commonTracksPerTriplets[0]);

      std::vector<size_t> vec_possibleTripletsSorted;
      for (size_t i = 0; i < vec_commonTracksPerTriplets.size(); ++i) {
        vec_possibleTripletsSorted.push_back( vec_possibleTriplets[packet_vec[i].index] );
      }
      vec_possibleTriplets.swap(vec_possibleTripletsSorted);

      // Try to solve the triplets
      // Search the possible triplet:
      for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
      {
        const graph::Triplet & triplet = vec_triplets[vec_possibleTriplets[i]];

        // If the triplet is already estimated by another thread; try a new one
        if (m_mutexSet.count(Pair(triplet.i, triplet.j)) &&
            m_mutexSet.count(Pair(triplet.i, triplet.k)) &&
            m_mutexSet.count(Pair(triplet.j, triplet.k)))
        {
          continue;
        }

        //--
        // Try to estimate this triplet.
        //--
        double dPrecision = 4.0; // upper bound of the pixel residual

        std::vector<Vec3> vec_tis(3);
        std::vector<size_t> vec_inliers;
        openMVG::tracks::STLMAPTracks pose_triplet_tracks;

        const std::string sOutDirectory = "./";
        const bool bTriplet_estimation = Estimate_T_triplet(
            sfm_data,
            map_globalR,
            normalized_features_provider,
            matches_provider,
            triplet,
            vec_tis,
            dPrecision,
            vec_inliers,
            pose_triplet_tracks,
            sOutDirectory);

        if (bTriplet_estimation)
        {

          if (m_mutexSet.count(Pair(triplet.i, triplet.j)) &&
              m_mutexSet.count(Pair(triplet.i, triplet.k)) &&
              m_mutexSet.count(Pair(triplet.j, triplet.k)))
          {
            // If triplet already estimated by another thread;
            //  stop, since the edge have already a translation estimate.
            break;
          }

          // Since the triplet is estimated, mark the edges as estimated
          m_mutexSet.insert(std::make_pair(triplet.i, triplet.j));
          m_mutexSet.insert(std::make_pair(triplet.j, triplet.k));
          m_mutexSet.insert(std::make_pair(triplet.i, triplet.k));

          // Compute the triplet relative motions (IJ, JK, IK)
          const Mat3 RI = map_globalR.at(triplet.i);
          const Mat3 RJ = map_globalR.at(triplet.j);
          const Mat3 RK = map_globalR.at(triplet.k);
          const Vec3 ti = vec_tis[0];
          const Vec3 tj = vec_tis[1];
          const Vec3 tk = vec_tis[2];

          //--- ATOMIC
          #ifdef OPENMVG_USE_OPENMP
             #pragma omp critical
          #endif
          {
            Mat3 RijGt;
            Vec3 tij;
            RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
            vec_initialEstimates.emplace_back(
              std::make_pair(triplet.i, triplet.j), std::make_pair(RijGt, tij));

            Mat3 RjkGt;
            Vec3 tjk;
            RelativeCameraMotion(RJ, tj, RK, tk, &RjkGt, &tjk);
            vec_initialEstimates.emplace_back(
              std::make_pair(triplet.j, triplet.k), std::make_pair(RjkGt, tjk));

            Mat3 RikGt;
            Vec3 tik;
            RelativeCameraMotion(RI, ti, RK, tk, &RikGt, &tik);
            vec_initialEstimates.emplace_back(
              std::make_pair(triplet.i, triplet.k), std::make_pair(RikGt, tik));

            // Add inliers as valid pairwise matches
            for (std::vector<size_t>::const_iterator iterInliers = vec_inliers.begin();
              iterInliers != vec_inliers.end(); ++iterInliers)
            {
              using namespace openMVG::tracks;
              const submapTrack & subTrack = pose_triplet_tracks.at(*iterInliers);

              // create pairwise matches from inlier track
              for (size_t index_I = 0; index_I < subTrack.size() ; ++index_I)
              { submapTrack::const_iterator iter_I = subTrack.begin();
                std::advance(iter_I, index_I);

                // extract camera indexes
                const size_t id_view_I = iter_I->first;
                const size_t id_feat_I = iter_I->second;

                // loop on subtracks
                for (size_t index_J = index_I+1; index_J < subTrack.size() ; ++index_J)
                { submapTrack::const_iterator iter_J = subTrack.begin();
                  std::advance(iter_J, index_J);

                  // extract camera indexes
                  const size_t id_view_J = iter_J->first;
                  const size_t id_feat_J = iter_J->second;

                  newpairMatches[std::make_pair(id_view_I, id_view_J)].emplace_back(id_feat_I, id_feat_J);
                }
              }
            }
          }
          // Since a relative translation have been found for the edge: vec_edges[k],
          //  we break and start to estimate the translations for some other edges.
          break;
        }
      }
    }
  }

  const double timeLP_triplet = timerLP_triplet.elapsed();
  std::cout << "TRIPLET COVERAGE TIMING: " << timeLP_triplet << " seconds" << std::endl;

  std::cout << "-------------------------------" << "\n"
      << "-- #Relative translations estimates: " << _vec_initialRijTijEstimates.size()/3
      << " computed from " << vec_triplets.size() << " triplets.\n"
      << "-- resulting in " << _vec_initialRijTijEstimates.size() << " translations estimation.\n"
      << "-- timing to obtain the relative translations: " << timeLP_triplet << " seconds.\n"
      << "-------------------------------" << std::endl;
}

// Robust estimation and refinement of a translation and 3D points of an image triplets.
bool GlobalSfM_Translation_AveragingSolver::Estimate_T_triplet(
  const SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const Features_Provider * normalized_features_provider,
  const Matches_Provider * matches_provider,
  const graph::Triplet & poses_id,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<size_t> & vec_inliers,
  openMVG::tracks::STLMAPTracks & tracks,
  const std::string & sOutDirectory) const
{
  // List matches that belong to the triplet of poses
  PairWiseMatches map_triplet_matches;
  std::set<IndexT> set_pose_ids;
  set_pose_ids.insert(poses_id.i);
  set_pose_ids.insert(poses_id.j);
  set_pose_ids.insert(poses_id.k);
  // List shared correspondences (pairs) between poses
  for (const auto & match_iterator : matches_provider->_pairWise_matches)
  {
    const Pair pair = match_iterator.first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();
    // Consider the pair iff it is supported by the triplet & rotation graph
    const bool b_different_pose_id = v1->id_pose != v2->id_pose;
    const int covered_pose =
      set_pose_ids.count(v1->id_pose) +
      set_pose_ids.count(v2->id_pose);
    // Different pose Id and the current edge cover the triplet edge
    if ( b_different_pose_id && covered_pose == 2 )
    {
      map_triplet_matches.insert( match_iterator );
    }
  }

  openMVG::tracks::TracksBuilder tracksBuilder;
  tracksBuilder.Build(map_triplet_matches);
  tracksBuilder.Filter(3);
  tracksBuilder.ExportToSTL(tracks);

  if (tracks.size() < 50)
    return false;

  // Convert data
  Mat x1(2, tracks.size());
  Mat x2(2, tracks.size());
  Mat x3(2, tracks.size());

  Mat* xxx[3] = {&x1, &x2, &x3};
  std::set<IndexT> intrinsic_ids;
  size_t cpt = 0;
  for (tracks::STLMAPTracks::const_iterator iterTracks = tracks.begin();
    iterTracks != tracks.end(); ++iterTracks, ++cpt) {
    const tracks::submapTrack & subTrack = iterTracks->second;
    size_t index = 0;
    for (tracks::submapTrack::const_iterator iter = subTrack.begin(); iter != subTrack.end(); ++iter, ++index) {
      const size_t idx_view = iter->first;
      const features::PointFeature pt = normalized_features_provider->getFeatures(idx_view)[iter->second];
      xxx[index]->col(cpt) = pt.coords().cast<double>();
      const View * view = sfm_data.views.at(idx_view).get();
      intrinsic_ids.insert(view->id_intrinsic);
    }
  }
  // Retrieve the smallest focal value, for threshold normalization
  double min_focal = std::numeric_limits<double>::max();
  for (const auto & ids : intrinsic_ids)
  {
    const cameras::IntrinsicBase * intrinsicPtr = sfm_data.GetIntrinsics().at(ids).get();
    const cameras::Pinhole_Intrinsic * intrinsic = dynamic_cast< const cameras::Pinhole_Intrinsic * > (intrinsicPtr);
    if (intrinsic)
    {
      min_focal = std::min(min_focal, intrinsic->focal());
    }
  }
  if (min_focal == std::numeric_limits<double>::max())
  {
    return false;
  }

  // Get rotations:
  const std::vector<Mat3> vec_global_R_Triplet =
    {map_globalR.at(poses_id.i), map_globalR.at(poses_id.j), map_globalR.at(poses_id.k)};

  using namespace openMVG::trifocal;
  using namespace openMVG::trifocal::kernel;

  typedef TranslationTripletKernel_ACRansac<
    translations_Triplet_Solver,
    translations_Triplet_Solver,
    TrifocalTensorModel> KernelType;
  const double ThresholdUpperBound = 1.0e-2; // upper bound of the pixel residual (normalized coordinates)
  KernelType kernel(x1, x2, x3, vec_global_R_Triplet, Mat3::Identity(), ThresholdUpperBound);

  const size_t ORSA_ITER = 320;  // max number of iterations of AC-RANSAC

  TrifocalTensorModel T;
  std::pair<double,double> acStat = robust::ACRANSAC(kernel, vec_inliers, ORSA_ITER, &T, dPrecision/min_focal, false);
  // If robust estimation fail => stop.
  if (dPrecision == std::numeric_limits<double>::infinity())
    return false;

  // Update output parameters
  dPrecision = acStat.first * min_focal;

  vec_tis.resize(3);
  Mat3 K, R;
  KRt_From_P(T.P1, &K, &R, &vec_tis[0]);
  KRt_From_P(T.P2, &K, &R, &vec_tis[1]);
  KRt_From_P(T.P3, &K, &R, &vec_tis[2]);

#ifdef DEBUG_TRIPLET
  // compute 3D scene base on motion estimation
  SfM_Data    tiny_scene;

  // intialize poses (which are now shared by a group of images)
  tiny_scene.poses[poses_id.i] = Pose3(vec_global_R_Triplet[0], -vec_global_R_Triplet[0].transpose() * vec_tis[0]);
  tiny_scene.poses[poses_id.j] = Pose3(vec_global_R_Triplet[1], -vec_global_R_Triplet[1].transpose() * vec_tis[1]);
  tiny_scene.poses[poses_id.k] = Pose3(vec_global_R_Triplet[2], -vec_global_R_Triplet[2].transpose() * vec_tis[2]);

  // insert views used by the relative pose pairs
  for (const auto & pairIterator : map_triplet_matches )
  {
    // initialize camera indexes
    const IndexT I = pairIterator.first.first;
    const IndexT J = pairIterator.first.second;

    // add views
    tiny_scene.views.insert(*sfm_data.GetViews().find(I));
    tiny_scene.views.insert(*sfm_data.GetViews().find(J));

    // add intrinsics
    const View * view_I = sfm_data.GetViews().at(I).get();
    const View * view_J = sfm_data.GetViews().at(J).get();
    tiny_scene.intrinsics.insert(*sfm_data.GetIntrinsics().find(view_I->id_intrinsic));
    tiny_scene.intrinsics.insert(*sfm_data.GetIntrinsics().find(view_J->id_intrinsic));
  }

  // Fill sfm_data with the inliers tracks. Feed image observations: no 3D yet.
  Landmarks & structure = tiny_scene.structure;
  for (size_t idx=0; idx < vec_inliers.size(); ++idx)
  {
    const size_t trackId = vec_inliers[idx];
    const tracks::submapTrack & track = tracks.at(trackId);
    Observations & obs = structure[idx].obs;
    for (tracks::submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
    {
      // get view Id and feat ID
      const size_t viewIndex = it->first;
      const size_t featIndex = it->second;

      // initialize view and get intrinsics
      const View * view = sfm_data.GetViews().at(viewIndex).get();
      const cameras::IntrinsicBase *  cam = sfm_data.GetIntrinsics().find(view->id_intrinsic)->second.get();
      const cameras::Pinhole_Intrinsic * intrinsicPtr = dynamic_cast< const cameras::Pinhole_Intrinsic * >(cam);
      const Vec2 principal_point = intrinsicPtr->principal_point();

      // get normalized feature
      const features::PointFeature & pt = normalized_features_provider->feats_per_view.at(viewIndex)[featIndex];
      const Vec2 pt_unnormalized (cam->cam2ima(pt.coords().cast<double>()));
      obs[viewIndex] = Observation(pt_unnormalized, featIndex);
    }
  }

  // Compute 3D landmark positions (triangulation of the observations)
  {
    SfM_Data_Structure_Computation_Blind structure_estimator(false);
    structure_estimator.triangulate(tiny_scene);
  }

  // Refine structure and poses (keep intrinsic constant)
  Bundle_Adjustment_Ceres::BA_options options(false, false);
  options._linear_solver_type = ceres::SPARSE_SCHUR;
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  if (bundle_adjustment_obj.Adjust(tiny_scene, false, true, false))
  {
    // export scene for visualization
    std::ostringstream os;
    os << poses_id.i << "_" << poses_id.j << "_" << poses_id.k << ".ply";
    Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));

    // Export refined translations
    vec_tis[0] = tiny_scene.poses[poses_id.i].translation();
    vec_tis[1] = tiny_scene.poses[poses_id.j].translation();
    vec_tis[2] = tiny_scene.poses[poses_id.k].translation();
  }
#endif

  // Keep model only if it has a sufficient inlier coverage
  const bool bTest = ( vec_inliers.size() > 30 && 0.33 * tracks.size() );

#ifdef DEBUG_TRIPLET
  {
    std::cout << "Triplet : status: " << bTest
      << " AC: " << dPrecision
      << " inliers % " << double(vec_inliers.size()) / tracks.size() * 100.0
      << " total putative " << tracks.size() << std::endl;
  }
#endif

  return bTest;
}

} // namespace sfm
} // namespace openMVG

