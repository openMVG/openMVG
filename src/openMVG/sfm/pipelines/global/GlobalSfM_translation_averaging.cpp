// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include "openMVG/graph/graph.hpp"
#include "openMVG/types.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/sfm/pipelines/global/mutexSet.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/sfm/pipelines/global/triplet_t_ACRansac_kernelAdaptator.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include <vector>

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::matching;

/// Use features in normalized camera frames
bool GlobalSfM_Translation_AveragingSolver::Run
(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  openMVG::sfm::SfM_Data & sfm_data,
  const openMVG::sfm::Features_Provider * features_provider,
  const openMVG::sfm::Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Compute the relative translations and save them to vec_initialRijTijEstimates:
  Compute_translations(
    sfm_data,
    features_provider,
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
  sfm::SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR)
{
  //-------------------
  //-- GLOBAL TRANSLATIONS ESTIMATION from initial triplets t_ij guess
  //-------------------

  // Keep the largest Biedge connected component graph of relative translations
  Pair_Set pairs;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_motion_)
  {
    for (const relativeInfo & rel : iter)
    {
      pairs.insert(rel.first);
    }
  }
  const std::set<IndexT> set_remainingIds =
    openMVG::graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
  KeepOnlyReferencedElement(set_remainingIds, vec_relative_motion_);

  {
    const std::set<IndexT> index = getIndexT(vec_relative_motion_);

    const size_t iNview = index.size();
    std::cout << "\n-------------------------------" << "\n"
      << " Global translations computation: " << "\n"
      << "   - Ready to compute " << iNview << " global translations." << "\n"
      << "     from #relative translations: " << vec_relative_motion_.size()*3 << std::endl;

    if (iNview < 3)
    {
      // Too tiny image set to perform motion averaging
      return false;
    }
    //-- Update initial estimates from [minId,maxId] to range [0->Ncam]
    std::vector<RelativeInfo_Vec> vec_relative_motion_cpy = vec_relative_motion_;
    const Pair_Set pairs = getPairs(vec_relative_motion_cpy);
    Hash_Map<IndexT,IndexT> reindex_forward, reindex_backward;
    reindex(pairs, reindex_forward, reindex_backward);
    for (openMVG::RelativeInfo_Vec & iter : vec_relative_motion_cpy)
    {
      for (relativeInfo & it : iter)
      {
        it.first = Pair(reindex_forward[it.first.first], reindex_forward[it.first.second]);
      }
    }

    openMVG::system::Timer timerLP_translation;

    switch (eTranslationAveragingMethod)
    {
      case TRANSLATION_AVERAGING_L1:
      {
        double gamma = -1.0;
        std::vector<double> vec_solution;
        {
          vec_solution.resize(iNview*3 + vec_relative_motion_cpy.size() + 1);
          using namespace openMVG::linearProgramming;
          OSI_CLP_SolverWrapper solverLP(vec_solution.size());

          lInfinityCV::Tifromtij_ConstraintBuilder cstBuilder(vec_relative_motion_cpy);

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
            << "-- #relative estimates: " << vec_relative_motion_cpy.size()
            << " converge with gamma: " << gamma << ".\n"
            << " timing (s): " << timeLP_translation << ".\n"
            << "-------------------------------" << "\n";
          std::cout << os.str() << std::endl;
        }

        std::cout << "Found solution:\n";
        std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

        std::vector<double> vec_camTranslation(iNview*3,0);
        std::copy(&vec_solution[0], &vec_solution[iNview*3], &vec_camTranslation[0]);

        std::vector<double> vec_camRelLambdas(&vec_solution[iNview*3], &vec_solution[iNview*3 + vec_relative_motion_cpy.size()]);
        std::cout << "\ncam position: " << std::endl;
        std::copy(vec_camTranslation.begin(), vec_camTranslation.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\ncam Lambdas: " << std::endl;
        std::copy(vec_camRelLambdas.begin(), vec_camRelLambdas.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;

        // Update the view poses according the found camera centers
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 t(vec_camTranslation[i*3], vec_camTranslation[i*3+1], vec_camTranslation[i*3+2]);
          const IndexT pose_id = reindex_backward[i];
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_SOFTL1:
      {
        std::vector<Vec3> vec_translations;
        if (!solve_translations_problem_softl1(
          vec_relative_motion_cpy, vec_translations))
        {
          std::cerr << "Compute global translations: failed" << std::endl;
          return false;
        }

        // A valid solution was found:
        // - Update the view poses according the found camera translations
        for (size_t i = 0; i < iNview; ++i)
        {
          const Vec3 & t = vec_translations[i];
          const IndexT pose_id = reindex_backward[i];
          const Mat3 & Ri = map_globalR.at(pose_id);
          sfm_data.poses[pose_id] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_L2_DISTANCE_CHORDAL:
      {
        std::vector<int> vec_edges;
        vec_edges.reserve(vec_relative_motion_cpy.size() * 2);
        std::vector<double> vec_poses;
        vec_poses.reserve(vec_relative_motion_cpy.size() * 3);
        std::vector<double> vec_weights;
        vec_weights.reserve(vec_relative_motion_cpy.size());

        for (const openMVG::RelativeInfo_Vec & iter : vec_relative_motion_cpy)
        {
          for (const relativeInfo & rel : iter)
          {
            vec_edges.push_back(rel.first.first);
            vec_edges.push_back(rel.first.second);
            // Since index have been remapped
            // (use the backward indexing to retrieve the second global rotation)
            const IndexT secondId = reindex_backward[rel.first.second];
            const View * view = sfm_data.views.at(secondId).get();
            const Mat3 & Ri = map_globalR.at(view->id_pose);
            const Vec3 direction = -(Ri.transpose() * rel.second.second.normalized());

            vec_poses.push_back(direction(0));
            vec_poses.push_back(direction(1));
            vec_poses.push_back(direction(2));

            vec_weights.push_back(1.0);
          }
        }

        const double function_tolerance = 1e-7, parameter_tolerance = 1e-8;
        const int max_iterations = 500;

        const double loss_width = 0.0; // No loss in order to compare with TRANSLATION_AVERAGING_L1

        std::vector<double> X(iNview*3, 0.0);
        if (!solve_translations_problem_l2_chordal(
          &vec_edges[0],
          &vec_poses[0],
          &vec_weights[0],
          vec_relative_motion_cpy.size()*3,
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
          const IndexT pose_id = reindex_backward[i]; // undo the reindexing
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

void GlobalSfM_Translation_AveragingSolver::Compute_translations
(
  const sfm::SfM_Data & sfm_data,
  const sfm::Features_Provider * features_provider,
  const sfm::Matches_Provider * matches_provider,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches &tripletWise_matches
)
{
  std::cout << "\n-------------------------------" << "\n"
    << " Relative translations computation: " << "\n"
    << "-------------------------------" << std::endl;

  // Compute relative translations over the graph of global rotations
  //  thanks to an edge coverage algorithm
  ComputePutativeTranslation_EdgesCoverage(
    sfm_data,
    map_globalR,
    features_provider,
    matches_provider,
    vec_relative_motion_,
    tripletWise_matches);
}

//-- Perform a trifocal estimation of the graph contain in vec_triplets with an
// edge coverage algorithm. Its complexity is sub-linear in term of edges count.
void GlobalSfM_Translation_AveragingSolver::ComputePutativeTranslation_EdgesCoverage
(
  const sfm::SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const sfm::Features_Provider * features_provider,
  const sfm::Matches_Provider * matches_provider,
  std::vector<RelativeInfo_Vec> & vec_triplet_relative_motion,
  matching::PairWiseMatches & newpairMatches
)
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
  std::transform(map_globalR.cbegin(), map_globalR.cend(),
    std::inserter(set_pose_ids, set_pose_ids.begin()), stl::RetrieveKey());
  // List shared correspondences (pairs) between poses
  for (const auto & match_iterator : matches_provider->pairWise_matches_)
  {
    const Pair pair = match_iterator.first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();

    if (// Consider the pair iff it is supported by the rotation graph
        (v1->id_pose != v2->id_pose)
        && set_pose_ids.count(v1->id_pose)
        && set_pose_ids.count(v2->id_pose))
    {
      rotation_pose_id_graph.insert({v1->id_pose, v2->id_pose});
    }
  }
  // List putative triplets (from global rotations Ids)
  const std::vector<graph::Triplet> vec_triplets =
    graph::TripletListing(rotation_pose_id_graph);
  std::cout << "#Triplets: " << vec_triplets.size() << std::endl;

  {
    // Compute triplets of translations
    // Avoid to cover each edge of the graph by using an edge coverage algorithm
    // An estimated triplets of translation mark three edges as estimated.

    //-- Alias (list triplet ids used per edges)
    using myEdge = Pair; // An edge between two pose id
    Hash_Map<myEdge, std::vector<uint32_t>> map_tripletIds_perEdge;
    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graph::Triplet & triplet = vec_triplets[i];
      map_tripletIds_perEdge[{triplet.i, triplet.j}].push_back(i);
      map_tripletIds_perEdge[{triplet.i, triplet.k}].push_back(i);
      map_tripletIds_perEdge[{triplet.j, triplet.k}].push_back(i);
    }

    //-- precompute the visibility count per triplets (sum of their 2 view matches)
    Hash_Map<IndexT, IndexT> map_tracksPerTriplets;
    for (const auto & match_iterator : matches_provider->pairWise_matches_)
    {
      const Pair pair = match_iterator.first;
      const View * v1 = sfm_data.GetViews().at(pair.first).get();
      const View * v2 = sfm_data.GetViews().at(pair.second).get();
      if (v1->id_pose != v2->id_pose)
      {
        // Consider the pair iff it is supported by 2 different pose id
        const myEdge edge(v1->id_pose,v2->id_pose);
        if (map_tripletIds_perEdge.count(edge) != 0)
        {
          const auto & edge_tripletIds = map_tripletIds_perEdge.at(edge);
          for (const auto & triplet_id : edge_tripletIds)
          {
            if (map_tracksPerTriplets.count(triplet_id) == 0)
              map_tracksPerTriplets[triplet_id] = match_iterator.second.size();
            else
              map_tracksPerTriplets[triplet_id] += match_iterator.second.size();
          }
        }
      }
    }
    // Collect edges that are covered by the triplets
    std::vector<myEdge> vec_edges;
    std::transform(map_tripletIds_perEdge.cbegin(),
                   map_tripletIds_perEdge.cend(),
                   std::back_inserter(vec_edges),
                   stl::RetrieveKey());

    openMVG::sfm::MutexSet<myEdge> m_mutexSet;

    C_Progress_display my_progress_bar(
      vec_edges.size(),
      std::cout,
      "\nRelative translations computation (edge coverage algorithm)\n");

#  ifdef OPENMVG_USE_OPENMP
    std::vector<std::vector<RelativeInfo_Vec>> initial_estimates(omp_get_max_threads());
#  else
    std::vector<std::vector<RelativeInfo_Vec>> initial_estimates(1);
#  endif

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int k = 0; k < static_cast<int>(vec_edges.size()); ++k)
    {
      const myEdge & edge = vec_edges[k];
      ++my_progress_bar;

      if (m_mutexSet.count(edge) == 0 && m_mutexSet.size() != vec_edges.size())
      {
        // Find the triplets that are supporting the given edge
        const auto & vec_possibleTripletIndexes = map_tripletIds_perEdge.at(edge);

        //-- Sort the triplets according their number of track
        std::vector<uint32_t> vec_commonTracksPerTriplets;
        for (const uint32_t triplet_index : vec_possibleTripletIndexes)
        {
          vec_commonTracksPerTriplets.push_back(map_tracksPerTriplets[triplet_index]);
        }

        using namespace stl::indexed_sort;
        std::vector<sort_index_packet_descend<uint32_t,uint32_t>> packet_vec(vec_commonTracksPerTriplets.size());
        sort_index_helper(packet_vec, &vec_commonTracksPerTriplets[0]);

        std::vector<uint32_t> vec_triplet_ordered(vec_commonTracksPerTriplets.size(), 0);
        for (size_t i = 0; i < vec_commonTracksPerTriplets.size(); ++i) {
          vec_triplet_ordered[i] = vec_possibleTripletIndexes[packet_vec[i].index];
        }

        // Try to solve a triplet of translations for the given edge
        for (const uint32_t triplet_index : vec_triplet_ordered)
        {
          const graph::Triplet & triplet = vec_triplets[triplet_index];

          // If the triplet is already estimated by another thread; try the next one
          if (m_mutexSet.count(Pair(triplet.i, triplet.j)) &&
              m_mutexSet.count(Pair(triplet.i, triplet.k)) &&
              m_mutexSet.count(Pair(triplet.j, triplet.k)))
          {
            continue;
          }

          //--
          // Try to estimate this triplet of translations
          //--
          double dPrecision = 4.0; // upper bound of the residual pixel reprojection error

          std::vector<Vec3> vec_tis(3);
          std::vector<uint32_t> vec_inliers;
          openMVG::tracks::STLMAPTracks pose_triplet_tracks;

          const std::string sOutDirectory = "./";

          const bool bTriplet_estimation = Estimate_T_triplet(
              sfm_data,
              map_globalR,
              features_provider,
              matches_provider,
              triplet,
              vec_tis,
              dPrecision,
              vec_inliers,
              pose_triplet_tracks,
              sOutDirectory);

          if (bTriplet_estimation)
          {
            // Since new translation edges have been computed, mark their corresponding edges as estimated
            #ifdef OPENMVG_USE_OPENMP
              #pragma omp critical
            #endif
            {
              m_mutexSet.insert({triplet.i, triplet.j});
              m_mutexSet.insert({triplet.j, triplet.k});
              m_mutexSet.insert({triplet.i, triplet.k});
            }

            // Compute the triplet relative motions (IJ, JK, IK)
            {
              const Mat3
                RI = map_globalR.at(triplet.i),
                RJ = map_globalR.at(triplet.j),
                RK = map_globalR.at(triplet.k);
              const Vec3
                ti = vec_tis[0],
                tj = vec_tis[1],
                tk = vec_tis[2];

              Mat3 Rij;
              Vec3 tij;
              RelativeCameraMotion(RI, ti, RJ, tj, &Rij, &tij);

              Mat3 Rjk;
              Vec3 tjk;
              RelativeCameraMotion(RJ, tj, RK, tk, &Rjk, &tjk);

              Mat3 Rik;
              Vec3 tik;
              RelativeCameraMotion(RI, ti, RK, tk, &Rik, &tik);

              #ifdef OPENMVG_USE_OPENMP
                const int thread_id = omp_get_thread_num();
              #else
                const int thread_id = 0;
              #endif

              RelativeInfo_Vec triplet_relative_motion;
              triplet_relative_motion.push_back(
                {{triplet.i, triplet.j}, {Rij, tij}});
              triplet_relative_motion.push_back(
                {{triplet.j, triplet.k}, {Rjk, tjk}});
              triplet_relative_motion.push_back(
                {{triplet.i, triplet.k}, {Rik, tik}});

              initial_estimates[thread_id].emplace_back(triplet_relative_motion);

              #ifdef OPENMVG_USE_OPENMP
                #pragma omp critical
              #endif
              {
                // Add inliers as valid pairwise matches
                for (const uint32_t & inlier_it : vec_inliers)
                {
                  tracks::STLMAPTracks::const_iterator it_tracks = pose_triplet_tracks.begin();
                  std::advance(it_tracks, inlier_it);
                  const tracks::submapTrack & track = it_tracks->second;

                  // create pairwise matches from the inlier track
                  tracks::submapTrack::const_iterator iter_I = track.begin();
                  tracks::submapTrack::const_iterator iter_J = track.begin();
                  std::advance(iter_J, 1);
                  while (iter_J != track.end())
                  { // matches(pair(view_id(I), view_id(J))) <= IndMatch(feat_id(I), feat_id(J))
                    newpairMatches[{iter_I->first, iter_J->first}]
                     .emplace_back(iter_I->second, iter_J->second);
                    ++iter_I;
                    ++iter_J;
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
    // Merge thread(s) estimates
    for (auto & vec : initial_estimates)
    {
      vec_triplet_relative_motion.insert( vec_triplet_relative_motion.end(),
        std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));
    }
  }

  const double timeLP_triplet = timerLP_triplet.elapsed();
  std::cout << "TRIPLET COVERAGE TIMING:\n"
    << "-------------------------------" << "\n"
    << "-- #Relative triplet of translations estimates: " << vec_triplet_relative_motion.size()
    << " computed from " << vec_triplets.size() << " triplets.\n"
    << "-- resulting in " << vec_triplet_relative_motion.size()*3 << " translations estimation.\n"
    << "-- time to compute triplets of relative translations: " << timeLP_triplet << " seconds.\n"
    << "-------------------------------" << std::endl;
}

// Robust estimation and refinement of a triplet of translations
bool GlobalSfM_Translation_AveragingSolver::Estimate_T_triplet
(
  const sfm::SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const sfm::Features_Provider * features_provider,
  const sfm::Matches_Provider * matches_provider,
  const graph::Triplet & poses_id,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<uint32_t> & vec_inliers,
  openMVG::tracks::STLMAPTracks & tracks,
  const std::string & sOutDirectory
) const
{
  // List matches that belong to the triplet of poses
  PairWiseMatches map_triplet_matches;
  const std::set<IndexT> set_pose_ids {poses_id.i, poses_id.j, poses_id.k};
  // List shared correspondences (pairs) between poses
  for (const auto & match_iterator : matches_provider->pairWise_matches_)
  {
    const Pair & pair = match_iterator.first;
    const View * v1 = sfm_data.GetViews().at(pair.first).get();
    const View * v2 = sfm_data.GetViews().at(pair.second).get();
    if (// Consider the pair iff it is supported by the triplet graph & 2 different pose id
        (v1->id_pose != v2->id_pose)
        && set_pose_ids.count(v1->id_pose)
        && set_pose_ids.count(v2->id_pose))
    {
      map_triplet_matches.insert( match_iterator );
    }
  }

  openMVG::tracks::TracksBuilder tracksBuilder;
  tracksBuilder.Build(map_triplet_matches);
  tracksBuilder.Filter(3);
  tracksBuilder.ExportToSTL(tracks);

  if (tracks.size() < 30)
    return false;

  // Data conversion
  // Fill image observations as a unique matrix array
  Mat
    x1(2, tracks.size()),
    x2(2, tracks.size()),
    x3(2, tracks.size());

  Mat* xxx[3] = {&x1, &x2, &x3};
  std::set<IndexT> intrinsic_ids;
  size_t cpt = 0;
  for (const auto & tracks_it : tracks)
  {
    const tracks::submapTrack & track = tracks_it.second;
    uint32_t index = 0;
    for (const auto & track_it : track)
    {
      const uint32_t idx_view = track_it.first;
      const View * view = sfm_data.views.at(idx_view).get();
      const IntrinsicBase * cam = sfm_data.intrinsics.at(view->id_intrinsic).get();
      intrinsic_ids.insert(view->id_intrinsic);
      const features::PointFeature pt = features_provider->getFeatures(idx_view)[track_it.second];
      xxx[index++]->col(cpt) = ((*cam)(cam->get_ud_pixel(pt.coords().cast<double>()))).colwise().hnormalized();
    }
    ++cpt;
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

  using KernelType =
    TranslationTripletKernel_ACRansac<
      translations_Triplet_Solver,
      translations_Triplet_Solver,
      TrifocalTensorModel>;
  const double ThresholdUpperBound = 1.0e-2; // upper bound of the pixel residual (normalized coordinates)
  KernelType kernel(x1, x2, x3, vec_global_R_Triplet, Mat3::Identity(), ThresholdUpperBound);

  const size_t ORSA_ITER = 320;  // max number of iterations of AC-RANSAC

  TrifocalTensorModel T;
  const std::pair<double,double> acStat =
    robust::ACRANSAC(kernel, vec_inliers, ORSA_ITER, &T, dPrecision/min_focal, false);
  // If robust estimation fails => stop.
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
  for (const auto & pairIterator : map_triplet_matches)
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
  for (const uint32_t & inlier_it : vec_inliers)
  {
    tracks::STLMAPTracks::const_iterator it_tracks = tracks.begin();
    std::advance(it_tracks, inlier_it);
    const tracks::submapTrack & track = it_tracks->second;
    Observations & obs = structure[inlier_it].obs;
    for (tracks::submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
    {
      // get view Id and feat ID
      const uint32_t viewIndex = it->first;
      const uint32_t featIndex = it->second;

      // get feature
      const features::PointFeature & pt = features_provider->getFeatures(viewIndex)[featIndex];
      obs[viewIndex] = Observation(pt.coords().cast<double>(), featIndex);
    }
  }

  // Compute 3D landmark positions (triangulation of the observations)
  {
    SfM_Data_Structure_Computation_Blind structure_estimator(false);
    structure_estimator.triangulate(tiny_scene);
  }

  // Refine structure and poses (keep intrinsic constant)
  Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  if (bundle_adjustment_obj.Adjust(tiny_scene,
        Optimize_Options(
          cameras::Intrinsic_Parameter_Type::NONE,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL)
        )
      )
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

  std::ostringstream os;
  os << poses_id.i << "_" << poses_id.j << "_" << poses_id.k << ".ply";
  Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));

#endif

  // Keep the model iff it has a sufficient inlier count
  const bool bTest = ( vec_inliers.size() > 30 && 0.33 * tracks.size() );

#ifdef DEBUG_TRIPLET
  {
    std::cout << "Triplet : status: " << bTest
      << " AC: " << std::sqrt(dPrecision)
      << " inliers % " << double(vec_inliers.size()) / tracks.size() * 100.0
      << " total putative " << tracks.size() << std::endl;
  }
#endif

  return bTest;
}

} // namespace sfm
} // namespace openMVG
