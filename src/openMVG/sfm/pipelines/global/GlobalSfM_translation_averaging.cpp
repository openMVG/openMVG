
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
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
namespace globalSfM{

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
  // Assert relative matches only share index between the relative poseId defined by the global rotations
  matching::PairWiseMatches map_Matches_Ecpy = matches_provider->_pairWise_matches;
  std::set<IndexT> set_remainingIds;
  std::transform(map_globalR.begin(), map_globalR.end(), std::inserter(set_remainingIds, set_remainingIds.begin()), RetrieveKey());
  KeepOnlyReferencedElement(set_remainingIds, map_Matches_Ecpy);

  // Compute the translations and save them to vec_initialRijTijEstimates:
  Compute_translations(
    sfm_data,
    normalized_features_provider,
    map_Matches_Ecpy,
    map_globalR,
    tripletWise_matches);

  // Keep the largest Biedge connected component graph of relative translations
  Pair_Set pairs;
  std::transform(vec_initialRijTijEstimates.begin(), vec_initialRijTijEstimates.end(), std::inserter(pairs, pairs.begin()), std::RetrieveKey());
  set_remainingIds = openMVG::graphUtils::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs, "./");
  KeepOnlyReferencedElement(set_remainingIds, map_Matches_Ecpy);
  KeepOnlyReferencedElement(set_remainingIds, vec_initialRijTijEstimates);
  KeepOnlyReferencedElement(set_remainingIds, tripletWise_matches);

  return Translation_averaging(
    eTranslationAveragingMethod,
    sfm_data,
    normalized_features_provider,
    map_Matches_Ecpy,
    map_globalR);
}

bool GlobalSfM_Translation_AveragingSolver::Translation_averaging(
  ETranslationAveragingMethod eTranslationAveragingMethod,
  SfM_Data & sfm_data,
  const Features_Provider * normalized_features_provider,
  const matching::PairWiseMatches & map_Matches_E,
  const Hash_Map<IndexT, Mat3> & map_globalR)
{
  //-------------------
  //-- GLOBAL TRANSLATIONS ESTIMATION from initial triplets t_ij guess
  //-------------------
  const std::string _sOutDirectory("./");
  {
    const std::set<IndexT> index = getIndexT(vec_initialRijTijEstimates);

    const size_t iNview = index.size();
    std::cout << "\n-------------------------------" << "\n"
      << " Global translations computation: " << "\n"
      << "   - Ready to compute " << iNview << " global translations." << "\n"
      << "     from #relative translations: " << vec_initialRijTijEstimates.size() << std::endl;

    if (iNview < 3)
    {
      // Too tiny image set to perform motion averaging
      return false;
    }
    //-- Update initial estimates from [minId,maxId] to range [0->Ncam]
    RelativeInfo_Vec vec_initialRijTijEstimates_cpy = vec_initialRijTijEstimates;
    const Pair_Set pairs = getPairs(vec_initialRijTijEstimates_cpy);
    Hash_Map<IndexT,IndexT> _reindexForward, _reindexBackward;
    openMVG::reindex(pairs, _reindexForward, _reindexBackward);
    for(size_t i = 0; i < vec_initialRijTijEstimates_cpy.size(); ++i)
    {
      openMVG::relativeInfo & rel = vec_initialRijTijEstimates_cpy[i];
      rel.first = Pair(_reindexForward[rel.first.first], _reindexForward[rel.first.second]);
    }

    openMVG::Timer timerLP_translation;

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
        //if (_bHtmlReport)
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
          const IndexT camNodeId = _reindexBackward[i];
          const View * view = sfm_data.views[camNodeId].get();
          const Mat3 & Ri = map_globalR.at(view->id_pose);
          sfm_data.poses[view->id_pose] = Pose3(Ri, -Ri.transpose()*t);
        }
      }
      break;

      case TRANSLATION_AVERAGING_L2:
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
        if(!solve_translations_problem(
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
          const IndexT camNodeId = _reindexBackward[i]; // undo the reindexing
          const View * view = sfm_data.views.at(camNodeId).get();
          const Mat3 & Ri = map_globalR.at(view->id_pose);
          sfm_data.poses[view->id_pose] = Pose3(Ri, C);
        }
      }
      break;
    }
  }
  return true;
}

void GlobalSfM_Translation_AveragingSolver::Compute_translations(
  const SfM_Data & sfm_data,
  const Features_Provider * normalized_features_provider,
  const matching::PairWiseMatches & map_Matches_E,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  matching::PairWiseMatches &tripletWise_matches)
{
  std::cout << "\n-------------------------------" << "\n"
    << " Relative translations computation: " << "\n"
    << "-------------------------------" << std::endl;

  // List putative triplets
  const Pair_Set pairs = getPairs(map_Matches_E);
  std::vector< graphUtils::Triplet > vec_triplets = graphUtils::tripletListing(pairs);

  std::cout << "#Triplets: " << vec_triplets.size() << std::endl;

  // Compute putative translations with an edge coverage algorithm

  openMVG::Timer timerLP_triplet;
  // Compute relative translations over the graph of putative triplets

  ComputePutativeTranslation_EdgesCoverage(
    sfm_data,
    map_globalR,
    normalized_features_provider,
    map_Matches_E,
    vec_triplets,
    vec_initialRijTijEstimates,
    tripletWise_matches);
  const double timeLP_triplet = timerLP_triplet.elapsed();
  std::cout << "TRIPLET COVERAGE TIMING: " << timeLP_triplet << " seconds" << std::endl;

  std::cout << "-------------------------------" << "\n"
      << "-- #Effective translations estimates: " << vec_initialRijTijEstimates.size()/3
      << " from " <<vec_triplets.size() << " triplets.\n"
      << "-- resulting in " <<vec_initialRijTijEstimates.size() << " translation estimation.\n"
      << "-- timing to obtain the relative translations: " << timeLP_triplet << " seconds.\n"
      << "-------------------------------" << std::endl;
}

//-- Perform a trifocal estimation of the graph contain in vec_triplets with an
// edge coverage algorithm. It's complexity is sub-linear in term of edges count.
void GlobalSfM_Translation_AveragingSolver::ComputePutativeTranslation_EdgesCoverage(
  const SfM_Data & sfm_data,
  const Hash_Map<IndexT, Mat3> & map_globalR,
  const Features_Provider * normalized_features_provider,
  const matching::PairWiseMatches & map_Matches_E,
  const std::vector< graphUtils::Triplet > & vec_triplets,
  RelativeInfo_Vec & vec_initialEstimates,
  matching::PairWiseMatches & newpairMatches)
{
  //-- Prepare tracks count per triplets:
  Hash_Map<IndexT, IndexT> map_tracksPerTriplets;
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < (int)vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const IndexT I = triplet.i, J = triplet.j , K = triplet.k;

    PairWiseMatches map_matchesIJK;
    if(map_Matches_E.find(std::make_pair(I,J)) != map_Matches_E.end())
      map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(I,J)));
    else
    if(map_Matches_E.find(std::make_pair(J,I)) != map_Matches_E.end())
      map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(J,I)));

    if(map_Matches_E.find(std::make_pair(I,K)) != map_Matches_E.end())
      map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(I,K)));
    else
    if(map_Matches_E.find(std::make_pair(K,I)) != map_Matches_E.end())
      map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(K,I)));

    if(map_Matches_E.find(std::make_pair(J,K)) != map_Matches_E.end())
      map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(J,K)));
    else
    if(map_Matches_E.find(std::make_pair(K,J)) != map_Matches_E.end())
      map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(K,J)));

    // Compute tracks:
    openMVG::tracks::STLMAPTracks map_tracks;
    openMVG::tracks::TracksBuilder tracksBuilder;
    {
      tracksBuilder.Build(map_matchesIJK);
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
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const IndexT I = triplet.i, J = triplet.j , K = triplet.k;
    // Add three edges
    set_edges.insert(std::make_pair(std::min(I,J), std::max(I,J)));
    set_edges.insert(std::make_pair(std::min(I,K), std::max(I,K)));
    set_edges.insert(std::make_pair(std::min(J,K), std::max(J,K)));
  }

  // Copy them in vector in order to try to compute them in parallel
  std::vector<myEdge > vec_edges(set_edges.begin(), set_edges.end());

  MutexSet<myEdge> m_mutexSet;

  C_Progress_display my_progress_bar(
    vec_edges.size(),
    std::cout,
    "\nComputation of the relative translations over the graph with an edge coverage algorithm\n");

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
    // Find the triplet that contain the given edge
    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graphUtils::Triplet & triplet = vec_triplets[i];
      if (triplet.contain(edge))
      {
        vec_possibleTriplets.push_back(i);
      }
    }

    //-- Sort the triplet according the number of matches they have on their edges
    std::vector<size_t> vec_commonTracksPerTriplets(vec_possibleTriplets.size(), 0);
    for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
    {
      vec_commonTracksPerTriplets[i] = map_tracksPerTriplets[vec_possibleTriplets[i]];
    }
    //-- If current edge already computed continue
    if (m_mutexSet.count(edge))
      continue;

    using namespace indexed_sort;
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
      const graphUtils::Triplet & triplet = vec_triplets[vec_possibleTriplets[i]];
      const IndexT I = triplet.i, J = triplet.j , K = triplet.k;
      {
        PairWiseMatches map_matchesIJK;
        if(map_Matches_E.find(std::make_pair(I,J)) != map_Matches_E.end())
          map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(I,J)));
        else
        if(map_Matches_E.find(std::make_pair(J,I)) != map_Matches_E.end())
          map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(J,I)));

        if(map_Matches_E.find(std::make_pair(I,K)) != map_Matches_E.end())
          map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(I,K)));
        else
        if(map_Matches_E.find(std::make_pair(K,I)) != map_Matches_E.end())
          map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(K,I)));

        if(map_Matches_E.find(std::make_pair(J,K)) != map_Matches_E.end())
          map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(J,K)));
        else
        if(map_Matches_E.find(std::make_pair(K,J)) != map_Matches_E.end())
          map_matchesIJK.insert(*map_Matches_E.find(std::make_pair(K,J)));

        // Select common point:
        openMVG::tracks::STLMAPTracks map_tracksCommon;
        openMVG::tracks::TracksBuilder tracksBuilder;
        {
          tracksBuilder.Build(map_matchesIJK);
          tracksBuilder.Filter(3);
          tracksBuilder.ExportToSTL(map_tracksCommon);
        }

        //--
        // Try to estimate this triplet.
        //--

        // Get rotations:
        std::vector<Mat3> vec_global_KR_Triplet(3);
        vec_global_KR_Triplet[0] = (map_globalR.at(I));
        vec_global_KR_Triplet[1] = (map_globalR.at(J));
        vec_global_KR_Triplet[2] = (map_globalR.at(K));

        const View * view_I = sfm_data.GetViews().at(I).get();
        const View * view_J = sfm_data.GetViews().at(J).get();
        const View * view_K = sfm_data.GetViews().at(K).get();

        // update pixel precision (from image to camera coordinate system, since features are normalized)
        const Pinhole_Intrinsic* cam_I = dynamic_cast<Pinhole_Intrinsic*>(sfm_data.GetIntrinsics().at(view_I->id_intrinsic).get());
        const Pinhole_Intrinsic* cam_J = dynamic_cast<Pinhole_Intrinsic*>(sfm_data.GetIntrinsics().at(view_J->id_intrinsic).get());
        const Pinhole_Intrinsic* cam_K = dynamic_cast<Pinhole_Intrinsic*>(sfm_data.GetIntrinsics().at(view_K->id_intrinsic).get());
        if (cam_I == NULL || cam_J == NULL || cam_K == NULL)
        {
          continue;
        }

        const double averageFocal = ( cam_I->focal() + cam_J->focal() + cam_K->focal() ) / 3.0 ;

        double dPrecision = 4.0 / averageFocal / averageFocal;
        const double ThresholdUpperBound = 0.5 / averageFocal;

        std::vector<Vec3> vec_tis(3);
        std::vector<size_t> vec_inliers;

        std::string _sOutDirectory = "./";
        if (map_tracksCommon.size() > 50 &&
            Estimate_T_triplet(
              map_tracksCommon, normalized_features_provider,  vec_global_KR_Triplet,
              vec_tis, dPrecision, vec_inliers, ThresholdUpperBound,
              I, J, K, _sOutDirectory))
        {
          if (bVerbose)
            std::cout << I << " " << J << " " << K << ": "
              << dPrecision * averageFocal << "\t" << vec_inliers.size() << std::endl;

          //-- Build the three camera:
          const Mat3 RI = map_globalR.find(I)->second;
          const Mat3 RJ = map_globalR.find(J)->second;
          const Mat3 RK = map_globalR.find(K)->second;
          const Vec3 ti = vec_tis[0];
          const Vec3 tj = vec_tis[1];
          const Vec3 tk = vec_tis[2];

          // Build the 3 relative translations estimations.
          // IJ, JK, IK

//--- ATOMIC
#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
          {
            Mat3 RijGt;
            Vec3 tij;
            RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
            vec_initialEstimates.push_back(
              std::make_pair(std::make_pair(I, J), std::make_pair(RijGt, tij)));

            Mat3 RjkGt;
            Vec3 tjk;
            RelativeCameraMotion(RJ, tj, RK, tk, &RjkGt, &tjk);
            vec_initialEstimates.push_back(
              std::make_pair(std::make_pair(J, K), std::make_pair(RjkGt, tjk)));

            Mat3 RikGt;
            Vec3 tik;
            RelativeCameraMotion(RI, ti, RK, tk, &RikGt, &tik);
            vec_initialEstimates.push_back(
              std::make_pair(std::make_pair(I, K), std::make_pair(RikGt, tik)));

            // Add trifocal inliers as valid 3D points
            for (std::vector<size_t>::const_iterator iterInliers = vec_inliers.begin();
              iterInliers != vec_inliers.end(); ++iterInliers)
            {
              openMVG::tracks::STLMAPTracks::const_iterator iterTracks = map_tracksCommon.begin();
              std::advance(iterTracks, *iterInliers);
              const openMVG::tracks::submapTrack & subTrack = iterTracks->second;
              openMVG::tracks::submapTrack::const_iterator iterI, iterJ, iterK;
              iterI = iterJ = iterK = subTrack.begin();
              std::advance(iterJ,1);
              std::advance(iterK,2);
              /*openMVG::tracks::submapTrack::const_iterator
                iterI = subTrack.find(I),
                iterJ = subTrack.find(J),
                iterK = subTrack.find(K);*/

              newpairMatches[std::make_pair(I,J)].push_back(IndMatch(iterI->second, iterJ->second));
              newpairMatches[std::make_pair(J,K)].push_back(IndMatch(iterJ->second, iterK->second));
              newpairMatches[std::make_pair(I,K)].push_back(IndMatch(iterI->second, iterK->second));
            }
          }
          //-- Set the 3 edges as already computed by the "trifocal tensor" (translations triplet)
          m_mutexSet.insert(std::make_pair(std::min(I,J), std::max(I,J)));
          m_mutexSet.insert(std::make_pair(std::min(I,K), std::max(I,K)));
          m_mutexSet.insert(std::make_pair(std::min(J,K), std::max(J,K)));
          break;
        }
      }
    }
  }
}

// Robust estimation and refinement of a translation and 3D points of an image triplets.
bool GlobalSfM_Translation_AveragingSolver::Estimate_T_triplet(
  const openMVG::tracks::STLMAPTracks & map_tracksCommon,
  const Features_Provider * features_provider,
  const std::vector<Mat3> & vec_global_R_Triplet,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<size_t> & vec_inliers,
  const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
  const size_t nI,
  const size_t nJ,
  const size_t nK,
  const std::string & sOutDirectory) const
{
  using namespace linearProgramming;
  using namespace lInfinityCV;

  // Convert data
  Mat x1(2, map_tracksCommon.size());
  Mat x2(2, map_tracksCommon.size());
  Mat x3(2, map_tracksCommon.size());

  Mat* xxx[3] = {&x1, &x2, &x3};

  size_t cpt = 0;
  for (tracks::STLMAPTracks::const_iterator iterTracks = map_tracksCommon.begin();
    iterTracks != map_tracksCommon.end(); ++iterTracks, ++cpt) {
    const tracks::submapTrack & subTrack = iterTracks->second;
    size_t index = 0;
    for (tracks::submapTrack::const_iterator iter = subTrack.begin(); iter != subTrack.end(); ++iter, ++index) {
      const size_t imaIndex = iter->first;
      const size_t featIndex = iter->second;
      const PointFeature pt = features_provider->getFeatures(imaIndex)[featIndex];
      xxx[index]->col(cpt) = pt.coords().cast<double>();
    }
  }

  using namespace openMVG::trifocal;
  using namespace openMVG::trifocal::kernel;

  typedef TranslationTripletKernel_ACRansac<
    translations_Triplet_Solver,
    translations_Triplet_Solver,
    TrifocalTensorModel> KernelType;
  KernelType kernel(x1, x2, x3, vec_global_R_Triplet, Mat3::Identity(), ThresholdUpperBound);

  const size_t ORSA_ITER = 320;

  TrifocalTensorModel T;
  // dPrecision = std::numeric_limits<double>::infinity();
  std::pair<double,double> acStat = robust::ACRANSAC(kernel, vec_inliers, ORSA_ITER, &T, dPrecision, false);
  dPrecision = acStat.first;

  //-- Export data in order to have an idea of the precision of the estimates
  vec_tis.resize(3);
  Mat3 K, R;
  KRt_From_P(T.P1, &K, &R, &vec_tis[0]);
  KRt_From_P(T.P2, &K, &R, &vec_tis[1]);
  KRt_From_P(T.P3, &K, &R, &vec_tis[2]);

  bool bTest(vec_inliers.size() > 30);
  //if (!bTest)
  //{
  //  std::cout << "Triplet rejected : AC: " << dPrecision
  //    << " median: " << median
  //    << " #inliers " << vec_inliers.size()
  //    << " #putatives " << map_tracksCommon.size() << std::endl;
  //}
  return bTest;
}


} // namespace globalSfM
} // namespace openMVG
