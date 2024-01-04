// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching_image_collection/Retrieval_Helpers.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::retrieval;

int main(int argc, char **argv) {
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sGTMatchFilename;
  std::string sInputPairFilename;
  std::string sInputSimFile;

  // required
  cmd.add(make_option('i', sSfM_Data_Filename, "input_file"));
  cmd.add(make_option('p', sInputPairFilename, "pair_file"));
  cmd.add(make_option('s', sInputSimFile, "sim_file"));
  cmd.add(make_option('g', sGTMatchFilename, "gt_matches"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string &s) {
    std::cerr
        << "Usage: " << argv[0] << '\n'
        << "[-i|--input_file] a SfM_Data file containing the GT "
           "reconstruction\n"
        << "[-s|--sim_file] use a pairwise similarity file to compute mAP\n"
        << "\n[Optional]\n"
        << "[-p|--pair_file] alternatively to the sim file, the pair file you "
           "want to compare to your GT pairs (mAP won't be computed)\n"
        << "[-g|--gt_matches] a matches.bin/.txt file (i.e "
           "matches.f.txt), else pairs from the input sfm_data are used\n"
        << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << " You called : "
            << "\n"
            << argv[0] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--sim_file " << sInputSimFile << "\n"
            << "--pair_file " << sInputPairFilename << "\n"
            << "--gt_matches " << sGTMatchFilename << "\n"
            << std::endl;

  // Are the gt pairs loaded from matches or are they extracted from sfm_data?
  const bool gt_from_matches = !sGTMatchFilename.empty();
  // Are the putative pairs loaded from sim file or from a pair file
  const bool putative_from_sim = !sInputSimFile.empty();

  if (gt_from_matches && !stlplus::file_exists(sGTMatchFilename)) {
    std::cerr << "It is an invalid --gt_matches file" << std::endl;
    return EXIT_FAILURE;
  }

  if (putative_from_sim && !stlplus::file_exists(sInputSimFile)) {
    std::cerr << "It is an invalid --sim_file similarity file" << std::endl;
    return EXIT_FAILURE;
  }

  if (!putative_from_sim && (sInputPairFilename.empty() ||
                             !stlplus::file_exists(sInputPairFilename))) {
    std::cerr << "It is an invalid --pair_file pair file" << std::endl;
    return EXIT_FAILURE;
  }

  PairWiseMatches gt_pairwise_matches;

  //---------------------------------------
  // Read SfM Scene (image view)
  //---------------------------------------
  ESfM_Data sfm_data_flags =
      gt_from_matches ? ESfM_Data(VIEWS) : ESfM_Data(VIEWS | STRUCTURE);

  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(sfm_data_flags))) {
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfM_Data_Filename
              << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  Pair_Set gt_pairs;
  if (gt_from_matches) {
    if (!(Load(gt_pairwise_matches, sGTMatchFilename))) {
      std::cerr << "Cannot load input GT matches file";
      return EXIT_FAILURE;
    }
    std::cout << "\t GT RESULTS LOADED;"
              << " #pairs: " << gt_pairwise_matches.size() << std::endl;
    gt_pairs = getPairs(gt_pairwise_matches);
    // clear some memory
    gt_pairwise_matches.clear();
  } else {
    for (auto &track : sfm_data.GetLandmarks()) {
      const auto &obs = track.second.obs;
      for (auto obs_it1 = obs.cbegin(); obs_it1 != obs.cend(); ++obs_it1) {
        for (auto obs_it2 = std::next(obs_it1, 1); obs_it2 != obs.cend();
             ++obs_it2) {
          gt_pairs.insert({std::min(obs_it1->first, obs_it2->first),
                           std::max(obs_it1->first, obs_it2->first)});
        }
      }
    }
  }

  Pair_Set putative_pairs;
  IndexedPairwiseSimilarity<std::greater<double>> ordered_results;
  if (putative_from_sim) {
    if (!loadPairwiseSimilarityScores(sInputSimFile, ordered_results)) {
      std::cerr << "Cannot load input sim file";
      return EXIT_FAILURE;
    }
    for (const auto &indexed_pairwise_score : ordered_results) {
      IndexT i = indexed_pairwise_score.first;
      for (const auto &tuple_score : indexed_pairwise_score.second) {
        const IndexT j = tuple_score.second;
        putative_pairs.insert({std::min(i, j), std::max(i, j)});
      }
    }
  } else {
    if (!loadPairs(sfm_data.GetViews().size(), sInputPairFilename,
                   putative_pairs)) {
      std::cerr << "Cannot load input pair file";
      return EXIT_FAILURE;
    }
  }
  std::cout << "\t PUTATIVE RESULTS LOADED;"
            << " #pairs: " << putative_pairs.size() << std::endl;

  // Computing retrieval statistics:
  // See https://en.wikipedia.org/wiki/Sensitivity_and_specificity
  //
  // TP: Pair count that must be found
  // FP: Pair count that are found that should not be found
  std::set<Pair> common_pairs;
  std::set_intersection(gt_pairs.cbegin(), gt_pairs.cend(),
                        putative_pairs.cbegin(), putative_pairs.cend(),
                        std::inserter(common_pairs, common_pairs.begin()));
  // TP: correctly identified
  const int true_positive_count = common_pairs.size();
  // FN: missed - not detected
  const int false_negative_count = gt_pairs.size() - common_pairs.size();
  // FP: detected by we should not have
  const int false_positive_count = putative_pairs.size() - common_pairs.size();
  // TN: correctly rejected
  // const int true_negative_count =

  // True positive = correctly identified
  // False positive = incorrectly identified
  // True negative = correctly rejected
  // False negative = incorrectly rejected

  // Recall, same as TPR(True Positive Rate)
  const double recall =
      true_positive_count /
      static_cast<float>(true_positive_count + false_negative_count);
  const double precision =
      true_positive_count /
      static_cast<float>(true_positive_count + false_positive_count);

  std::cout << "Retrieval statistics:\n"
            << "\t#GT pairs: " << gt_pairs.size() << "\n"
            << "\t#Putative pairs: " << putative_pairs.size() << "\n"
            << "\n"
            << "\tTrue positive count: " << true_positive_count << "\n"
            << "\tFalse negative count: " << false_negative_count << "\n"
            << "\n"
            //<< "\tFalse positive count: " << false_positive_count << "\n"
            << "\tRecall: " << recall << "\n"
            << "\tprecision: " << precision << "\n"
            << std::endl;

  // mAP computation
  if (putative_from_sim) {
    double average_precision = 0.0;

    for (const auto &indexed_pairwise_score : ordered_results) {
      const IndexT i = indexed_pairwise_score.first;
      uint32_t counter = 0;
      uint32_t counter_TP = 0;
      double cur_AP = 0.;
      for (const auto &tuple_score : indexed_pairwise_score.second) {
        const IndexT j = tuple_score.second;
        bool is_GT = (gt_pairs.count({std::min(i, j), std::max(i, j)}) > 0);
        counter++;
        if (is_GT) {
          counter_TP++;
          cur_AP += counter_TP / static_cast<double>(counter);
        }
      }
      uint32_t local_num_gt = 0;
      for (const auto &pair : gt_pairs) {  // TODO: C++17 reduce
        local_num_gt = (pair.first == i || pair.second == i) ? local_num_gt + 1
                                                             : local_num_gt;
      }
      cur_AP = cur_AP / local_num_gt;
      average_precision += cur_AP;
    }

    double mAP = average_precision / ordered_results.size();
    std::cout << "\tmAP: " << mAP << std::endl;
  }
  return EXIT_SUCCESS;
}