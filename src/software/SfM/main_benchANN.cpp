// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Romain Janvier, Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_preemptive_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching/regions_matcher.hpp"

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using namespace openMVG;
using namespace matching;
using namespace openMVG::matching_image_collection;
using namespace openMVG::sfm;

struct Counter
{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type&) { ++count; }
  size_t count = 0;
};

template<typename T1, typename T2>
size_t intersection_size(const T1& s1, const T2& s2)
{
  Counter c;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
  return c.count;
}

int main(int argc, char **argv)
{
  CmdLine cmd;
  float fDistRatio = 0.8;

  std::string sSfM_Data_Filename;
  std::string sFeatDirectory = "";
  int max_feature_count_per_image = -1;

  //required
  cmd.add(make_option('i', sSfM_Data_Filename, "input_file"));
  cmd.add(make_option('f', sFeatDirectory, "feat_dir"));
  // optional
  cmd.add(make_option('c', max_feature_count_per_image, "max_feature_count"));


  try
  {
    if (argc == 1)
      throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string &s)
  {
    OPENMVG_LOG_ERROR << "Usage: " << argv[0] << '\n'
              << "[-i|--input_file] a SfM_Data file\n"
              << "[-f|--feat_dir path] output path where features are stored\n"
              << "--- Optional ---\n"
              << "[-c|--max_feature_count number] max_feature count per image (i.e 1000).";
    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO << " You called : "
            << "\n"
            << argv[0] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--feat_dir " << sFeatDirectory << "\n";

  if (sFeatDirectory.empty() || !stlplus::is_folder(sFeatDirectory))
  {
    OPENMVG_LOG_ERROR << "It is an invalid feature directory";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS)))
  {
    OPENMVG_LOG_ERROR <<
      "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sFeatDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    OPENMVG_LOG_ERROR << "Invalid: "
              << sImage_describer << " regions type file.";
    return EXIT_FAILURE;
  }

  system::LoggerProgress progress;

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider =
    (cmd.used('c') && max_feature_count_per_image != -1) ?
      std::make_shared<Preemptive_Regions_Provider>(max_feature_count_per_image)
      : std::make_shared<Regions_Provider>();

  if (!regions_provider->load(sfm_data, sFeatDirectory, regions_type, &progress))
  {
    OPENMVG_LOG_ERROR << "Invalid regions.";
    return EXIT_FAILURE;
  }

  // Select the pairs
  const Pair_Set pairs = exhaustivePairs(sfm_data.GetViews().size());

  // Compute matches for a reference implementation (Brute Force L2)
  // - and compare the accuracy and timing of some other method
  // - accuracy is defined as the median percentage of similar index retrieved
  const std::vector<std::string> matcher_to_evaluate = {
    "brute_force_l2",
    "hnsw_l1",
    "hnsw_l2",
    "ann_l2",
    "cascade_l2",
    "fast_cascade_l2"
  };

  OPENMVG_LOG_INFO << "Going to bench: ";
  for (const auto & method : matcher_to_evaluate)
  {
    OPENMVG_LOG_INFO << "- " << method;
  }

  struct collected_data
  {
    double time = 0.0;
    double accuracy = 0.0;
  };

  std::map<std::string, collected_data> collected_stats;
  std::unique_ptr<Matcher> collectionMatcher;
  PairWiseMatches reference_matches;
  for (const auto & method : matcher_to_evaluate)
  {
    if (method == "brute_force_l2")
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_L2));
    else if (method == "hnsw_l1")
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L1));
    else if (method == "hnsw_l2")
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L2));
    else if (method == "ann_l2")
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, ANN_L2));
    else if (method == "cascade_l2")
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, CASCADE_HASHING_L2));
    else if (method == "fast_cascade_l2")
      collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    else
      OPENMVG_LOG_ERROR << "Invalid Regions Matcher: " << method;

    collected_data data;
    if (method == "brute_force_l2") // Populate the reference matches
    {
      system::Timer timer;
      collectionMatcher->Match(regions_provider, pairs, reference_matches, &progress);
      data.time = timer.elapsed();
      collected_stats[method] = data;
    }
    else // Compute corresponding indexes and compare them against the "GT"
    {
      system::Timer timer;
      PairWiseMatches matches;
      collectionMatcher->Match(regions_provider, pairs, matches, &progress);
      data.time = timer.elapsed();

      // Compute accuracy
      // - percentage of correct retrieved index to the "GT"
      std::vector<double> accuracy_per_pair;
      for (const auto ref_matches_it : reference_matches)
      {
        const Pair pair = ref_matches_it.first;
        const auto & matches_ref = ref_matches_it.second;
        if (matches.count(pair) != 0)
        {
          const auto & matches_to_compare = matches[pair];
          const std::set<IndMatch> set_bf( matches_ref.cbegin(), matches_ref.cend());
          const std::set<IndMatch> set_compare( matches_to_compare.cbegin(), matches_to_compare.cend());
          const size_t repetability_count = intersection_size(set_bf, set_compare);
          accuracy_per_pair.emplace_back( repetability_count / static_cast<double>(set_bf.size()));
        }
      }
      double min, max, mean, median;
      minMaxMeanMedian(accuracy_per_pair.cbegin(), accuracy_per_pair.cend(),
                        min, max, mean, median);


      data.accuracy = median;
      collected_stats[method] = data;
    }
  }

  // Print results:
  for (const auto & method : matcher_to_evaluate)
  {
    OPENMVG_LOG_INFO << "Method: " << method << "\n"
      << "time(seconds): " << collected_stats[method].time;
     if (method != "brute_force_l2")
      OPENMVG_LOG_INFO << "accuracy(percent): " << collected_stats[method].accuracy;
  }

  return EXIT_SUCCESS;
}
