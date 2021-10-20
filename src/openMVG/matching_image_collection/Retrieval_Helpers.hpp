// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.
// Copyright (c) 2020 Romain JANVIER.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_RETRIEVAL_HELPERS_HPP
#define OPENMVG_RETRIEVAL_HELPERS_HPP

#include <string>
namespace openMVG {
namespace retrieval {

template <class Order>
using IndexedPairwiseSimilarity =
    std::map<IndexT, std::map<double, IndexT, Order>>;

template <class Order>
inline bool savePairwiseSimilarityScores(
    const std::string &sFileName,
    const IndexedPairwiseSimilarity<Order> &simScores) {
  std::ofstream outStream(sFileName.c_str());
  if (!outStream.is_open()) {
    std::cerr << std::endl
              << "Pairwise Similarity: Impossible to open the output specified "
                 "file: \""
              << sFileName << "\"." << std::endl;
    return false;
  }
  for (const auto &cur_view : simScores) {
    for (const auto &cur_pairwise_scores : cur_view.second) {
      outStream << cur_view.first << ' ' << cur_pairwise_scores.second << ' '
                << cur_pairwise_scores.first << '\n';
    }
  }
  const bool bOk = !outStream.bad();
  outStream.close();
  return bOk;
}

template <class Order>
inline bool loadPairwiseSimilarityScores(
    const std::string &sFileName, IndexedPairwiseSimilarity<Order> &simScores) {
  std::ifstream inStream(sFileName.c_str());
  if (!inStream.is_open()) {
    std::cerr << std::endl
              << "Pairwise Similarity: Impossible to open the input specified "
                 "file: \""
              << sFileName << "\"." << std::endl;
    return false;
  }

  std::string sValue;
  std::vector<std::string> vec_str;
  while (std::getline(inStream, sValue)) {
    stl::split(sValue, ' ', vec_str);
    if (vec_str.size() != 3) {
      std::cerr << "Pairwise Similarity: Invalid input file: \"" << sFileName
                << "\"." << std::endl;
      return false;
    }

    IndexT i = 0, j = 0;
    double score = 0.;
    try {
      i = static_cast<IndexT>(std::stoul(vec_str[0]));
      j = std::stoul(vec_str[1]);
      score = std::stod(vec_str[2]);
    } catch (const std::exception &exp) {
      return false;
    }

    if (i == j) continue;

    simScores[i].insert({score, j});
  }

  inStream.close();
  return true;
}

}  // namespace retrieval
}  // namespace openMVG

#endif  // OPENMVG_RETRIEVAL_HELPERS_HPP
