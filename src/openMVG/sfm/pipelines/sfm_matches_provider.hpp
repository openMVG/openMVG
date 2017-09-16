// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_MATCHES_PROVIDER_HPP
#define OPENMVG_SFM_SFM_MATCHES_PROVIDER_HPP

#include <string>

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

namespace openMVG {
namespace sfm {

/// Return the matches loaded from a provided matches file
struct Matches_Provider
{
  matching::PairWiseMatches pairWise_matches_;

  virtual ~Matches_Provider() = default;

  // Load matches from the provided matches file
  virtual bool load(const SfM_Data & sfm_data, const std::string & matchesfile)
  {
    if (!stlplus::is_file(matchesfile))
    {
      return false;
    }
    if (!matching::Load(pairWise_matches_, matchesfile)) {
      std::cerr<< "Unable to read the matches file:" << matchesfile << std::endl;
      return false;
    }
    // Filter to keep only the one defined in SfM_Data
    {
      const Views & views = sfm_data.GetViews();
      matching::PairWiseMatches matches_saved;
      for (matching::PairWiseMatches::const_iterator iter = pairWise_matches_.begin();
        iter != pairWise_matches_.end();
        ++iter)
      {
        if (views.find(iter->first.first) != views.end() &&
          views.find(iter->first.second) != views.end())
        {
          matches_saved.insert(*iter);
        }
      }
      pairWise_matches_.swap(matches_saved);
    }
    return true;
  }

  /// Return the pairs used by the visibility graph defined by the pairwiser matches
  virtual Pair_Set getPairs() const
  {
    return matching::getPairs(pairWise_matches_);
  }
}; // Features_Provider

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_MATCHES_PROVIDER_HPP
