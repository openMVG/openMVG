
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_MATCHES_PROVIDER_HPP
#define OPENMVG_SFM_MATCHES_PROVIDER_HPP

#include <openMVG/types.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/matching/indMatch.hpp>
#include <openMVG/matching/indMatch_utils.hpp>

namespace openMVG {
namespace sfm {

/// Return the matches loaded from a provided matches file
struct Matches_Provider
{
  matching::PairWiseMatches _pairWise_matches;

  // Load matches from the provided matches file
  virtual bool load(const SfM_Data & sfm_data, const std::string & folder, const std::string & matchesMode)
  {
    if (!matching::Load(_pairWise_matches, sfm_data.GetViewsKeys(), folder, matchesMode))
    {
      std::cerr<< "Unable to read the matches file: " << folder << "/" << matchesMode << std::endl;
      return false;
    }
    return true;
  }

  /// Return the pairs used by the visibility graph defined by the pairwiser matches
  virtual Pair_Set getPairs() const
  {
    return matching::getPairs(_pairWise_matches);
  }
}; // Features_Provider

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_MATCHES_PROVIDER_HPP
