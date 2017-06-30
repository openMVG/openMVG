// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_IO_HPP
#define OPENMVG_MATCHING_IND_MATCH_IO_HPP

#include "openMVG/matching/indMatch.hpp"

#include <cereal/cereal.hpp> // Serialization

// Serialization
template <class Archive>
void openMVG::matching::IndMatch::serialize( Archive & ar )  {
  ar(i_, j_);
}

namespace cereal
{
  // This struct specialization will tell cereal which is the right way to serialize PairWiseMatches
  template <class Archive>
  struct specialize<
    Archive,
    openMVG::matching::PairWiseMatches,
    cereal::specialization::member_serialize>
  {};
} // namespace cereal

#endif // OPENMVG_MATCHING_IND_MATCH_IO_HPP