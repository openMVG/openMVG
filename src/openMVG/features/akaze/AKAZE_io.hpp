// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_IO_HPP
#define OPENMVG_FEATURES_AKAZE_IO_HPP

#include "openMVG/features/akaze/AKAZE.hpp"

#include <cereal/cereal.hpp>

template<class Archive>
void openMVG::features::AKAZE::Params::serialize(Archive & ar)
{
  ar(
    cereal::make_nvp("iNbOctave", iNbOctave),
    cereal::make_nvp("iNbSlicePerOctave", iNbSlicePerOctave),
    cereal::make_nvp("fSigma0", fSigma0),
    cereal::make_nvp("fThreshold", fThreshold),
    cereal::make_nvp("fDesc_factor", fDesc_factor));
}

#endif // OPENMVG_FEATURES_AKAZE_HPP
