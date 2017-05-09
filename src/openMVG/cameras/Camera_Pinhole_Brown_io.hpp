// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Sida Li, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_IO_HPP

#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Brown_T2::save( Archive & ar ) const
{
    ar(cereal::base_class<Pinhole_Intrinsic>(this));
    ar( cereal::make_nvp( "disto_t2", params_ ) );
}

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Brown_T2::load( Archive & ar )
{
    ar(cereal::base_class<Pinhole_Intrinsic>(this));
    ar( cereal::make_nvp( "disto_t2", params_ ) );
}

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Brown_T2, "pinhole_brown_t2" );
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Brown_T2)

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_IO_HPP

