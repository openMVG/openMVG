// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romain Janvier <romain.janvier~AT~univ-orleans.fr> for the given adaptation

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_FISHEYE_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_FISHEYE_IO_HPP

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Fisheye::save( Archive & ar ) const
{
    ar(cereal::base_class<Pinhole_Intrinsic>(this));
    ar( cereal::make_nvp( "fisheye", params_ ) );
}

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Fisheye::load( Archive & ar )
{
    ar(cereal::base_class<Pinhole_Intrinsic>(this));
    ar( cereal::make_nvp( "fisheye", params_ ) );
}

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Fisheye, "fisheye" );
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Fisheye)

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_FISHEYE_IO_HPP
