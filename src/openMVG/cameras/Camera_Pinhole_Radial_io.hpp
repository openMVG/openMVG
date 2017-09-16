// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_RADIAL_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_RADIAL_IO_HPP

#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Radial_K1::save( Archive & ar ) const
{
    Pinhole_Intrinsic::save( ar );
    ar( cereal::make_nvp( "disto_k1", params_ ) );
}

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Radial_K1::load( Archive & ar )
{
    Pinhole_Intrinsic::load(ar);
    ar( cereal::make_nvp( "disto_k1", params_ ) );
}

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Radial_K3::save( Archive & ar ) const
{
    Pinhole_Intrinsic::save(ar);
    ar( cereal::make_nvp( "disto_k3", params_ ) );
}

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Radial_K3::load( Archive & ar )
{
    Pinhole_Intrinsic::load(ar);
    ar( cereal::make_nvp( "disto_k3", params_ ) );
}

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_Radial_K1, "pinhole_radial_k1");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Radial_K1);
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_Radial_K3, "pinhole_radial_k3");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Radial_K3);

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_RADIAL_IO_HPP
