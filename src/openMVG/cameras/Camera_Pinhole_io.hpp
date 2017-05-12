// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_IO_HPP

#include "openMVG/cameras/Camera_Pinhole.hpp"

#include <cereal/types/polymorphic.hpp>

template <class Archive>
void openMVG::cameras::Pinhole_Intrinsic::save( Archive & ar ) const
{
    IntrinsicBase::save(ar);
    ar( cereal::make_nvp( "focal_length", K_( 0, 0 ) ) );
    const std::vector<double> pp = {K_( 0, 2 ), K_( 1, 2 )};
    ar( cereal::make_nvp( "principal_point", pp ) );
}


/**
* @brief  Serialization in
* @param ar Archive
*/
template <class Archive>
void openMVG::cameras::Pinhole_Intrinsic::load( Archive & ar )
{
    IntrinsicBase::load(ar);
    double focal_length;
    ar( cereal::make_nvp( "focal_length", focal_length ) );
    std::vector<double> pp( 2 );
    ar( cereal::make_nvp( "principal_point", pp ) );
    *this = Pinhole_Intrinsic( w_, h_, focal_length, pp[0], pp[1] );
}

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic, "pinhole");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic);

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_IO_HPP
