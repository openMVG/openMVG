// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_INTRINSICS_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_INTRINSICS_IO_HPP

#include <cereal/types/polymorphic.hpp>

template <class Archive>
void openMVG::cameras::IntrinsicBase::save( Archive & ar ) const
{
  ar( cereal::make_nvp( "width", w_ ) );
  ar( cereal::make_nvp( "height", h_ ) );
}

template <class Archive>
void openMVG::cameras::IntrinsicBase::load( Archive & ar )
{
  ar( cereal::make_nvp( "width", w_ ) );
  ar( cereal::make_nvp( "height", h_ ) );
}

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_INTRINSICS_IO_HPP
