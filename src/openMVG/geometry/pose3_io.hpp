// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_POSE3_IO_HPP
#define OPENMVG_GEOMETRY_POSE3_IO_HPP

#include <cereal/cereal.hpp> // Serialization

#include <vector>

template <class Archive>
inline void openMVG::geometry::Pose3::save( Archive & ar ) const
{
    const std::vector<std::vector<double>> mat =
    {
    { rotation_( 0, 0 ), rotation_( 0, 1 ), rotation_( 0, 2 ) },
    { rotation_( 1, 0 ), rotation_( 1, 1 ), rotation_( 1, 2 ) },
    { rotation_( 2, 0 ), rotation_( 2, 1 ), rotation_( 2, 2 ) }
    };

    ar( cereal::make_nvp( "rotation", mat ) );

    const std::vector<double> vec = { center_( 0 ), center_( 1 ), center_( 2 ) };
    ar( cereal::make_nvp( "center", vec ) );
}

template <class Archive>
inline void openMVG::geometry::Pose3::load( Archive & ar )
{
    std::vector<std::vector<double>> mat( 3, std::vector<double>( 3 ) );
    ar( cereal::make_nvp( "rotation", mat ) );
    // copy back to the rotation
    rotation_.row( 0 ) = Eigen::Map<const Vec3>( &( mat[0][0] ) );
    rotation_.row( 1 ) = Eigen::Map<const Vec3>( &( mat[1][0] ) );
    rotation_.row( 2 ) = Eigen::Map<const Vec3>( &( mat[2][0] ) );

    std::vector<double> vec( 3 );
    ar( cereal::make_nvp( "center", vec ) );
    center_ = Eigen::Map<const Vec3>( &vec[0] );
}

#endif  // OPENMVG_GEOMETRY_POSE3_IO_HPP