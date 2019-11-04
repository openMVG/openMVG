#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_2_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_2_IO_HPP

#include "openMVG/cameras/Camera_Pinhole_2.hpp"

#include <cereal/types/polymorphic.hpp>

template <class Archive>
void openMVG::cameras::Pinhole_Intrinsic_2::save( Archive & ar ) const
{
    IntrinsicBase::save(ar);
    const std::vector<double> fl {K_( 0, 0 ), K_( 1, 1 )};
    ar( cereal::make_nvp( "focal_length", fl ) );
    const std::vector<double> pp {K_( 0, 2 ), K_( 1, 2 )};
    ar( cereal::make_nvp( "principal_point", pp ) );
}


/**
* @brief  Serialization in
* @param ar Archive
*/
template <class Archive>
void openMVG::cameras::Pinhole_Intrinsic_2::load( Archive & ar )
{
    IntrinsicBase::load(ar);
    std::vector<double> focal_length( 2 );
    ar( cereal::make_nvp( "focal_length", focal_length ) );
    std::vector<double> pp( 2 );
    ar( cereal::make_nvp( "principal_point", pp ) );
    *this = Pinhole_Intrinsic_2( w_, h_, focal_length[0], focal_length[1], pp[0], pp[1] );
}

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_2, "pinhole_2");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_2);

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_2_IO_HPP
