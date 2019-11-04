#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_2_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_2_IO_HPP

#include "openMVG/cameras/Camera_Pinhole_Brown_2.hpp"

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Brown_T2_2::save( Archive & ar ) const
{
    ar(cereal::base_class<Pinhole_Intrinsic_2>(this));
    ar( cereal::make_nvp( "disto_t2_2", params_ ) );
}

template <class Archive>
inline void openMVG::cameras::Pinhole_Intrinsic_Brown_T2_2::load( Archive & ar )
{
    ar(cereal::base_class<Pinhole_Intrinsic_2>(this));
    ar( cereal::make_nvp( "disto_t2_2", params_ ) );
}

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Brown_T2_2, "pinhole_brown_t2_2" );
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Brown_T2_2)

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_BROWN_2_IO_HPP

