// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_POSE3_H_
#define OPENMVG_GEOMETRY_POSE3_H_

#include "openMVG/multiview/projection.hpp"
#include <cereal/cereal.hpp> // Serialization

namespace openMVG
{
namespace geometry
{


/**
* @brief Defines a pose in 3d space
* [R|C] t = -RC
*/
class Pose3
{
  protected:

    /// Orientation matrix
    Mat3 _rotation;

    /// Center of rotation
    Vec3 _center;

  public:

    /**
    * @brief Default constructor
    * @note This defines a Null transform (aligned with cartesian frame, centered at origin)
    */
    Pose3()
      : _rotation( Mat3::Identity() ),
        _center( Vec3::Zero() )
    {

    }
    /**
    * @brief Constructor
    * @param r Rotation
    * @param c Center
    */
    Pose3( const Mat3& r, const Vec3& c ) : _rotation( r ), _center( c ) {}

    /**
    * @brief Get Rotation matrix
    * @return Rotation matrix
    */
    const Mat3& rotation() const
    {
      return _rotation;
    }

    /**
    * @brief Get Rotation matrix
    * @return Rotation matrix
    */
    Mat3& rotation()
    {
      return _rotation;
    }

    /**
    * @brief Get center of rotation
    * @return center of rotation
    */
    const Vec3& center() const
    {
      return _center;
    }

    /**
    * @brief Get center of rotation
    * @return Center of rotation
    */
    Vec3& center()
    {
      return _center;
    }

    /**
    * @brief Get translation vector
    * @return translation vector
    * @note t = -RC
    */
    inline Vec3 translation() const
    {
      return -( _rotation * _center );
    }


    /**
    * @brief Apply pose
    * @param p Point
    * @return transformed point
    */
    inline Mat3X operator () ( const Mat3X& p ) const
    {
      return _rotation * ( p.colwise() - _center );
    }


    /**
    * @brief Composition of poses
    * @param P a Pose
    * @return Composition of current pose and parameter pose
    */
    Pose3 operator * ( const Pose3& P ) const
    {
      return Pose3( _rotation * P._rotation, P._center + P._rotation.transpose() * _center );
    }


    /**
    * @brief Get inverse of the pose
    * @return Inverse of the pose
    */
    Pose3 inverse() const
    {
      return Pose3( _rotation.transpose(),  -( _rotation * _center ) );
    }


    /**
    * @brief Return the depth (distance) of a point respect to the camera center
    * @param X Input point
    * @return Distance to center
    */
    double depth( const Vec3 &X ) const
    {
      return ( _rotation * ( X - _center ) )[2];
    }

    /**
    * Serialization out
    * @param ar Archive
    */
    template <class Archive>
    void save( Archive & ar ) const
    {
      const std::vector<std::vector<double>> mat =
      {
        { _rotation( 0, 0 ), _rotation( 0, 1 ), _rotation( 0, 2 ) },
        { _rotation( 1, 0 ), _rotation( 1, 1 ), _rotation( 1, 2 ) },
        { _rotation( 2, 0 ), _rotation( 2, 1 ), _rotation( 2, 2 ) }
      };

      ar( cereal::make_nvp( "rotation", mat ) );

      const std::vector<double> vec = { _center( 0 ), _center( 1 ), _center( 2 ) };
      ar( cereal::make_nvp( "center", vec ) );
    }

    /**
    * @brief Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar )
    {
      std::vector<std::vector<double>> mat( 3, std::vector<double>( 3 ) );
      ar( cereal::make_nvp( "rotation", mat ) );
      // copy back to the rotation
      _rotation.row( 0 ) = Eigen::Map<const Vec3>( &( mat[0][0] ) );
      _rotation.row( 1 ) = Eigen::Map<const Vec3>( &( mat[1][0] ) );
      _rotation.row( 2 ) = Eigen::Map<const Vec3>( &( mat[2][0] ) );

      std::vector<double> vec( 3 );
      ar( cereal::make_nvp( "center", vec ) );
      _center = Eigen::Map<const Vec3>( &vec[0] );
    }
};
} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_POSE3_H_
