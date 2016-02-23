// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_SIMILARITY3_H_
#define OPENMVG_GEOMETRY_SIMILARITY3_H_


namespace openMVG
{
namespace geometry
{

/**
* @brief Define a 3D Similarity transform encoded as a 3D pose plus a scale
*/
struct Similarity3
{
  /// Pose
  Pose3 _pose;

  /// Scale
  double _scale;

  /**
  * @brief Default constructor
  * @note This define identity transformation centered at origin of a cartesian frame
  */
  Similarity3()
    : _pose( Pose3() ),
      _scale( 1.0 )
  {

  }

  /**
  * @brief Constructor
  * @param pose a 3d pose
  * @param scale a scale factor
  */
  Similarity3( const Pose3 & pose, const double scale )
    : _pose( pose ),
      _scale( scale )
  {

  }

  /**
  * @brief Apply transformation to a point
  * @param point Input point
  * @return transformed point
  */
  Vec3 operator () ( const Vec3 & point ) const
  {
    return _scale * _pose( point );
  }

  /**
  * @brief Concatenation of pose
  * @param pose Pose to be concatenated with the current one
  * @return Concatenation of poses
  */
  Pose3 operator () ( const Pose3 & pose ) const
  {
    return Pose3( pose.rotation() * _pose.rotation().transpose(), this->operator()( pose.center() ) );
  }
};

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_SIMILARITY3_H_
