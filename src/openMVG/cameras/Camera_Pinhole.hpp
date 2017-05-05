// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_HPP

#include "openMVG/cameras/Camera_Intrinsics.hpp"

namespace openMVG
{
namespace cameras
{

/**
* @brief Define an ideal Pinhole camera intrinsics (store a K 3x3 matrix)
* with intrinsic parameters defining the K calibration matrix
*
* Intrinsic camera matrix is \f$ K = \begin{pmatrix} f & 0 & u_0 \\ 0 & f & v_0 \\ 0 & 0 & 1 \end{pmatrix} \f$
*
* @note This is an ideal Pinhole camera because it doesn't handle skew and distortion
* @note The camera does only handle one focal length (ie: \f$ f_x = f_y = f \f$ )
*/
class Pinhole_Intrinsic : public IntrinsicBase
{
  using class_type = Pinhole_Intrinsic;

  protected:

    /// Intrinsic matrix : Focal & principal point are embed into the calibration matrix K
    Mat3 K_;

    /// Inverse of intrinsic matrix
    Mat3 Kinv_;

  public:

    /**
    * @brief Constructor
    * @param w Width of the image plane
    * @param h Height of the image plane
    * @param focal_length_pix Focal length (in pixel) of the camera
    * @param ppx Principal point on x-axis
    * @param ppy Principal point on y-axis
    */
    Pinhole_Intrinsic(
      unsigned int w = 0, unsigned int h = 0,
      double focal_length_pix = 0.0,
      double ppx = 0.0, double ppy = 0.0 );

    /**
    * @brief Constructor
    * @param w Width of the image plane
    * @param h Height of the image plane
    * @param K Intrinsic Matrix (3x3) {f,0,ppx; 0,f,ppy; 0,0,1}
    */
    Pinhole_Intrinsic(
      unsigned int w,
      unsigned int h,
      const Mat3& K);

    /**
    * @brief Destructor
    */
    ~Pinhole_Intrinsic() override = default;

    /**
    * @brief Get type of the intrinsic
    * @retval PINHOLE_CAMERA
    */
    EINTRINSIC getType() const override;

    /**
    * @brief Get the intrinsic matrix
    * @return 3x3 intrinsic matrix
    */
    const Mat3& K() const;

    /**
    * @brief Get the inverse of the intrinsic matrix
    * @return Inverse of intrinsic matrix
    */
    const Mat3& Kinv() const;


    /**
    * @brief Return the value of the focal in pixels
    * @return Focal of the camera (in pixel)
    */
    double focal() const;

    /**
    * @brief Get principal point of the camera
    * @return Principal point of the camera
    */
    Vec2 principal_point() const;


    /**
    * @brief Get bearing vector of a point given an image coordinate
    * @return bearing vector
    */
    Vec3 operator () ( const Vec2& p ) const override;

    /**
    * @brief Transform a point from the camera plane to the image plane
    * @param p Camera plane point
    * @return Point on image plane
    */
    Vec2 cam2ima( const Vec2& p ) const override;

    /**
    * @brief Transform a point from the image plane to the camera plane
    * @param p Image plane point
    * @return camera plane point
    */
    Vec2 ima2cam(const Vec2& p) const override;

    /**
    * @brief Does the camera model handle a distortion field?
    * @retval false if intrinsic does not hold distortion
    */
    bool have_disto() const override;

    /**
    * @brief Add the distortion field to a point (that is in normalized camera frame)
    * @param p Point before distortion computation (in normalized camera frame)
    * @return point with distortion
    */
    Vec2 add_disto( const Vec2& p ) const override;

    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    */
    Vec2 remove_disto( const Vec2& p ) const override;

    /**
    * @brief Normalize a given unit pixel error to the camera plane
    * @param value Error in image plane
    * @return error of passing from the image plane to the camera plane
    */
    double imagePlane_toCameraPlaneError( double value ) const override;

    /**
    * @brief Return the projection matrix (interior & exterior) as a simplified projective projection
    * @param pose Extrinsic matrix
    * @return Concatenation of intrinsic matrix and extrinsic matrix
    */
    Mat34 get_projective_equivalent( const geometry::Pose3 & pose ) const override;


    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    std::vector<double> getParams() const override;


    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    bool updateFromParams(const std::vector<double> & params) override;

    /**
    * @brief Return the list of parameter indexes that must be held constant
    * @param parametrization The given parametrization
    */
    std::vector<int> subsetParameterization(
      const Intrinsic_Parameter_Type & parametrization) const override;

    /**
    * @brief Return the un-distorted pixel (with removed distortion)
    * @param p Input distorted pixel
    * @return Point without distortion
    */
    Vec2 get_ud_pixel( const Vec2& p ) const override;

    /**
    * @brief Return the distorted pixel (with added distortion)
    * @param p Input pixel
    * @return Distorted pixel
    */
    Vec2 get_d_pixel( const Vec2& p ) const override;

    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    void save( Archive & ar ) const;


    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar );

    /**
    * @brief Clone the object
    * @return A clone (copy of the stored object)
    */
    IntrinsicBase * clone(void) const override;
};

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_HPP
