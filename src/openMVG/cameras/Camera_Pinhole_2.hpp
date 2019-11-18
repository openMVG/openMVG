#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_2_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_2_HPP

#include <vector>

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{
namespace cameras
{

/**
* @brief Define an ideal Pinhole camera intrinsics (store a K 3x3 matrix)
* with intrinsic parameters defining the K calibration matrix
*
* Intrinsic camera matrix is \f$ K = \begin{pmatrix} fx & 0 & u_0 \\ 0 & fy & v_0 \\ 0 & 0 & 1 \end{pmatrix} \f$
*
* @note This is an ideal Pinhole camera because it doesn't handle skew and distortion
* @note The camera does only handle two focal length (ie: \f$ f_x and f_y can be different \f$ )
*/
class Pinhole_Intrinsic_2 : public IntrinsicBase
{
  using class_type = Pinhole_Intrinsic_2;

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
    * @param focal_length_x_pix Focal length x (in pixel) of the camera
    * @param focal_length_y_pix Focal length y (in pixel) of the camera
    * @param ppx Principal point on x-axis
    * @param ppy Principal point on y-axis
    */
    Pinhole_Intrinsic_2(
      unsigned int w = 0, unsigned int h = 0,
      double focal_length_x_pix = 0.0, double focal_length_y_pix = 0.0,
      double ppx = 0.0, double ppy = 0.0 )
      : IntrinsicBase( w, h )
    {
      K_ << focal_length_x_pix, 0., ppx, 0., focal_length_y_pix, ppy, 0., 0., 1.;
      Kinv_ = K_.inverse();
    }

    /**
    * @brief Constructor
    * @param w Width of the image plane
    * @param h Height of the image plane
    * @param K Intrinsic Matrix (3x3) {f_x,0,ppx; 0,f_y,ppy; 0,0,1}
    */
    Pinhole_Intrinsic_2(
      unsigned int w,
      unsigned int h,
      const Mat3& K)
      : IntrinsicBase( w, h ), K_(K)
    {
      K_(0,0) = K(0,0);
      K_(1,1) = K(1,1);
      Kinv_ = K_.inverse();
    }

    /**
    * @brief Destructor
    */
    ~Pinhole_Intrinsic_2() override = default;

    /**
    * @brief Get type of the intrinsic
    * @retval PINHOLE_CAMERA_2
    */
    EINTRINSIC getType() const override
    {
      return PINHOLE_CAMERA_2;
    }

    /**
    * @brief Get the intrinsic matrix
    * @return 3x3 intrinsic matrix
    */
    const Mat3& K() const
    {
      return K_;
    }

    /**
    * @brief Get the inverse of the intrinsic matrix
    * @return Inverse of intrinsic matrix
    */
    const Mat3& Kinv() const
    {
      return Kinv_;
    }


    /**
    * @brief Return the value of the focalx in pixels
    * @return Focal x of the camera (in pixel)
    */
    inline double focalx() const
    {
      return K_( 0, 0 );
    }

    /**
    * @brief Return the value of the focaly in pixels
    * @return Focal y of the camera (in pixel)
    */
    inline double focaly() const
    {
      return K_( 1, 1 );
    }

    /**
    * @brief Get principal point of the camera
    * @return Principal point of the camera
    */
    inline Vec2 principal_point() const
    {
      return {K_( 0, 2 ), K_( 1, 2 )};
    }

    /**
    * @brief Get bearing vectors from image coordinates
    * @return bearing vectors
    */
    Mat3X operator () ( const Mat2X& points ) const override
    {
      return Kinv_ * points.colwise().homogeneous();
    }

    /**
    * @brief Transform a point from the camera plane to the image plane
    * @param p Camera plane point
    * @return Point on image plane
    */
    Vec2 cam2ima( const Vec2& p ) const override
    {
      double p0 = focalx() * p(0);
      double p1 = focaly() * p(1);
      return Vec2(p0, p1) + principal_point();
    }

    /**
    * @brief Transform a point from the image plane to the camera plane
    * @param p Image plane point
    * @return camera plane point
    */
    Vec2 ima2cam( const Vec2& p ) const override
    {
      Vec2 p_im = p -  principal_point();
      return Vec2( p_im(0) / focalx(), p_im(1) / focaly() );
    }

    /**
    * @brief Does the camera model handle a distortion field?
    * @retval false if intrinsic does not hold distortion
    */
    bool have_disto() const override
    {
      return false;
    }

    /**
    * @brief Add the distortion field to a point (that is in normalized camera frame)
    * @param p Point before distortion computation (in normalized camera frame)
    * @return point with distortion
    */
    Vec2 add_disto( const Vec2& p ) const override
    {
      return p;
    }

    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    */
    Vec2 remove_disto( const Vec2& p ) const override
    {
      return p;
    }

    /**
    * @brief Normalize a given unit pixel error to the camera plane
    * @param value Error in image plane
    * @return error of passing from the image plane to the camera plane
    */
    double imagePlane_toCameraPlaneError( double value ) const override
    {
      return value * 2.0 / (focalx() + focaly());
    }

    /**
    * @brief Return the projection matrix (interior & exterior) as a simplified projective projection
    * @param pose Extrinsic matrix
    * @return Concatenation of intrinsic matrix and extrinsic matrix
    */
    Mat34 get_projective_equivalent( const geometry::Pose3 & pose ) const override
    {
      Mat34 P;
      P_From_KRt( K(), pose.rotation(), pose.translation(), &P );
      return P;
    }


    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    std::vector<double> getParams() const override
    {
      return  { K_(0, 0), K_(1, 1), K_(0, 2), K_(1, 2) };
    }


    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    bool updateFromParams(const std::vector<double> & params) override
    {
      if ( params.size() == 3 )
      {
        *this = Pinhole_Intrinsic_2( w_, h_, params[0], params[1], params[2], params[3] );
        return true;
      }
      else
      {
        return false;
      }
    }

    /**
    * @brief Return the list of parameter indexes that must be held constant
    * @param parametrization The given parametrization
    */
    std::vector<int> subsetParameterization
    (
      const Intrinsic_Parameter_Type & parametrization) const override
    {
      std::vector<int> constant_index;
      const int param = static_cast<int>(parametrization);
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH)
           || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(0);
        constant_index.push_back(1);
      }
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT)
          || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(2);
        constant_index.push_back(3);
      }
      return constant_index;
    }

    /**
    * @brief Return the un-distorted pixel (with removed distortion)
    * @param p Input distorted pixel
    * @return Point without distortion
    */
    Vec2 get_ud_pixel( const Vec2& p ) const override
    {
      return p;
    }

    /**
    * @brief Return the distorted pixel (with added distortion)
    * @param p Input pixel
    * @return Distorted pixel
    */
    Vec2 get_d_pixel( const Vec2& p ) const override
    {
      return p;
    }

    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    inline void save( Archive & ar ) const;


    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    inline void load( Archive & ar );

    /**
    * @brief Clone the object
    * @return A clone (copy of the stored object)
    */
    IntrinsicBase * clone( void ) const override
    {
      return new class_type( *this );
    }
};

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_2_HPP
