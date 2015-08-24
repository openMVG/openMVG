// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP
#define OPENMVG_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP

#include "openMVG/cameras/cameras.hpp"
#include "ceres/rotation.h"

//--
//- Define ceres Cost_functor for each OpenMVG camera model
//--

namespace openMVG {
namespace sfm {

/**
 * @brief Ceres functor to use a Pinhole_Intrinsic (pinhole camera model K[R[t]) and a 3D point.
 *
 *  Data parameter blocks are the following <2,3,6,3>
 *  - 2 => dimension of the residuals,
 *  - 3 => the intrinsic data block [focal, principal point x, principal point y],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic
{
  ResidualErrorFunctor_Pinhole_Intrinsic(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_K: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_Rt: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_K,
    const T* const cam_Rt,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_Rt;
    const T * cam_t = &cam_Rt[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_K[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_K[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_K[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  double m_pos_2dpoint[2]; // The 2D observation
};

/**
 * @brief Ceres functor to use a Pinhole_Intrinsic_Radial_K1
 *
 *  Data parameter blocks are the following <2,4,6,3>
 *  - 2 => dimension of the residuals,
 *  - 4 => the intrinsic data block [focal, principal point x, principal point y, K1],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3
  };

  /**
   * @param[in] cam_K: Camera intrinsics( focal, principal point [x,y], K1 )
   * @param[in] cam_Rt: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_K,
    const T* const cam_Rt,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_Rt;
    const T * cam_t = &cam_Rt[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_K[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_K[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_K[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_K[OFFSET_DISTO_K1];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r_coeff = (T(1) + k1*r2);
    const T x_d = x_u * r_coeff;
    const T y_d = y_u * r_coeff;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  double m_pos_2dpoint[2]; // The 2D observation
};

/**
 * @brief Ceres functor to use a Pinhole_Intrinsic_Radial_K3
 *
 *  Data parameter blocks are the following <2,6,6,3>
 *  - 2 => dimension of the residuals,
 *  - 6 => the intrinsic data block [focal, principal point x, principal point y, K1, K2, K3],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
  };

  /**
   * @param[in] cam_K: Camera intrinsics( focal, principal point [x,y], k1, k2, k3 )
   * @param[in] cam_Rt: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_K,
    const T* const cam_Rt,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_Rt;
    const T * cam_t = &cam_Rt[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_K[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_K[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_K[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_K[OFFSET_DISTO_K1];
    const T& k2 = cam_K[OFFSET_DISTO_K2];
    const T& k3 = cam_K[OFFSET_DISTO_K3];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (T(1) + k1*r2 + k2*r4 + k3*r6);
    const T x_d = x_u * r_coeff;
    const T y_d = y_u * r_coeff;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  double m_pos_2dpoint[2]; // The 2D observation
};

/**
 * @brief Ceres functor with constrained 3D points to use a Pinhole_Intrinsic_Brown_T2
 *
 *  Data parameter blocks are the following <2,8,6,3>
 *  - 2 => dimension of the residuals,
 *  - 8 => the intrinsic data block [focal, principal point x, principal point y, K1, K2, K3, T1, T2],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
    OFFSET_DISTO_T1 = 6,
    OFFSET_DISTO_T2 = 7,
  };

  /**
   * @param[in] cam_K: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_Rt: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_K,
    const T* const cam_Rt,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_Rt;
    const T * cam_t = &cam_Rt[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_K[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_K[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_K[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_K[OFFSET_DISTO_K1];
    const T& k2 = cam_K[OFFSET_DISTO_K2];
    const T& k3 = cam_K[OFFSET_DISTO_K3];
    const T& t1 = cam_K[OFFSET_DISTO_T1];
    const T& t2 = cam_K[OFFSET_DISTO_T2];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (T(1) + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + T(2) * x_u*x_u) + T(2) * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + T(2) * y_u*y_u) + T(2) * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  double m_pos_2dpoint[2]; // The 2D observation
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP
