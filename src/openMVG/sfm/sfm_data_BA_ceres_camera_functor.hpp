// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP
#define OPENMVG_SFM_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"
#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include "openMVG/cameras/Camera_Spherical.hpp"

//--
//- Define ceres Cost_functor for each OpenMVG camera model
//--

namespace openMVG {
namespace sfm {


/// Decorator used to Weight a given cost camera functor
/// i.e useful to weight GCP (Ground Control Points)
template <typename CostFunctor>
struct WeightedCostFunction
{
  WeightedCostFunction() :weight_(1.0) {}

  explicit WeightedCostFunction
  (
    CostFunctor * func,
    const double weight
  )
    :functor_(func), weight_(weight)
  {}

  template <typename T>
  bool operator()
  (
    const T* const cam_intrinsic,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  ) const
  {
    if (functor_->operator()(cam_intrinsic, cam_extrinsics, pos_3dpoint, out_residuals))
    {
      // Reweight the residual values
      for (int i = 0; i < CostFunctor::num_residuals(); ++i)
      {
        out_residuals[i] *= T(weight_);
      }
      return true;
    }
    return false;
  }

  template <typename T>
  bool operator()
  (
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  ) const
  {
    if (functor_->operator()(cam_extrinsics, pos_3dpoint, out_residuals))
    {
      // Reweight the residual values
      for (int i = 0; i < CostFunctor::num_residuals(); ++i)
      {
        out_residuals[i] *= T(weight_);
      }
      return true;
    }
    return false;
  }

  ceres::internal::scoped_ptr<CostFunctor> functor_;
  const double weight_;
};

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
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

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

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic, 2, 3, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>, 2, 3, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic>
            (new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
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
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], K1 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

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

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r_coeff = 1.0 + k1*r2;
    const T x_d = x_u * r_coeff;
    const T y_d = y_u * r_coeff;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1, 2, 4, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1>, 2, 4, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
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
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

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

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T x_d = x_u * r_coeff;
    const T y_d = y_u * r_coeff;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3>, 2, 6, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
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
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
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
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, t1, t2 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

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

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& t1 = cam_intrinsics[OFFSET_DISTO_T1];
    const T& t2 = cam_intrinsics[OFFSET_DISTO_T2];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r4 = r2 * r2;
    const T r6 = r4 * r2;
    const T r_coeff = (1.0 + k1*r2 + k2*r4 + k3*r6);
    const T t_x = t2 * (r2 + 2.0 * x_u*x_u) + 2.0 * t1 * x_u * y_u;
    const T t_y = t1 * (r2 + 2.0 * y_u*y_u) + 2.0 * t2 * x_u * y_u;
    const T x_d = x_u * r_coeff + t_x;
    const T y_d = y_u * r_coeff + t_y;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2, 2, 8, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2>, 2, 8, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};


/**
 * @brief Ceres functor with constrained 3D points to use a Pinhole_Intrinsic_Fisheye
 *
 *  Data parameter blocks are the following <2,8,6,3>
 *  - 2 => dimension of the residuals,
 *  - 7 => the intrinsic data block [focal, principal point x, principal point y, K1, K2, K3, K4],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */

struct ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum : uint8_t {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2,
    OFFSET_DISTO_K1 = 3,
    OFFSET_DISTO_K2 = 4,
    OFFSET_DISTO_K3 = 5,
    OFFSET_DISTO_K4 = 6,
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y], k1, k2, k3, k4 )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

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
    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1 = cam_intrinsics[OFFSET_DISTO_K1];
    const T& k2 = cam_intrinsics[OFFSET_DISTO_K2];
    const T& k3 = cam_intrinsics[OFFSET_DISTO_K3];
    const T& k4 = cam_intrinsics[OFFSET_DISTO_K4];

    // Apply distortion (xd,yd) = disto(x_u,y_u)
    const T r2 = x_u*x_u + y_u*y_u;
    const T r = sqrt(r2);
    const T
      theta = atan(r),
      theta2 = theta*theta,
      theta3 = theta2*theta,
      theta4 = theta2*theta2,
      theta5 = theta4*theta,
      theta7 = theta3*theta3*theta, //thetha6*theta
      theta8 = theta4*theta4,
      theta9 = theta8*theta;
    const T theta_dist = theta + k1*theta3 + k2*theta5 + k3*theta7 + k4*theta9;
    const T inv_r = r > T(1e-8) ? T(1.0)/r : T(1.0);
    const T cdist = r > T(1e-8) ? theta_dist * inv_r : T(1.0);

    const T x_d = x_u * cdist;
    const T y_d = y_u * cdist;

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_d;
    const T projected_y = principal_point_y + focal * y_d;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye, 2, 7, 6, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye>, 2, 7, 6, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct ResidualErrorFunctor_Intrinsic_Spherical
{
  ResidualErrorFunctor_Intrinsic_Spherical
  (
    const double* const pos_2dpoint,
    const size_t imageSize_w,
    const size_t imageSize_h
  )
  : m_pos_2dpoint(pos_2dpoint),
    m_imageSize{imageSize_w, imageSize_h}
  {
  }

  /**
  * @param[in] cam_extrinsics: Camera parameterized using one block of 6 parameters [R;t]:
  *   - 3 for rotation(angle axis), 3 for translation
  * @param[in] pos_3dpoint
  * @param[out] out_residuals
  */
  template <typename T>
  bool operator()
  (
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  )
  const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the coord in is Image space
    const T lon = ceres::atan2(pos_proj[0] , pos_proj[2]); // Horizontal normalization of the  X-Z component
    const T lat = ceres::atan2(-pos_proj[1], ceres::sqrt(pos_proj[0] * pos_proj[0]  + pos_proj[2] * pos_proj[2])); // Tilt angle
    const T coord[] = {lon / (2 * M_PI), lat / (2 * M_PI)}; // normalization

    const T size ( std::max(m_imageSize[0], m_imageSize[1]) );
    const T projected_x = coord[0] * size - 0.5 + m_imageSize[0] / 2.0;
    const T projected_y = coord[1] * size - 0.5 + m_imageSize[1] / 2.0;

    out_residuals[0] = projected_x - m_pos_2dpoint[0];
    out_residuals[1] = projected_y - m_pos_2dpoint[1];

    return true;
  }

  static const int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const cameras::IntrinsicBase * cameraInterface,
    const Vec2 & observation,
    const double weight = 0.0
  )
  {

    if (weight == 0.0)
    {
      return
          (new ceres::AutoDiffCostFunction
              <ResidualErrorFunctor_Intrinsic_Spherical, 2, 6, 3>(
              new ResidualErrorFunctor_Intrinsic_Spherical(
                observation.data(), cameraInterface->w(), cameraInterface->h())));
    }
    else
    {
      return
          (new ceres::AutoDiffCostFunction
              <WeightedCostFunction<ResidualErrorFunctor_Intrinsic_Spherical>, 2, 6, 3>
              (new WeightedCostFunction<ResidualErrorFunctor_Intrinsic_Spherical>
                   (new ResidualErrorFunctor_Intrinsic_Spherical(
                     observation.data(), cameraInterface->w(), cameraInterface->h()),
                     weight)));
    }
  }

  const double * m_pos_2dpoint;  // The 2D observation
  size_t         m_imageSize[2]; // The image width and height
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_BA_CERES_CAMERA_FUNCTOR_HPP
