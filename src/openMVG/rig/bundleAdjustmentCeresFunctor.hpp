#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <openMVG/geometry/pose3.hpp>

class ResidualErrorMainCameraFunctor
{
public :
  ResidualErrorMainCameraFunctor(const openMVG::Mat3 & K, const openMVG::Vec2 & x, const openMVG::Vec3 & X) // const double* const pos_2dpoint
  {
    // Set the intrinsics
     _K = K;
    
    // Set the observation
    _observation[0] = x[0];
    _observation[1] = x[1];
    
    // Set the 3D point
    _point(0) = X(0);
    _point(1) = X(1);
    _point(2) = X(2);
    
    // 
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  //enum {
  //  OFFSET_FOCAL_X_LENGTH = 0,
  //  OFFSET_FOCAL_Y_LENGTH = 1,
  //  OFFSET_PRINCIPAL_POINT_X = 2,
  //  OFFSET_PRINCIPAL_POINT_Y = 3
  //};

  /**
   * @param[in] cam_K: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_Rt: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_Rt,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--
    
    const T * cam_R = cam_Rt;
    const T * cam_t = &cam_Rt[3];
    
    T pos_3dpoint[3];
    pos_3dpoint[0]= T(_point(0)); // better solution ? todo@L
    pos_3dpoint[1]= T(_point(1));
    pos_3dpoint[2]= T(_point(2));

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

    const T& focal_x = T(_K(0,0));//.data()[OFFSET_FOCAL_X_LENGTH];
    const T& focal_y = T(_K(1,1));//[OFFSET_FOCAL_Y_LENGTH];
    const T& principal_point_x = T(_K(0,2));//[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = T(_K(1,2));//[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal_x * x_u;
    const T projected_y = principal_point_y + focal_y * y_u;

    const T img_x = principal_point_x + focal_x * T(_observation[0]);
    const T img_y = principal_point_y + focal_y * T(_observation[1]);  
    
    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - img_x;
    out_residuals[1] = projected_y - img_y;
 
#if 0
    std::size_t semiWidth = 3.0;
    cv::rectangle(imgRes, 
           cvPoint(projected_x-T(semiWidth),projected_y-T(semiWidth)),
           cvPoint(projected_x+T(semiWidth),projected_y+T(semiWidth)),
           cv::Scalar(255,0,0));

   cv::rectangle(imgRes,
           cvPoint(_observation[0]-semiWidth,_observation[1]-semiWidth),
           cvPoint(_observation[0]+semiWidth,_observation[1]+semiWidth),
           cv::Scalar(0,255,0));
     cv::imshow("Reprojection", imgRes);
     popart::vision::cvpause();
#endif
    
    return true;
  }

private :
  //double m_pos_2dpoint[2]; // The 2D observation
  openMVG::Mat3 _K;
  openMVG::Vec2 _observation;
  openMVG::Vec3 _point;
};

class ResidualErrorSecondaryCameraFunctor
{
public :
  ResidualErrorSecondaryCameraFunctor(const openMVG::Mat3 & K, const openMVG::Vec2 & x, const openMVG::Vec3 & X) // const double* const pos_2dpoint
  {
    // Set the intrinsics
     _K = K;
    
    // Set the observation
    _observation[0] = x[0];
    _observation[1] = x[1];
    
    // Set the 3D point
    _point(0) = X(0);
    _point(1) = X(1);
    _point(2) = X(2);
    
    // 
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  //enum {
  //  OFFSET_FOCAL_X_LENGTH = 0,
  //  OFFSET_FOCAL_Y_LENGTH = 1,
  //  OFFSET_PRINCIPAL_POINT_X = 2,
  //  OFFSET_PRINCIPAL_POINT_Y = 3
  //};

  /**
   * @param[in] cam_K: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_Rt: Camera parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_Rt_main,
    const T* const cam_Rt_relative,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * RMain = cam_Rt_main;
    const T * tMain = &cam_Rt_main[3];
    
    const T * RRelative = cam_Rt_relative;
    const T * tRelative = &cam_Rt_relative[3];  
    
    //getAbsolutePose(cam_Rt_main, cam_Rt_relative, cam_R, cam_t);
    
    T pos_3dpoint[3];
    pos_3dpoint[0]= T(_point(0)); // better solution ? todo@L
    pos_3dpoint[1]= T(_point(1));
    pos_3dpoint[2]= T(_point(2));

    T pos_tmp[3];
    // Rotate the point according the relative rotation first
    ceres::AngleAxisRotatePoint(RMain, pos_3dpoint, pos_tmp);

    // Apply the relative translation first
    pos_tmp[0] += tMain[0];
    pos_tmp[1] += tMain[1];
    pos_tmp[2] += tMain[2];
    
    T pos_proj[3];
    // Rotate the point according the main camera rotation
    ceres::AngleAxisRotatePoint(RRelative, pos_tmp, pos_proj);
    
    // Apply the main camera translation
    pos_proj[0] += tRelative[0];
    pos_proj[1] += tRelative[1];
    pos_proj[2] += tRelative[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal_x = T(_K(0,0));//.data()[OFFSET_FOCAL_X_LENGTH];
    const T& focal_y = T(_K(1,1));//[OFFSET_FOCAL_Y_LENGTH];
    const T& principal_point_x = T(_K(0,2));//[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = T(_K(1,2));//[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal_x * x_u;
    const T projected_y = principal_point_y + focal_y * y_u;

    const T img_x = principal_point_x + focal_x * T(_observation[0]);
    const T img_y = principal_point_y + focal_y * T(_observation[1]);  
    
    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - img_x;
    out_residuals[1] = projected_y - img_y;
    
    
    //POPART_COUT("projected : " << projected_x << "," << projected_y);
    //POPART_COUT("img : " << img_x << "," << img_y);
    
    return true;
  }

private :
  //double m_pos_2dpoint[2]; // The 2D observation
  openMVG::Mat3 _K;
  openMVG::Vec2 _observation;
  openMVG::Vec3 _point;
};