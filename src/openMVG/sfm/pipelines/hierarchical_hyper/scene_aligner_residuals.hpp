#pragma once

#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace openMVG {
namespace sfm {

void printResiduals(const ceres::Problem & problem, const std::vector<ceres::ResidualBlockId> & residuals_vector);


// !! \\ WHEN CONSTRUCTING : by default the first submap is the one with the fixed base node. We will only modify
//  the position of the second submap's basenode ! TODO : make this less error-prone somehow
struct ResidualErrorFunctor_BaseNode_Separators
{
  ResidualErrorFunctor_BaseNode_Separators(const double * const pos_3dmeasurement_fixed_submap, const double * const pos_3dmeasurement_moving_submap)
    :m_pos_3dmeasurement_A(pos_3dmeasurement_fixed_submap), m_pos_3dmeasurement_B(pos_3dmeasurement_moving_submap)
  {}

  /**
   * @param[in] pose_base_node : base node parameterized using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[in] scaling factor
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const pose_base_node, // the coordinates of the second base node
    const T* const pos_3dpoint, // the position of the landmarks (remember that the first base node is the origin)
    const T* const scaling_factor,
    T* out_residuals) const
  {
    // here pos_3dpoint is the position of the 3d point as described in the
    // first base node referential. pose_base_node is the posiiton and orientation
    // of the second base node referential.
    // moreover, m_pos_3dmeasurement contains the measurement of the 3d point in the second
    // base node referential. we want to minimize the error between this measurement and the
    // projected position of the 3d point in the second base node referential.
    
    // find projection of pos_3dpoint in second referential
    const T * base_node_R = pose_base_node;
    const T * base_node_t = &pose_base_node[3];
    
    T new_position[3];
    T translated_only[3];

    // translate to new origin
    translated_only[0] = pos_3dpoint[0] - base_node_t[0];
    translated_only[1] = pos_3dpoint[1] - base_node_t[1];
    translated_only[2] = pos_3dpoint[2] - base_node_t[2];

    // in angle-axis reprasentation, the inverse rotation is represented by negating the vector
    T inverse_rotation[3];
    inverse_rotation[0] = -base_node_R[0];
    inverse_rotation[1] = -base_node_R[1];
    inverse_rotation[2] = -base_node_R[2];

    // rotate around new origin
    ceres::AngleAxisRotatePoint(inverse_rotation, translated_only, new_position);

    // scale distance to origin
    new_position[0] *= *scaling_factor;
    new_position[1] *= *scaling_factor;
    new_position[2] *= *scaling_factor;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    // residuals corresponding to the second submap
    out_residuals[0] = new_position[0] - T(m_pos_3dmeasurement_B[0]);
    out_residuals[1] = new_position[1] - T(m_pos_3dmeasurement_B[1]);
    out_residuals[2] = new_position[2] - T(m_pos_3dmeasurement_B[2]);

    // residuals corresponding to the first submap
    out_residuals[3] = pos_3dpoint[0] - T(m_pos_3dmeasurement_A[0]);
    out_residuals[4] = pos_3dpoint[1] - T(m_pos_3dmeasurement_A[1]);
    out_residuals[5] = pos_3dpoint[2] - T(m_pos_3dmeasurement_A[2]);

    return true;
  }

  static int num_residuals() { return 6; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec3 & measurement_A,
    const Vec3 & measurement_B,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_BaseNode_Separators, 6, 6, 3, 1>(
            new ResidualErrorFunctor_BaseNode_Separators(measurement_A.data(), measurement_B.data())));
    }
  }

  const double * m_pos_3dmeasurement_A; // the 3d measurement in the first submap
  const double * m_pos_3dmeasurement_B; // the 3d measurement in the second submap
};

}
}
