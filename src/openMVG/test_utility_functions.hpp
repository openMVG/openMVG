#pragma once
#include <random>
#include "third_party/ceres-solver/include/ceres/rotation.h"
#include "openMVG/sfm/sfm_data.hpp"

inline int randomPositiveInteger(int boundary = 100)
{
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_int_distribution<> distr(0, boundary);
  return distr(eng);
}

inline double randomPositiveDouble(double boundary = 100.0)
{
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_real_distribution<double> distr(1e-5, boundary);
  return distr(eng);
}

inline openMVG::Vec3 randomVector(double boundary = 100.0)
{
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_real_distribution<double> distr(-boundary, boundary);

  return openMVG::Vec3(distr(eng), distr(eng), distr(eng));
}

inline openMVG::Mat3 randomRotationMatrix()
{
  // generate random unitary vector
  openMVG::Vec3 rotation_axis = randomVector(1.0);
  rotation_axis /= rotation_axis.norm();

  // generate random angle
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_real_distribution<double> distr(-2*M_PI, 2*M_PI);

  // create rotation matrix from angle axis representation
  rotation_axis *= distr(eng);
  openMVG::Mat3 rotation_matrix;
  ceres::AngleAxisToRotationMatrix(rotation_axis.data(), rotation_matrix.data());

  return rotation_matrix;
}

inline std::vector<double> randomQuaternion()
{
  // generate random unitary vector
  openMVG::Vec3 rotation_axis = randomVector();
  rotation_axis /= rotation_axis.norm();

  // generate random angle
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_real_distribution<double> distr(-2*M_PI, 2*M_PI);
  rotation_axis *= distr(eng);

  // create quaternion from angle axis
  std::vector<double> quaternion(4);
  ceres::AngleAxisToQuaternion(rotation_axis.data(), &quaternion[0]);

  return quaternion;
}

inline openMVG::sfm::Poses generateRandomPoses(const int n_poses)
{
  openMVG::sfm::Poses poses;
  for (int i(0); i < n_poses; i++)
  {
    const openMVG::geometry::Pose3 pose = openMVG::geometry::Pose3(randomRotationMatrix(), randomVector(20.0));
    poses.insert({i, pose});
  }

  return poses;
}

// we use an eigen matrix as return value in order to be able
// to use the EXPECT_MATRIX_NEAR macro during the test
inline openMVG::Mat computeDistancesBetweenPosesAndLandmarks(const openMVG::sfm::SfM_Data & sfm_data)
{
  enum EPointType
  {
    pLandmark, pPose
  };

  // put all positions of all poses and landmarks in a single vector
  std::vector<openMVG::Vec3> positions;
  // store indices and point type in a second vector
  std::vector<std::pair<openMVG::IndexT, EPointType>> indices;

  int max_pose_id(-1);
  for (const auto & pose : sfm_data.poses)
  {
    const openMVG::IndexT & pose_id = pose.first;
    positions.push_back(pose.second.center());
    indices.emplace_back(pose_id, pPose);
    max_pose_id = std::max((int)pose_id, max_pose_id);
  }

  int max_landmark_id(-1);
  for (const auto & landmark : sfm_data.structure)
  {
    const openMVG::IndexT & landmark_id = landmark.first;
    positions.push_back(landmark.second.X);
    indices.emplace_back(landmark.first, pLandmark);
    max_landmark_id = std::max((int)landmark_id, max_landmark_id);
  }

  // compute distances between all poses and landmarks
  const int n_total = positions.size();
  openMVG::Mat distances(max_landmark_id + max_pose_id + 2, max_landmark_id + max_pose_id + 2);
  for (int i(0); i < n_total; i++)
  {
    for (int j(0); j < n_total; j++)
    {
      const std::pair<openMVG::IndexT, EPointType> corresponding_index_pair_i = indices.at(i);
      const std::pair<openMVG::IndexT, EPointType> corresponding_index_pair_j = indices.at(j);
      openMVG::IndexT corresponding_i = corresponding_index_pair_i.second == pPose ? corresponding_index_pair_i.first : corresponding_index_pair_i.first + max_pose_id + 1;
      openMVG::IndexT corresponding_j = corresponding_index_pair_j.second == pPose ? corresponding_index_pair_j.first : corresponding_index_pair_j.first + max_pose_id + 1;
      distances(corresponding_i,corresponding_j) = (positions[i] - positions[j]).norm();
    }
  }
  return distances;
}

inline openMVG::sfm::SfM_Data generate_random_poses_and_landmarks_in_scene(const int n_poses, const int n_landmarks)
{
  openMVG::sfm::SfM_Data sfm_data;

  for (openMVG::IndexT i(0); i < n_poses;i++)
  {
    // random orientation
    openMVG::Mat3 orientation = randomRotationMatrix();
    openMVG::Vec3 center = randomVector();
    openMVG::geometry::Pose3 pose(orientation, center);

    sfm_data.poses[i] = pose;
  }

  for (openMVG::IndexT i(0); i < n_landmarks; i++)
  {
    openMVG::sfm::Landmark landmark;
    openMVG::Vec3 position = randomVector();
    landmark.X = position;
    sfm_data.structure[i] = landmark;
  }

  return sfm_data;
}
