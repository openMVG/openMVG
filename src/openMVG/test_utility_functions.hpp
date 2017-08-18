#pragma once
#include <random>
#include "third_party/ceres-solver/include/ceres/rotation.h"
#include "openMVG/sfm/sfm_data.hpp"

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
    poses.insert(std::pair<openMVG::IndexT, openMVG::geometry::Pose3>(i, pose));
  }

  return poses;
}

// we use an eigen matrix as return value in order to be able
// to use the EXPECT_MATRIX_NEAR macro during the test
inline openMVG::Mat computeDistances(const openMVG::sfm::SfM_Data & sfm_data)
{
  // put all positions of all poses and landmarks in a single vector
  std::vector<openMVG::Vec3> positions;

  for (const auto & pose : sfm_data.poses)
    positions.push_back(pose.second.center());
  for (const auto & landmark : sfm_data.GetLandmarks())
    positions.push_back(landmark.second.X);

  // compute distances between all poses and landmarks
  const int n_total = positions.size();

  openMVG::Mat distances(n_total, n_total);
  for (int i(0); i < n_total; i++)
  {
    for (int j(0); j < n_total; j++)
    {
      distances(i,j) = (positions[j] - positions[i]).norm();
    }
  }
  return distances;
}

inline openMVG::sfm::SfM_Data generate_random_poses_and_landmarks_in_scene(const int n_poses, const int n_landmarks)
{
  openMVG::sfm::SfM_Data sfm_data;

  for (int i(0); i < n_poses;i++)
  {
    // random orientation
    openMVG::Mat3 orientation = randomRotationMatrix();
    openMVG::Vec3 center = randomVector();
    openMVG::geometry::Pose3 pose(orientation, center);

    sfm_data.poses[i] = pose;
  }

  for (int i(0); i < n_landmarks; i++)
  {
    openMVG::sfm::Landmark landmark;
    openMVG::Vec3 position = randomVector();
    landmark.X = position;
    sfm_data.structure[i] = landmark;
  }

  return sfm_data;
}
