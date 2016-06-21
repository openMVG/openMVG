// Copyright (c) 2016 Pierre MOULON, Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "opengv/types.hpp"
#include "opengv/relative_pose/methods.hpp"
#include "opengv/relative_pose/NoncentralRelativeAdapter.hpp"
#include "opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp"

#include "openMVG/geometry/pose3.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <fstream>

using namespace openMVG;
using namespace opengv;


struct Rig_Dataset2
{
  std::vector<geometry::Pose3> sub_poses; // The rigidly coupled camera (local coordinates system)

  // Poses on which the RIG will be placed on
  geometry::Pose3 pose_1;
  geometry::Pose3 pose_2;

  // 3D observed structure
  Mat3X X;

  void generate2D2DCorrespondences
  (
    bearingVectors_t &bearingVectors1,
    bearingVectors_t &bearingVectors2,
    std::vector<int> &camCorrespondences1,
    std::vector<int> &camCorrespondences2
  )
  {
    //create the 2D3D-correspondences by looping through the cameras
    size_t numberCams = sub_poses.size();
    size_t camCorrespondence = 0;

    for( size_t i = 0; i < (size_t) X.cols(); ++i)
    {
      //get the global camera pose (sub_pose * pose)
      const geometry::Pose3 current_pose1 (sub_poses[camCorrespondence] * pose_1);
      const geometry::Pose3 current_pose2 (sub_poses[camCorrespondence] * pose_2);

      // Compute the normalized bearing vector in the current view pose
      bearingVectors1.emplace_back( current_pose1( X.col(i)).normalized() );
      bearingVectors2.emplace_back( current_pose2( X.col(i)).normalized() );

      //push back the camera correspondences
      camCorrespondences1.push_back(camCorrespondence);
      camCorrespondences2.push_back(camCorrespondence++);
      if (camCorrespondence > (numberCams - 1) )
        camCorrespondence = 0;
    }
  }

  Mat3 generateRandomRotation(const double maxAngle = 2 * M_PI)
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, maxAngle);
    return RotationAroundX(dis(rd)) * RotationAroundY(dis(rd)) * RotationAroundY (dis(rd));
  }

  // Generate N random camera at a given distance to viewpoint origin
  void generateRandomCameraPose
  (
    const unsigned int nb_camera_in_a_rig,
    const unsigned int nb_points
  )
  {
    // Generate two camera pose
    pose_1 = geometry::Pose3(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0, 0, 2.0));
    pose_2 = geometry::Pose3(generateRandomRotation(0.5).transpose(), 2 * Vec3::Random());

    // generate a static rig (that will be put over the two previously defined pose)
    sub_poses.resize(nb_camera_in_a_rig);
    for (unsigned int i = 0; i < nb_camera_in_a_rig; ++i)
    {
      const double offset = 0.5; //this is the distance from the viewpoint origin
      sub_poses[i] = geometry::Pose3(generateRandomRotation(2*M_PI).transpose(), offset * Vec3::Random());
    }

    // Generate random points
    X.resize(3, nb_points);

    //initialize point-cloud
    const double minimumDepth = 4, maximumDepth = 8;

    X.setRandom();
    for( size_t i = 0; i < (size_t) X.cols(); i++ )
    {
      X.col(i) = (maximumDepth-minimumDepth) * X.col(i) + minimumDepth * X.col(i).normalized();
    }

    {
      std::ofstream outfile;
      outfile.open("final.ply", std::ios_base::out);
      if (outfile.is_open()) {
        outfile << "ply"
          << std::endl << "format ascii 1.0"
          << std::endl << "element vertex " << X.cols() + nb_camera_in_a_rig * 2
          << std::endl << "property float x"
          << std::endl << "property float y"
          << std::endl << "property float z"
          << std::endl << "property uchar red"
          << std::endl << "property uchar green"
          << std::endl << "property uchar blue"
          << std::endl << "end_header" << std::endl;

        //-- Export 3D point cloud
        for (int i = 0; i < X.cols(); ++i)
        {
          outfile << X.col(i).transpose() << " 255 255 255" << std::endl;
        }

        for (int i = 0; i < sub_poses.size(); ++i)
        {
          outfile << (sub_poses[i] * pose_2).center().transpose() << " 0 255 0" << std::endl;
        }

        for (int i = 0; i < sub_poses.size(); ++i)
        {
          outfile << (sub_poses[i] * pose_1).center().transpose() << " 255 0 0" << std::endl;
        }
      }
    }
  }
};

void
extractRelativePose(
    const geometry::Pose3 & pose1,
    const geometry::Pose3 & pose2,
    geometry::Pose3 & relativePose
)
{
  relativePose = geometry::Pose3(
    pose1.rotation().transpose() * pose2.rotation(),
    pose1.rotation().transpose() * (pose2.center() - pose1.center())
    );
}

TEST(SixPtRelativePoseTest, SixPtRelativePoseTest_Kernel ) {
{
    Rig_Dataset2 dataset;
    dataset.generateRandomCameraPose(6, 100);

    bearingVectors_t bearingVectors1;
    bearingVectors_t bearingVectors2;
    std::vector<int> camCorrespondences1;
    std::vector<int> camCorrespondences2;
    dataset.generate2D2DCorrespondences
    (
      bearingVectors1,  bearingVectors2,
      camCorrespondences1, camCorrespondences2
    );

    translations_t camOffsets(dataset.sub_poses.size());
    rotations_t camRotations(dataset.sub_poses.size());
    for (unsigned int i = 0; i < dataset.sub_poses.size(); ++i)
    {
      camOffsets[i] = dataset.sub_poses[i].center();
      camRotations[i] = dataset.sub_poses[i].rotation().transpose();

      std::cout <<  camOffsets[i] << "\n\n" << camRotations[i] << std::endl;
    }

    geometry::Pose3 relative_pose;
    extractRelativePose(dataset.pose_1, dataset.pose_2, relative_pose);
    std::cout << "GT\n"
     << relative_pose.rotation() << "\n"
     << relative_pose.center().transpose() << "\n";

    //create non-central relative adapter
    relative_pose::NoncentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2,
      camCorrespondences1,
      camCorrespondences2,
      camOffsets,
      camRotations
      //, relative_pose.center(), relative_pose.rotation().transpose()
      );

    {
      // create non central relative sac problem
      sac_problems::relative_pose::NoncentralRelativePoseSacProblem
                problem(adapter,
                        //sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SIXPT,
                        sac_problems::relative_pose::NoncentralRelativePoseSacProblem::GE,
                        //sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SEVENTEENPT,
                        false);

      // solve pose problem
      transformation_t relativePose;
      std::vector<int> idx(dataset.X.cols());
      std::iota(idx.begin(), idx.end(), 0);
      problem.computeModelCoefficients(idx, relativePose);

      std::cout << relativePose << "\n" << std::endl;

      std::cout << "Rotation error: " << FrobeniusDistance((Mat3)relative_pose.rotation().transpose(), Mat3(relativePose.block<3,3>(0,0))) << std::endl;
      std::cout << "Translation error: " << (relativePose.col(3)/relativePose.col(3).norm() - relative_pose.center()/relative_pose.center().norm()).norm() << std::endl;

      adapter.sett12(relativePose.col(3));
      adapter.setR12(relativePose.block<3,3>(0,0));
      relativePose = relative_pose::optimize_nonlinear(adapter);

      std::cout << relativePose << "\n" << std::endl;

      std::cout << "Rotation error: " << FrobeniusDistance((Mat3)relative_pose.rotation().transpose(), Mat3(relativePose.block<3,3>(0,0))) << std::endl;
      std::cout << "Translation error: " << (relativePose.col(3)/relativePose.col(3).norm() - relative_pose.center()/relative_pose.center().norm()).norm() << std::endl;
    }
  }
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
