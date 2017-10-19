// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: joydeepb@ri.cmu.edu (Joydeep Biswas)
//
// This example demonstrates how to use the DynamicAutoDiffCostFunction
// variant of CostFunction. The DynamicAutoDiffCostFunction is meant to
// be used in cases where the number of parameter blocks or the sizes are not
// known at compile time.
//
// This example simulates a robot traversing down a 1-dimension hallway with
// noise odometry readings and noisy range readings of the end of the hallway.
// By fusing the noisy odometry and sensor readings this example demonstrates
// how to compute the maximum likelihood estimate (MLE) of the robot's pose at
// each timestep.
//
// The robot starts at the origin, and it is travels to the end of a corridor of
// fixed length specified by the "--corridor_length" flag. It executes a series
// of motion commands to move forward a fixed length, specified by the
// "--pose_separation" flag, at which pose it receives relative odometry
// measurements as well as a range reading of the distance to the end of the
// hallway. The odometry readings are drawn with Gaussian noise and standard
// deviation specified by the "--odometry_stddev" flag, and the range readings
// similarly with standard deviation specified by the "--range-stddev" flag.
//
// There are two types of residuals in this problem:
// 1) The OdometryConstraint residual, that accounts for the odometry readings
//    between successive pose estimatess of the robot.
// 2) The RangeConstraint residual, that accounts for the errors in the observed
//    range readings from each pose.
//
// The OdometryConstraint residual is modeled as an AutoDiffCostFunction with
// a fixed parameter block size of 1, which is the relative odometry being
// solved for, between a pair of successive poses of the robot. Differences
// between observed and computed relative odometry values are penalized weighted
// by the known standard deviation of the odometry readings.
//
// The RangeConstraint residual is modeled as a DynamicAutoDiffCostFunction
// which sums up the relative odometry estimates to compute the estimated
// global pose of the robot, and then computes the expected range reading.
// Differences between the observed and expected range readings are then
// penalized weighted by the standard deviation of readings of the sensor.
// Since the number of poses of the robot is not known at compile time, this
// cost function is implemented as a DynamicAutoDiffCostFunction.
//
// The outputs of the example are the initial values of the odometry and range
// readings, and the range and odometry errors for every pose of the robot.
// After computing the MLE, the computed poses and corrected odometry values
// are printed out, along with the corresponding range and odometry errors. Note
// that as an MLE of a noisy system the errors will not be reduced to zero, but
// the odometry estimates will be updated to maximize the joint likelihood of
// all odometry and range readings of the robot.
//
// Mathematical Formulation
// ======================================================
//
// Let p_0, .., p_N be (N+1) robot poses, where the robot moves down the
// corridor starting from p_0 and ending at p_N. We assume that p_0 is the
// origin of the coordinate system.
// Odometry u_i is the observed relative odometry between pose p_(i-1) and p_i,
// and range reading y_i is the range reading of the end of the corridor from
// pose p_i. Both odometry as well as range readings are noisy, but we wish to
// compute the maximum likelihood estimate (MLE) of corrected odometry values
// u*_0 to u*_(N-1), such that the Belief is optimized:
//
// Belief(u*_(0:N-1) | u_(0:N-1), y_(0:N-1))                                  1.
//   =        P(u*_(0:N-1) | u_(0:N-1), y_(0:N-1))                            2.
//   \propto  P(y_(0:N-1) | u*_(0:N-1), u_(0:N-1)) P(u*_(0:N-1) | u_(0:N-1))  3.
//   =       \prod_i{ P(y_i | u*_(0:i)) P(u*_i | u_i) }                       4.
//
// Here, the subscript "(0:i)" is used as shorthand to indicate entries from all
// timesteps 0 to i for that variable, both inclusive.
//
// Bayes' rule is used to derive eq. 3 from 2, and the independence of
// odometry observations and range readings is expolited to derive 4 from 3.
//
// Thus, the Belief, up to scale, is factored as a product of a number of
// terms, two for each pose, where for each pose term there is one term for the
// range reading, P(y_i | u*_(0:i) and one term for the odometry reading,
// P(u*_i | u_i) . Note that the term for the range reading is dependent on all
// odometry values u*_(0:i), while the odometry term, P(u*_i | u_i) depends only
// on a single value, u_i. Both the range reading as well as odoemtry
// probability terms are modeled as the Normal distribution, and have the form:
//
// p(x) \propto \exp{-((x - x_mean) / x_stddev)^2}
//
// where x refers to either the MLE odometry u* or range reading y, and x_mean
// is the corresponding mean value, u for the odometry terms, and y_expected,
// the expected range reading based on all the previous odometry terms.
// The MLE is thus found by finding those values x* which minimize:
//
// x* = \arg\min{((x - x_mean) / x_stddev)^2}
//
// which is in the nonlinear least-square form, suited to being solved by Ceres.
// The non-linear component arise from the computation of x_mean. The residuals
// ((x - x_mean) / x_stddev) for the residuals that Ceres will optimize. As
// mentioned earlier, the odometry term for each pose depends only on one
// variable, and will be computed by an AutoDiffCostFunction, while the term
// for the range reading will depend on all previous odometry observations, and
// will be computed by a DynamicAutoDiffCostFunction since the number of
// odoemtry observations will only be known at run time.

#include <cstdio>
#include <math.h>
#include <vector>

#include "ceres/ceres.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "random.h"

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::LossFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::examples::RandNormal;
using std::min;
using std::vector;

DEFINE_double(corridor_length, 30.0, "Length of the corridor that the robot is "
              "travelling down.");

DEFINE_double(pose_separation, 0.5, "The distance that the robot traverses "
              "between successive odometry updates.");

DEFINE_double(odometry_stddev, 0.1, "The standard deviation of "
              "odometry error of the robot.");

DEFINE_double(range_stddev, 0.01, "The standard deviation of range readings of "
              "the robot.");

// The stride length of the dynamic_autodiff_cost_function evaluator.
static const int kStride = 10;

struct OdometryConstraint {
  typedef AutoDiffCostFunction<OdometryConstraint, 1, 1> OdometryCostFunction;

  OdometryConstraint(double odometry_mean, double odometry_stddev) :
      odometry_mean(odometry_mean), odometry_stddev(odometry_stddev) {}

  template <typename T>
  bool operator()(const T* const odometry, T* residual) const {
    *residual = (*odometry - odometry_mean) / odometry_stddev;
    return true;
  }

  static OdometryCostFunction* Create(const double odometry_value) {
    return new OdometryCostFunction(
        new OdometryConstraint(odometry_value, FLAGS_odometry_stddev));
  }

  const double odometry_mean;
  const double odometry_stddev;
};

struct RangeConstraint {
  typedef DynamicAutoDiffCostFunction<RangeConstraint, kStride>
      RangeCostFunction;

  RangeConstraint(
      int pose_index,
      double range_reading,
      double range_stddev,
      double corridor_length) :
      pose_index(pose_index), range_reading(range_reading),
      range_stddev(range_stddev), corridor_length(corridor_length) {}

  template <typename T>
  bool operator()(T const* const* relative_poses, T* residuals) const {
    T global_pose(0);
    for (int i = 0; i <= pose_index; ++i) {
      global_pose += relative_poses[i][0];
    }
    residuals[0] = (global_pose + range_reading - corridor_length) /
        range_stddev;
    return true;
  }

  // Factory method to create a CostFunction from a RangeConstraint to
  // conveniently add to a ceres problem.
  static RangeCostFunction* Create(const int pose_index,
                                   const double range_reading,
                                   vector<double>* odometry_values,
                                   vector<double*>* parameter_blocks) {
    RangeConstraint* constraint = new RangeConstraint(
        pose_index, range_reading, FLAGS_range_stddev, FLAGS_corridor_length);
    RangeCostFunction* cost_function = new RangeCostFunction(constraint);
    // Add all the parameter blocks that affect this constraint.
    parameter_blocks->clear();
    for (int i = 0; i <= pose_index; ++i) {
      parameter_blocks->push_back(&((*odometry_values)[i]));
      cost_function->AddParameterBlock(1);
    }
    cost_function->SetNumResiduals(1);
    return (cost_function);
  }

  const int pose_index;
  const double range_reading;
  const double range_stddev;
  const double corridor_length;
};

void SimulateRobot(vector<double>* odometry_values,
                   vector<double>* range_readings) {
  const int num_steps = static_cast<int>(
      ceil(FLAGS_corridor_length / FLAGS_pose_separation));

  // The robot starts out at the origin.
  double robot_location = 0.0;
  for (int i = 0; i < num_steps; ++i) {
    const double actual_odometry_value = min(
        FLAGS_pose_separation, FLAGS_corridor_length - robot_location);
    robot_location += actual_odometry_value;
    const double actual_range = FLAGS_corridor_length - robot_location;
    const double observed_odometry =
        RandNormal() * FLAGS_odometry_stddev + actual_odometry_value;
    const double observed_range =
        RandNormal() * FLAGS_range_stddev + actual_range;
    odometry_values->push_back(observed_odometry);
    range_readings->push_back(observed_range);
  }
}

void PrintState(const vector<double>& odometry_readings,
                const vector<double>& range_readings) {
  CHECK_EQ(odometry_readings.size(), range_readings.size());
  double robot_location = 0.0;
  printf("pose: location     odom    range  r.error  o.error\n");
  for (int i = 0; i < odometry_readings.size(); ++i) {
    robot_location += odometry_readings[i];
    const double range_error =
        robot_location + range_readings[i] - FLAGS_corridor_length;
    const double odometry_error =
        FLAGS_pose_separation - odometry_readings[i];
    printf("%4d: %8.3f %8.3f %8.3f %8.3f %8.3f\n",
           static_cast<int>(i), robot_location, odometry_readings[i],
           range_readings[i], range_error, odometry_error);
  }
}

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  CERES_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  // Make sure that the arguments parsed are all positive.
  CHECK_GT(FLAGS_corridor_length, 0.0);
  CHECK_GT(FLAGS_pose_separation, 0.0);
  CHECK_GT(FLAGS_odometry_stddev, 0.0);
  CHECK_GT(FLAGS_range_stddev, 0.0);

  vector<double> odometry_values;
  vector<double> range_readings;
  SimulateRobot(&odometry_values, &range_readings);

  printf("Initial values:\n");
  PrintState(odometry_values, range_readings);
  ceres::Problem problem;

  for (int i = 0; i < odometry_values.size(); ++i) {
    // Create and add a DynamicAutoDiffCostFunction for the RangeConstraint from
    // pose i.
    vector<double*> parameter_blocks;
    RangeConstraint::RangeCostFunction* range_cost_function =
        RangeConstraint::Create(
            i, range_readings[i], &odometry_values, &parameter_blocks);
    problem.AddResidualBlock(range_cost_function, NULL, parameter_blocks);

    // Create and add an AutoDiffCostFunction for the OdometryConstraint for
    // pose i.
    problem.AddResidualBlock(OdometryConstraint::Create(odometry_values[i]),
                             NULL,
                             &(odometry_values[i]));
  }

  ceres::Solver::Options solver_options;
  solver_options.minimizer_progress_to_stdout = true;

  Solver::Summary summary;
  printf("Solving...\n");
  Solve(solver_options, &problem, &summary);
  printf("Done.\n");
  std::cout << summary.FullReport() << "\n";
  printf("Final values:\n");
  PrintState(odometry_values, range_readings);
  return 0;
}
