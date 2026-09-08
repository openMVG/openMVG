// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/parameter_block_ordering.h"

#include <cstddef>
#include <memory>
#include <unordered_set>
#include <vector>

#include "ceres/cost_function.h"
#include "ceres/graph.h"
#include "ceres/problem_impl.h"
#include "ceres/program.h"
#include "ceres/sized_cost_function.h"
#include "ceres/stl_util.h"
#include "gtest/gtest.h"

namespace ceres::internal {

using VertexSet = std::unordered_set<ParameterBlock*>;

template <int M, int... Ns>
class DummyCostFunction : public SizedCostFunction<M, Ns...> {
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    return true;
  }
};

class SchurOrderingTest : public ::testing::Test {
 protected:
  void SetUp() final {
    // The explicit calls to AddParameterBlock are necessary because
    // the below tests depend on the specific numbering of the
    // parameter blocks.
    problem_.AddParameterBlock(x_, 3);
    problem_.AddParameterBlock(y_, 4);
    problem_.AddParameterBlock(z_, 5);
    problem_.AddParameterBlock(w_, 6);

    problem_.AddResidualBlock(new DummyCostFunction<2, 3>, nullptr, x_);
    problem_.AddResidualBlock(new DummyCostFunction<6, 5, 4>, nullptr, z_, y_);
    problem_.AddResidualBlock(new DummyCostFunction<3, 3, 5>, nullptr, x_, z_);
    problem_.AddResidualBlock(new DummyCostFunction<7, 5, 3>, nullptr, z_, x_);
    problem_.AddResidualBlock(
        new DummyCostFunction<1, 5, 3, 6>, nullptr, z_, x_, w_);
  }

  ProblemImpl problem_;
  double x_[3], y_[4], z_[5], w_[6];
};

TEST_F(SchurOrderingTest, NoFixed) {
  const Program& program = problem_.program();
  const std::vector<ParameterBlock*>& parameter_blocks =
      program.parameter_blocks();
  auto graph = CreateHessianGraph(program);

  const VertexSet& vertices = graph->vertices();
  EXPECT_EQ(vertices.size(), 4);

  for (int i = 0; i < 4; ++i) {
    EXPECT_TRUE(vertices.find(parameter_blocks[i]) != vertices.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[0]);
    EXPECT_EQ(neighbors.size(), 2);
    EXPECT_TRUE(neighbors.find(parameter_blocks[2]) != neighbors.end());
    EXPECT_TRUE(neighbors.find(parameter_blocks[3]) != neighbors.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[1]);
    EXPECT_EQ(neighbors.size(), 1);
    EXPECT_TRUE(neighbors.find(parameter_blocks[2]) != neighbors.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[2]);
    EXPECT_EQ(neighbors.size(), 3);
    EXPECT_TRUE(neighbors.find(parameter_blocks[0]) != neighbors.end());
    EXPECT_TRUE(neighbors.find(parameter_blocks[1]) != neighbors.end());
    EXPECT_TRUE(neighbors.find(parameter_blocks[3]) != neighbors.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[3]);
    EXPECT_EQ(neighbors.size(), 2);
    EXPECT_TRUE(neighbors.find(parameter_blocks[0]) != neighbors.end());
    EXPECT_TRUE(neighbors.find(parameter_blocks[2]) != neighbors.end());
  }
}

TEST_F(SchurOrderingTest, AllFixed) {
  problem_.SetParameterBlockConstant(x_);
  problem_.SetParameterBlockConstant(y_);
  problem_.SetParameterBlockConstant(z_);
  problem_.SetParameterBlockConstant(w_);

  const Program& program = problem_.program();
  auto graph = CreateHessianGraph(program);
  EXPECT_EQ(graph->vertices().size(), 0);
}

TEST_F(SchurOrderingTest, OneFixed) {
  problem_.SetParameterBlockConstant(x_);

  const Program& program = problem_.program();
  const std::vector<ParameterBlock*>& parameter_blocks =
      program.parameter_blocks();
  auto graph = CreateHessianGraph(program);

  const VertexSet& vertices = graph->vertices();

  EXPECT_EQ(vertices.size(), 3);
  EXPECT_TRUE(vertices.find(parameter_blocks[0]) == vertices.end());

  for (int i = 1; i < 3; ++i) {
    EXPECT_TRUE(vertices.find(parameter_blocks[i]) != vertices.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[1]);
    EXPECT_EQ(neighbors.size(), 1);
    EXPECT_TRUE(neighbors.find(parameter_blocks[2]) != neighbors.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[2]);
    EXPECT_EQ(neighbors.size(), 2);
    EXPECT_TRUE(neighbors.find(parameter_blocks[1]) != neighbors.end());
    EXPECT_TRUE(neighbors.find(parameter_blocks[3]) != neighbors.end());
  }

  {
    const VertexSet& neighbors = graph->Neighbors(parameter_blocks[3]);
    EXPECT_EQ(neighbors.size(), 1);
    EXPECT_TRUE(neighbors.find(parameter_blocks[2]) != neighbors.end());
  }

  // The constant parameter block is at the end.
  std::vector<ParameterBlock*> ordering;
  ComputeSchurOrdering(program, &ordering);
  EXPECT_EQ(ordering.back(), parameter_blocks[0]);
}

}  // namespace ceres::internal
